"""KEGG REST access: compound formula/charge/name and reaction fetching.

The notebook only fetched *compounds* (``get_compound_info``).  This module adds
the *reaction* fetcher required by the workbench (CLAUDE.md §5): given KEGG
reaction ids (or a module id), it returns rows for the two canonical
dataframes, handling KEGG's quirks explicitly (R-groups/polymers, missing
H+/H2O, unbalanced-as-written) by surfacing them rather than guessing.

Both compound and reaction entries are cached on disk via :mod:`pipeline.cache`.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional

import pandas as pd
import requests

from . import cache
from .balance import has_nonnumeric_coefficient, normalize_arrow, parse_equation

KEGG_BASE = "http://rest.kegg.jp/get/"
TIMEOUT = 20


class KeggError(RuntimeError):
    """Raised when a KEGG entry cannot be fetched or parsed."""


def _fetch_text(entry: str) -> str:
    """GET a KEGG flat-file entry, raising :class:`KeggError` on failure."""
    url = f"{KEGG_BASE}{entry}"
    try:
        r = requests.get(url, timeout=TIMEOUT)
        r.raise_for_status()
    except requests.RequestException as e:
        raise KeggError(f"KEGG request failed for {entry}: {e}") from e
    if not r.text.strip():
        raise KeggError(f"KEGG returned an empty entry for {entry}")
    return r.text


# -------- Compounds --------
@lru_cache(maxsize=10000)
def get_compound_info(cid: str):
    """Return ``(formula_str | None, charge_int)`` for a Cxxxxx id.

    Cached in-process (lru) and on disk.  On network error the compound is
    treated as unknown (``(None, 0)``) so balance-checking degrades gracefully
    offline rather than crashing the whole analysis.
    """
    key = f"cpd:{cid}"
    cached = cache.get(key)
    if cached is not None:
        return cached[0], cached[1]

    try:
        text = _fetch_text(cid)
    except KeggError:
        return None, 0

    formula, charge = None, 0
    for line in text.splitlines():
        m = re.match(r"^FORMULA\s+(.+)$", line)
        if m:
            formula = m.group(1).strip()
            continue
        m = re.match(r"^CHARGE\s+(-?\d+)$", line)
        if m:
            try:
                charge = int(m.group(1))
            except ValueError:
                pass
    cache.put(key, [formula, charge])
    return formula, charge


@lru_cache(maxsize=10000)
def get_compound_name(cid: str) -> Optional[str]:
    """Return the first NAME of a Cxxxxx compound (trimmed of trailing ';')."""
    key = f"cpdname:{cid}"
    cached = cache.get(key)
    if cached is not None:
        return cached
    try:
        text = _fetch_text(cid)
    except KeggError:
        return None
    name = None
    for line in text.splitlines():
        m = re.match(r"^NAME\s+(.+)$", line)
        if m:
            name = m.group(1).strip().rstrip(";").strip()
            break
    cache.put(key, name)
    return name


# -------- Reactions --------
def get_reaction_entry(rid: str):
    """Return ``(name, equation)`` for a Rxxxxx reaction id."""
    key = f"rn:{rid}"
    cached = cache.get(key)
    if cached is not None:
        return cached[0], cached[1]

    text = _fetch_text(rid)
    name, equation = None, None
    for line in text.splitlines():
        m = re.match(r"^NAME\s+(.+)$", line)
        if m and name is None:
            name = m.group(1).strip().rstrip(";").strip()
            continue
        m = re.match(r"^EQUATION\s+(.+)$", line)
        if m:
            equation = m.group(1).strip()
    if equation is None:
        raise KeggError(f"No EQUATION line in KEGG entry {rid}")
    cache.put(key, [name, equation])
    return name, equation


def get_module_reactions(mid: str):
    """Return the ordered list of reaction ids referenced by a KEGG module."""
    text = _fetch_text(mid)
    rxns: list[str] = []
    in_section = False
    for line in text.splitlines():
        if line.startswith("REACTION"):
            in_section = True
            payload = line[len("REACTION"):]
        elif in_section and line.startswith(" "):
            payload = line
        else:
            in_section = False
            continue
        for rid in re.findall(r"R\d{5}", payload):
            if rid not in rxns:
                rxns.append(rid)
    if not rxns:
        raise KeggError(f"No REACTION block found in module {mid}")
    return rxns


def fetch_reactions(ids, existing_metabolites=None, existing_reactions=None):
    """Fetch KEGG reaction (or module) ids into the two canonical dataframes.

    Parameters
    ----------
    ids : iterable of str
        KEGG reaction ids (``Rxxxxx``) and/or module ids (``Mxxxxx``).
    existing_metabolites, existing_reactions : pandas.DataFrame, optional
        Current workbench tables to append to (avoids duplicate ids).

    Returns
    -------
    (df_metabolites, df_reactions, messages)
        The two canonical dataframes (schema per CLAUDE.md §5) and a list of
        human-readable status/warning strings.

    The KEGG ID *is* the working compound ID here, so no ID↔KEGG mapping step is
    needed.  Reactions with R-group/polymer coefficients are flagged and skipped
    (routed to manual entry).  Atom/charge consistency is left to the live
    balance loop.
    """
    met_cols = ["ID", "Name", "KEGG ID"]
    rxn_cols = ["ID", "Name", "Reaction stoichiometry", "Reversibility"]

    df_met = (existing_metabolites.copy() if existing_metabolites is not None
              else pd.DataFrame(columns=met_cols))
    df_rxn = (existing_reactions.copy() if existing_reactions is not None
              else pd.DataFrame(columns=rxn_cols))
    for c in met_cols:
        if c not in df_met.columns:
            df_met[c] = pd.Series(dtype=object)
    for c in rxn_cols:
        if c not in df_rxn.columns:
            df_rxn[c] = pd.Series(dtype=object)

    messages: list[str] = []
    known_met = set(df_met["ID"].dropna().astype(str))
    known_rxn = set(df_rxn["ID"].dropna().astype(str))

    # Expand any module ids into their reaction lists.
    expanded: list[str] = []
    for raw in ids:
        rid = str(raw).strip().upper()
        if not rid:
            continue
        if rid.startswith("M"):
            try:
                module_rxns = get_module_reactions(rid)
                messages.append(f"Module {rid}: {len(module_rxns)} reactions")
                expanded.extend(module_rxns)
            except KeggError as e:
                messages.append(f"⚠ {e}")
        else:
            expanded.append(rid)

    new_met_rows, new_rxn_rows = [], []

    def ensure_metabolite(cid):
        if cid in known_met:
            return
        name = get_compound_name(cid) or cid
        new_met_rows.append({"ID": cid, "Name": name, "KEGG ID": cid})
        known_met.add(cid)

    for rid in expanded:
        if rid in known_rxn:
            messages.append(f"• {rid} already present — skipped")
            continue
        try:
            name, equation = get_reaction_entry(rid)
        except KeggError as e:
            messages.append(f"⚠ {e}")
            continue

        if has_nonnumeric_coefficient(equation):
            messages.append(
                f"⚠ {rid}: R-group / polymer coefficient — enter manually "
                f"({equation})")
            continue

        try:
            subs, prods = parse_equation(equation)
        except Exception as e:
            messages.append(f"⚠ {rid}: could not parse equation ({e})")
            continue

        for _coeff, cid in subs + prods:
            ensure_metabolite(cid)

        # KEGG writes everything reversible; reversibility/direction is a
        # user decision (CLAUDE.md §5) — default to reversible (1).
        new_rxn_rows.append({
            "ID": rid,
            "Name": name or rid,
            "Reaction stoichiometry": normalize_arrow(equation),
            "Reversibility": 1,
        })
        known_rxn.add(rid)
        messages.append(f"✓ {rid} added")

    if new_met_rows:
        df_met = pd.concat([df_met, pd.DataFrame(new_met_rows)], ignore_index=True)
    if new_rxn_rows:
        df_rxn = pd.concat([df_rxn, pd.DataFrame(new_rxn_rows)], ignore_index=True)

    return df_met, df_rxn, messages
