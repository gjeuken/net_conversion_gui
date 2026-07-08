"""KEGG-compound -> BiGG-id translation for human-readable working ids.

This is a *convenience* layer, never ground truth (CLAUDE.md keeps the KEGG id
load-bearing for balance + eQuilibrator).  It:

1. downloads BiGG's official ``bigg_models_metabolites.txt`` once and builds a
   ``{KEGG compound id: [BiGG ids]}`` map, cached on disk like the KEGG cache;
2. picks a deterministic BiGG id per KEGG compound (shortest, then alphabetical
   — the base form ``h2o`` beats the variant ``oh1``);
3. rewrites the working ids of a metabolite/reaction pair to those BiGG ids,
   **keeping the ``KEGG ID`` column intact** so thermodynamics and balance still
   resolve compounds by KEGG.

Reliability is partial (measured ~87% clean 1:1, ~6% ambiguous, ~8% no match on
central-carbon metabolism).  Compounds with no clean BiGG id fall back to a
slugified KEGG name and are flagged in the report so the user can hand-edit.
"""

from __future__ import annotations

import re
import urllib.request

from . import cache
from .idmap import apply_id_map
from .kegg import get_compound_name

BIGG_URL = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"
_CACHE_KEY = "bigg:kegg2bigg"
TIMEOUT = 90

_MAP = None  # in-process cache of {kegg_id: [bigg_id, ...]}


class BiggError(RuntimeError):
    """Raised when the BiGG namespace file can't be fetched."""


# --------------------------------------------------------------------------- #
# The KEGG -> BiGG map
# --------------------------------------------------------------------------- #
def get_kegg_to_bigg(force_refresh: bool = False):
    """Return ``{kegg_id: [bigg_id, ...]}``, building + caching it on first use.

    Raises :class:`BiggError` if the file must be downloaded but can't be
    reached (offline) and nothing is cached yet.
    """
    global _MAP
    if _MAP is not None and not force_refresh:
        return _MAP
    if not force_refresh:
        cached = cache.get(_CACHE_KEY)
        if cached:
            _MAP = cached
            return _MAP
    _MAP = _download_and_parse()
    cache.put(_CACHE_KEY, _MAP)
    return _MAP


def _download_and_parse():
    try:
        raw = urllib.request.urlopen(BIGG_URL, timeout=TIMEOUT).read().decode(
            "utf-8", "replace")
    except Exception as e:  # network, HTTP, timeout, ...
        raise BiggError(f"Could not download the BiGG namespace file: {e}") from e

    lines = raw.splitlines()
    if not lines:
        raise BiggError("BiGG namespace file was empty")
    header = lines[0].split("\t")
    try:
        li = header.index("database_links")
        bi = header.index("universal_bigg_id")
    except ValueError as e:
        raise BiggError(f"Unexpected BiGG file format: {e}") from e

    mapping: dict[str, list[str]] = {}
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) <= li:
            continue
        bigg = parts[bi].strip()
        if not bigg:
            continue
        for cid in re.findall(r"kegg\.compound/(C\d{5})", parts[li]):
            bag = mapping.setdefault(cid, [])
            if bigg not in bag:
                bag.append(bigg)
    return mapping


def choose_bigg(kegg_id, mapping):
    """Deterministic pick: shortest BiGG id, ties broken alphabetically."""
    ids = mapping.get(kegg_id)
    if not ids:
        return None
    return sorted(ids, key=lambda s: (len(s), s))[0]


# --------------------------------------------------------------------------- #
# Slugify (fallback when there's no clean BiGG id)
# --------------------------------------------------------------------------- #
def slugify(name, fallback):
    """Turn a compound name into a compact, readable, id-safe token."""
    if not name or not str(name).strip():
        return str(fallback)
    s = str(name).strip()
    s = s.split(";")[0]                      # first synonym only
    s = re.sub(r"[^\w\s-]", "", s)           # drop punctuation (keep - and _)
    s = re.sub(r"\s+", "_", s.strip())       # spaces -> underscore
    s = s.strip("_-")
    return s or str(fallback)


# --------------------------------------------------------------------------- #
# Translate a metabolite / reaction pair
# --------------------------------------------------------------------------- #
def translate_to_bigg(df_metabolites, df_reactions, mapping=None):
    """Rewrite working ids to BiGG ids (KEGG id kept in the ``KEGG ID`` column).

    Returns ``(df_met, df_rxn, report)``.  ``report`` is a list of dicts, one per
    metabolite: ``{old, new, kegg, status}`` with ``status`` in
    ``{"bigg", "bigg-ambiguous", "name-fallback", "kegg-fallback", "collision"}``.
    Reaction stoichiometry strings and ``EX{id}`` reaction ids are rewritten to
    match; the reaction that carries no mapped id is left untouched.
    """
    mapping = get_kegg_to_bigg() if mapping is None else mapping
    df_met = df_metabolites.copy()
    df_rxn = df_reactions.copy()

    id_map, report = {}, []
    used = set()
    for _, row in df_met.iterrows():
        old = str(row.get("ID", "")).strip()
        kegg = str(row.get("KEGG ID", "")).strip()
        if not old:
            continue
        kegg_valid = bool(re.fullmatch(r"C\d{5}", kegg))
        biggs = mapping.get(kegg, []) if kegg_valid else []
        if biggs:
            new = choose_bigg(kegg, mapping)
            status = "bigg" if len(biggs) == 1 else "bigg-ambiguous"
        elif kegg_valid:
            new = slugify(get_compound_name(kegg), kegg)
            status = "name-fallback"
        else:
            new = old
            status = "kegg-fallback"

        # Keep ids unique: if the chosen id collides with one already assigned,
        # suffix it and flag the row.
        base, n = new, 2
        while new in used and new != old:
            new = f"{base}_{n}"
            n += 1
            status = "collision"
        used.add(new)
        if new != old:
            id_map[old] = new
        report.append({"old": old, "new": new, "kegg": kegg, "status": status})

    df_met, df_rxn = apply_id_map(df_met, df_rxn, id_map)
    return df_met, df_rxn, report
