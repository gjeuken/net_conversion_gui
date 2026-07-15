"""Consistently rename metabolite working ids across the two canonical frames.

Renaming a metabolite's working ``ID`` must also rewrite every reaction
stoichiometry string that mentions it and any ``EX{id}`` exchange-reaction id
that encodes it — otherwise the reaction/metabolite links break.  Both the BiGG
translation (:mod:`pipeline.bigg`) and user-supplied custom ids go through
:func:`apply_id_map`.
"""

from __future__ import annotations

import re

import pandas as pd


def apply_id_map(df_metabolites, df_reactions, id_map):
    """Return ``(df_met, df_rxn)`` with ``id_map`` (``{old: new}``) applied.

    Rewrites the metabolite ``ID`` column, the reaction stoichiometry strings,
    and ``EX{old}`` exchange-reaction ids.  A no-op when ``id_map`` is empty.

    If the rewrite makes two rows share the same final id (e.g. two rows for
    the same compound — one already translated, one still on its raw KEGG id),
    they are merged into a single row rather than left as duplicates.
    """
    df_met = df_metabolites.copy()
    df_rxn = df_reactions.copy()
    if not id_map:
        return df_met, df_rxn
    df_met["ID"] = df_met["ID"].map(lambda x: id_map.get(str(x).strip(), x))
    df_rxn["Reaction stoichiometry"] = df_rxn["Reaction stoichiometry"].map(
        lambda eq: rewrite_equation(eq, id_map))
    df_rxn["ID"] = df_rxn["ID"].map(lambda r: rewrite_exchange_id(r, id_map))
    df_met = merge_duplicate_ids(df_met)
    df_rxn = drop_duplicate_ids(df_rxn)
    return df_met, df_rxn


def _is_empty(v) -> bool:
    return v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() == ""


def merge_duplicate_ids(df, id_col="ID"):
    """Merge rows that share the same ``id_col`` value.

    Keeps the first occurrence's position; for every other column, a later
    row's value fills in only if the first row's value there is empty.  Rows
    with a blank id are left untouched (not merged into each other).
    """
    ids = df[id_col].astype(str).str.strip()
    if not (ids[ids != ""].duplicated()).any():
        return df

    merged: dict[str, dict] = {}
    order: list = []
    for i, row in df.iterrows():
        mid = str(row.get(id_col, "")).strip()
        if not mid:
            order.append(("_row", i))
            continue
        if mid not in merged:
            merged[mid] = row.to_dict()
            order.append(("_id", mid))
        else:
            existing = merged[mid]
            for col in df.columns:
                if col == id_col:
                    continue
                if _is_empty(existing.get(col)) and not _is_empty(row.get(col)):
                    existing[col] = row.get(col)

    out_rows = []
    for kind, key in order:
        out_rows.append(df.loc[key].to_dict() if kind == "_row" else merged[key])
    return pd.DataFrame(out_rows, columns=df.columns)


def drop_duplicate_ids(df, id_col="ID"):
    """Drop later rows that share an already-seen non-blank ``id_col`` value."""
    ids = df[id_col].astype(str).str.strip()
    if not (ids[ids != ""].duplicated()).any():
        return df
    keep = ~(ids.duplicated() & (ids != ""))
    return df[keep].reset_index(drop=True)


def rewrite_equation(eq, id_map):
    """Replace whole-token metabolite ids in a stoichiometry string."""
    if not isinstance(eq, str) or not id_map:
        return eq
    pattern = re.compile(
        r"(?<![A-Za-z0-9_])(" +
        "|".join(re.escape(k) for k in sorted(id_map, key=len, reverse=True)) +
        r")(?![A-Za-z0-9_])")
    return pattern.sub(lambda m: id_map[m.group(1)], eq)


def rewrite_exchange_id(rid, id_map):
    """Rewrite ``EX{old}`` -> ``EX{new}`` when the metabolite id is remapped."""
    rid = str(rid)
    if rid.startswith("EX") and rid[2:] in id_map:
        return "EX" + id_map[rid[2:]]
    return rid
