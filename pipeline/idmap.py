"""Consistently rename metabolite working ids across the two canonical frames.

Renaming a metabolite's working ``ID`` must also rewrite every reaction
stoichiometry string that mentions it and any ``EX{id}`` exchange-reaction id
that encodes it — otherwise the reaction/metabolite links break.  Both the BiGG
translation (:mod:`pipeline.bigg`) and user-supplied custom ids go through
:func:`apply_id_map`.
"""

from __future__ import annotations

import re


def apply_id_map(df_metabolites, df_reactions, id_map):
    """Return ``(df_met, df_rxn)`` with ``id_map`` (``{old: new}``) applied.

    Rewrites the metabolite ``ID`` column, the reaction stoichiometry strings,
    and ``EX{old}`` exchange-reaction ids.  A no-op when ``id_map`` is empty.
    """
    df_met = df_metabolites.copy()
    df_rxn = df_reactions.copy()
    if not id_map:
        return df_met, df_rxn
    df_met["ID"] = df_met["ID"].map(lambda x: id_map.get(str(x).strip(), x))
    df_rxn["Reaction stoichiometry"] = df_rxn["Reaction stoichiometry"].map(
        lambda eq: rewrite_equation(eq, id_map))
    df_rxn["ID"] = df_rxn["ID"].map(lambda r: rewrite_exchange_id(r, id_map))
    return df_met, df_rxn


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
