"""Auto-generate exchange + transport reactions for chosen boundary metabolites.

KEGG has no compartments, so a modelled pathway needs, for each metabolite that
crosses the system boundary, three things the paper's examples spell out (see
any ``Example*.xlsx``):

* an **extracellular** counterpart metabolite ``{X}ex`` (same KEGG id/formula),
* a **transport** reaction ``T{X}`` : ``X <=> {X}ex``,
* an **exchange** reaction ``EX{X}`` : ``{X}ex <=>``  (the ``EX`` prefix is
  load-bearing — the thermo step strips it to map back to KEGG).

This module builds those rows in the two canonical dataframes.  Direction is a
modelling choice made later (in the analysis app's SUBSTRATES/PRODUCTS/
REV_ALLOWED selectors), so everything defaults to reversible here.
"""

from __future__ import annotations

import pandas as pd

from .io import METABOLITE_COLS, REACTION_COLS


def add_exchange_transport(df_metabolites, df_reactions, metabolite_ids,
                           reversible=True):
    """Add ``{X}ex`` metabolite + ``T{X}`` transport + ``EX{X}`` exchange rows.

    Parameters
    ----------
    metabolite_ids : iterable of str
        Working ids of the (intracellular) boundary metabolites.
    reversible : bool
        Reversibility for the created transport/exchange rows (default True).

    Returns ``(df_met, df_rxn, messages)``.  Existing rows are never duplicated.
    """
    df_met = df_metabolites.copy()
    df_rxn = df_reactions.copy()
    for c in METABOLITE_COLS:
        if c not in df_met.columns:
            df_met[c] = pd.Series(dtype=object)
    for c in REACTION_COLS:
        if c not in df_rxn.columns:
            df_rxn[c] = pd.Series(dtype=object)

    met_by_id = {str(r["ID"]).strip(): r for _, r in df_met.iterrows()}
    known_met = set(met_by_id)
    known_rxn = {str(x).strip() for x in df_rxn["ID"].dropna()}
    rev = 1 if reversible else 0

    new_met, new_rxn, messages = [], [], []
    for raw in metabolite_ids or []:
        base = str(raw).strip()
        if not base:
            continue
        if base.endswith("ex"):
            messages.append(f"• {base} looks extracellular already — skipped")
            continue
        if base not in known_met:
            messages.append(f"⚠ {base} not in the metabolite table — skipped")
            continue

        ex_id = f"{base}ex"
        t_id = f"T{base}"
        x_id = f"EX{base}"
        src = met_by_id[base]

        if ex_id not in known_met:
            name = str(src.get("Name") or base)
            new_met.append({
                "ID": ex_id,
                "Name": f"{name}, extracellular",
                "KEGG ID": src.get("KEGG ID"),
                "Chemical formula": src.get("Chemical formula"),
            })
            known_met.add(ex_id)

        added = []
        if t_id not in known_rxn:
            new_rxn.append({"ID": t_id, "Name": f"{base} transport",
                            "Reaction stoichiometry": f"{base} <=> {ex_id}",
                            "Reversibility": rev})
            known_rxn.add(t_id)
            added.append(t_id)
        if x_id not in known_rxn:
            new_rxn.append({"ID": x_id, "Name": f"{base} exchange",
                            "Reaction stoichiometry": f"{ex_id} <=>",
                            "Reversibility": rev})
            known_rxn.add(x_id)
            added.append(x_id)

        if added:
            messages.append(f"✓ {base}: added {', '.join(added)}")
        else:
            messages.append(f"• {base}: exchange + transport already present")

    if new_met:
        df_met = pd.concat([df_met, pd.DataFrame(new_met)], ignore_index=True)
    if new_rxn:
        df_rxn = pd.concat([df_rxn, pd.DataFrame(new_rxn)], ignore_index=True)
    return df_met, df_rxn, messages
