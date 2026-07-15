"""COBRApy model building, exchange configuration, FBA sanity checks, pruning.

Ported from the notebook's model cells.  Pure functions: dataframes / config in,
``cobra.Model`` and structured results out (no prints, no globals).
"""

from __future__ import annotations

import math

import cobra
import pandas as pd


def build_model(df_metabolites, df_reactions, model_name="model"):
    """Create a COBRApy model from the two canonical dataframes.

    Irreversible reactions (``Reversibility == 0``) get ``lower_bound = 0``.
    Returns ``(model, irreversible_ids)``.
    """
    model = cobra.Model(name=model_name)
    model.id = model_name

    mets = [
        cobra.Metabolite(
            id=str(df_metabolites.loc[i, "ID"]),
            name=str(df_metabolites.loc[i, "Name"]),
            compartment="c",
        )
        for i in df_metabolites.index
    ]
    model.add_metabolites(mets)

    irreversible = []
    for i in df_reactions.index:
        r = cobra.Reaction(
            id=str(df_reactions.loc[i, "ID"]),
            name=str(df_reactions.loc[i, "Name"]),
            lower_bound=-1000,
            upper_bound=1000,
        )
        model.add_reactions([r])
        r.reaction = str(df_reactions.loc[i, "Reaction stoichiometry"])
        if int(df_reactions.loc[i, "Reversibility"]) == 0:
            r.lower_bound = 0
            irreversible.append(r.id)

    return model, irreversible


def configure_exchanges(model, substrates, products, rev_allowed):
    """Restrict/remove exchange reactions per the substrate/product/rev sets.

    Exchanges not in ``rev_allowed``/``substrates`` lose uptake (lb=0); those in
    none of the three sets are removed entirely.  Returns the list of removed
    reaction ids.
    """
    substrates = set(substrates or [])
    products = set(products or [])
    rev_allowed = set(rev_allowed or [])
    removed = []
    for r in list(model.reactions):
        if not r.id.startswith("EX"):
            continue
        if r.id not in rev_allowed and r.id not in substrates:
            r.lower_bound = 0
        if (r.id not in products and r.id not in substrates
                and r.id not in rev_allowed):
            model.remove_reactions([r])
            removed.append(r.id)
    return removed


def check_atp_without_substrate(model, energy_product, substrates):
    """Maximise the energy product with all substrate uptake blocked.

    A nonzero optimum implies an energy-generating cycle.  Returns
    ``(objective_value, {reaction_id: flux})`` for the offending fluxes.
    """
    with model:
        model.objective = model.reactions.get_by_id(energy_product)
        for s in substrates:
            model.reactions.get_by_id(s).lower_bound = 0
        sol = model.optimize()
        obj = sol.objective_value if sol.status == "optimal" else float("nan")
        fluxes = {}
        if sol.status == "optimal":
            for rid, val in sol.fluxes.items():
                if abs(val) > 1e-6:
                    fluxes[rid] = float(val)
    return obj, fluxes


def max_product_yields(model, substrate, carbon_products, uptake=10.0):
    """Max FBA yield of each carbon product per ``uptake`` of substrate."""
    yields = {}
    with model:
        model.reactions.get_by_id(substrate).bounds = (-uptake, 1000)
        for p in carbon_products:
            model.objective = model.reactions.get_by_id(p)
            sol = model.optimize()
            yields[p] = (float(sol.objective_value) / uptake
                         if sol.status == "optimal" and sol.objective_value is not None
                         else None)
    return yields


def sanitize_model(model):
    """Replace NaN/None bounds, charges, formula weights, names with defaults."""
    for r in model.reactions:
        if r.lower_bound is None or (isinstance(r.lower_bound, float) and math.isnan(r.lower_bound)):
            r.lower_bound = 0.0
        if r.upper_bound is None or (isinstance(r.upper_bound, float) and math.isnan(r.upper_bound)):
            r.upper_bound = 0.0
        if not isinstance(r.name, str) or r.name is None:
            r.name = r.id
    for m in model.metabolites:
        if m.charge is None or (isinstance(m.charge, float) and math.isnan(m.charge)):
            m.charge = 0
    return model


def prune_blocked(model):
    """Remove reactions that cannot carry flux in any steady state.

    Returns the list of removed reaction ids.
    """
    blocked = cobra.flux_analysis.find_blocked_reactions(model)
    for rid in blocked:
        model.remove_reactions([model.reactions.get_by_id(rid)])
    return blocked


def find_blocked_reactions(df_metabolites, df_reactions, model_name="model"):
    """Ids of reactions that can never carry flux, with default bounds only.

    Builds a model straight from the two canonical dataframes — every
    reaction reversible except those marked ``Reversibility == 0`` — with no
    substrate/product/energy-product configuration (those run-time selectors
    are set later, in the analysis app). A dead end found here is usually a
    missing exchange, transport, or connecting reaction: useful right after
    adding exchange/transport reactions, before the network ever reaches FBA.
    """
    model, _irrev = build_model(df_metabolites, df_reactions, model_name)
    sanitize_model(model)
    return cobra.flux_analysis.find_blocked_reactions(model)
