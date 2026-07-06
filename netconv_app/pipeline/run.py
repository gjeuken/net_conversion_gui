"""End-to-end orchestration used by both the CLI and the Dash app.

Each stage is a thin wrapper that returns structured data so callers can render
or inspect intermediate results.  Nothing here prints.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from . import balance, efm, model as model_mod, thermo


@dataclass
class PipelineResult:
    config: dict
    df_balance: Optional[pd.DataFrame] = None
    balance_ok: bool = False
    balance_errors: list = field(default_factory=list)
    removed_exchanges: list = field(default_factory=list)
    atp_without_substrate: Optional[float] = None
    atp_cycle_fluxes: dict = field(default_factory=dict)
    product_yields: dict = field(default_factory=dict)
    blocked_reactions: list = field(default_factory=list)
    n_efms: int = 0
    df_efms: Optional[pd.DataFrame] = None
    df_normalized: Optional[pd.DataFrame] = None
    rank_failures: dict = field(default_factory=dict)
    per_efm_thermo: Optional[pd.DataFrame] = None
    net_conversions: Optional[pd.DataFrame] = None
    notes: list = field(default_factory=list)


def run_balance(df_metabolites, df_reactions):
    """Stage 1: per-reaction atom/charge balance check.

    Returns ``(df_balance, ok, errors)`` where ``errors`` are non-exchange
    reactions that fail (exchanges are expected to be unbalanced).
    """
    df_bal = balance.analyze_custom_equations(
        list(df_reactions["Reaction stoichiometry"]),
        df_metabolites,
        reaction_ids=list(df_reactions["ID"]),
    )
    errors = []
    for _, row in df_bal.iterrows():
        if row.get("is_exchange"):
            continue
        if row.get("atom_balanced") is False or row.get("charge_balanced") is False \
                or (row.get("atom_balanced") is None and "Empty" not in str(row.get("notes"))):
            errors.append({"reaction_id": row.get("reaction_id"),
                           "imbalance": row.get("atom_imbalance"),
                           "notes": row.get("notes")})
    return df_bal, (len(errors) == 0), errors


def build_and_check(df_metabolites, df_reactions, config):
    """Stages 3-7: build model, configure exchanges, FBA checks, prune, S/rev.

    Returns ``(model, result_dict)``.
    """
    model, _irrev = model_mod.build_model(
        df_metabolites, df_reactions, config.get("MODEL_NAME", "model"))
    removed = model_mod.configure_exchanges(
        model, config.get("SUBSTRATES"), config.get("PRODUCTS"),
        config.get("REV_ALLOWED"))

    atp_obj, atp_fluxes = model_mod.check_atp_without_substrate(
        model, config["ENERGY_PRODUCT"], config.get("SUBSTRATES", []))
    yields = model_mod.max_product_yields(
        model, config["SUBSTRATES"][0], config.get("CARBON_PRODUCTS", []))

    model_mod.sanitize_model(model)
    blocked = model_mod.prune_blocked(model)

    return model, {
        "removed_exchanges": removed,
        "atp_without_substrate": atp_obj,
        "atp_cycle_fluxes": atp_fluxes,
        "product_yields": yields,
        "blocked_reactions": blocked,
    }


def enumerate_efms(model, config):
    """Stages 7-10: S matrix, efmtool, rank check, normalisation.

    Returns ``(df_efms, df_normalized, rank_failures)``.
    """
    S, rev, rxn_ids, met_ids = efm.build_stoichiometry(
        model, config.get("REV_ALLOWED"), config.get("SUBSTRATES"))
    df_efms = efm.calculate_efms(S, rev, rxn_ids, met_ids)
    rank_failures = efm.validate_efms(df_efms, S)
    df_norm = efm.normalize_efms(df_efms, config["SUBSTRATES"][0])
    counts, _excluded = efm.count_metabolic_reactions(
        df_norm, model, config.get("PSEUDO_METS"))
    df_norm.loc["n_metabolic_reactions"] = pd.Series(
        {c: counts.get(c) for c in df_norm.columns})
    return df_efms, df_norm, rank_failures


def run_all(df_metabolites, df_reactions, config, with_thermo=True, which="dGm",
            cc=None):
    """Run every stage and return a :class:`PipelineResult`."""
    res = PipelineResult(config=config)

    res.df_balance, res.balance_ok, res.balance_errors = run_balance(
        df_metabolites, df_reactions)

    model, checks = build_and_check(df_metabolites, df_reactions, config)
    res.removed_exchanges = checks["removed_exchanges"]
    res.atp_without_substrate = checks["atp_without_substrate"]
    res.atp_cycle_fluxes = checks["atp_cycle_fluxes"]
    res.product_yields = checks["product_yields"]
    res.blocked_reactions = checks["blocked_reactions"]

    df_efms, df_norm, rank_failures = enumerate_efms(model, config)
    res.df_efms = df_efms
    res.df_normalized = df_norm
    res.n_efms = df_efms.shape[1]
    res.rank_failures = rank_failures

    if with_thermo:
        per_efm, net_df = thermo.compute_thermodynamics(
            df_norm, df_metabolites, model, config.get("PSEUDO_METS"),
            which=which, cc=cc)
        res.per_efm_thermo = per_efm
        res.net_conversions = net_df

    return res
