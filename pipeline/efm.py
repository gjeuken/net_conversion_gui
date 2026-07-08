"""Stoichiometric matrix, efmtool enumeration, rank validation, normalisation.

Ported from the notebook's EFM cells.  ``calculate_efms`` shells out to efmtool
(JVM required — see :func:`pipeline.efm.java_available`).
"""

from __future__ import annotations

import shutil

import numpy as np
import pandas as pd


def java_available() -> bool:
    """True if a JVM is on PATH (efmtool needs it)."""
    return shutil.which("java") is not None


def build_stoichiometry(model, rev_allowed, substrates):
    """Build ``(S, reversibility_array, reaction_ids, metabolite_ids)``.

    All-zero metabolite rows are dropped.  A reaction is reversible (1) iff its
    lower bound is negative; exchanges outside ``rev_allowed``/``substrates`` are
    forced irreversible first, mirroring the notebook.
    """
    import cobra

    S = cobra.util.array.create_stoichiometric_matrix(model)
    keep = ~np.all(S == 0, axis=1)
    metabolite_ids = np.array([m.id for m in model.metabolites])[keep]
    S = S[keep]

    rev_allowed = set(rev_allowed or [])
    substrates = set(substrates or [])
    reversibility = []
    for r in model.reactions:
        if r.id.startswith("EX") and r.id not in rev_allowed and r.id not in substrates:
            r.lower_bound = 0
        reversibility.append(1 if r.lower_bound < 0 else 0)

    reaction_ids = [r.id for r in model.reactions]
    return S, reversibility, reaction_ids, metabolite_ids


def calculate_efms(S, reversibility, reaction_ids, metabolite_ids):
    """Run efmtool and return a DataFrame (reactions x EFMs)."""
    import efmtool

    efms = efmtool.calculate_efms(
        np.asarray(S, dtype=float),
        list(reversibility),
        list(reaction_ids),
        list(metabolite_ids),
    )
    cols = [f"EFM_{i + 1}" for i in range(efms.shape[1])]
    return pd.DataFrame(efms, index=list(reaction_ids), columns=cols)


def validate_efms(df_efms, S):
    """Rank check: a valid EFM has ``n_active - rank(S_active) == 1``.

    Returns ``{efm_column: (n_active - rank)}`` for any that fail.
    """
    S = np.asarray(S, dtype=float)
    failures = {}
    for col in df_efms.columns:
        active = np.where(df_efms[col].values != 0)[0]
        if active.size == 0:
            failures[col] = -1
            continue
        S_active = S[:, active]
        S_active = S_active[~np.all(S_active == 0, axis=1), :]
        ncols = S_active.shape[1]
        rank = np.linalg.matrix_rank(S_active)
        if ncols - rank != 1:
            failures[col] = int(ncols - rank)
    return failures


def normalize_efms(df_efms, substrate_id):
    """Normalise every EFM by ``|flux through substrate_id|``."""
    denom = np.abs(df_efms.loc[substrate_id].values.astype(float))
    mask = denom != 0
    df = df_efms.copy().astype(float)
    df.loc[:, mask] = df_efms.loc[:, mask].astype(float) / denom[mask]
    return df


def pseudo_transport_reactions(model, pseudo_mets):
    """Non-exchange reactions whose whole metabolite set is currency mets."""
    pseudo_mets = set(pseudo_mets or [])
    return {
        r.id for r in model.reactions
        if not r.id.startswith("EX")
        and set(m.id for m in r.metabolites).issubset(pseudo_mets)
    }


def count_metabolic_reactions(df_normalized, model, pseudo_mets):
    """Active metabolic-reaction count per EFM (excludes EX + pseudo-transport).

    Returns ``(counts_series, sorted_excluded_pseudo_transport_ids)``.
    """
    pseudo_rxns = pseudo_transport_reactions(model, pseudo_mets)
    thermo_rows = {"dG0prime value", "dG0prime error", "dGm value", "dGm error",
                   "n_metabolic_reactions"}
    metabolic_only = df_normalized[
        ~df_normalized.index.str.startswith("EX")
        & ~df_normalized.index.isin(pseudo_rxns)
        & ~df_normalized.index.isin(thermo_rows)
    ]
    counts = (metabolic_only != 0).sum().rename("active_metabolic_reactions")
    return counts, sorted(pseudo_rxns)


def flux_sum_per_efm(df_normalized, model, pseudo_mets):
    """Σ|vᵢ| over metabolic reactions (for the Ω flux-sum denominator).

    Same exclusion set as :func:`count_metabolic_reactions`, but sums the
    absolute normalised fluxes rather than counting them.
    """
    pseudo_rxns = pseudo_transport_reactions(model, pseudo_mets)
    thermo_rows = {"dG0prime value", "dG0prime error", "dGm value", "dGm error",
                   "n_metabolic_reactions"}
    metabolic_only = df_normalized[
        ~df_normalized.index.str.startswith("EX")
        & ~df_normalized.index.isin(pseudo_rxns)
        & ~df_normalized.index.isin(thermo_rows)
    ].astype(float)
    return metabolic_only.abs().sum()
