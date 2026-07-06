"""eQuilibrator thermodynamics: net-conversion ΔG°′/ΔGm′ and the Ω measure.

``ComponentContribution`` downloads a large cache on first use (minutes,
network).  It is instantiated lazily and exactly once (CLAUDE.md §8) via
:func:`get_cc`; pass that instance into :func:`compute_thermodynamics`.

The Ω measure (CLAUDE.md §7), for an EFM normalised so the net-conversion rate
v_NC = 1::

    Ω = −ΔG_CAT / Σᵢ (vᵢ / v_NC)

where the flux sum runs over genuine metabolic reactions (exchanges and the
pseudo-transport reactions of currency metabolites excluded).  ΔG_CAT is the
net-conversion ΔG (ΔG°′ or ΔGm′, user's choice).  The separate
``equal_flux_bound`` column is ``−ΔG / n_metabolic_reactions`` — the upper bound
on Ω that coincides with Ω only when all fluxes are equal.
"""

from __future__ import annotations

import pandas as pd

from .efm import count_metabolic_reactions, flux_sum_per_efm
from .io import dict_to_kegg_reaction

_CC = None


def get_cc():
    """Return a process-wide ``ComponentContribution`` (instantiated once)."""
    global _CC
    if _CC is None:
        from equilibrator_api import ComponentContribution
        _CC = ComponentContribution()
    return _CC


def net_conversion_dict(efm_col, df_normalized, metabolite_map):
    """KEGG ``{kegg_id: signed_coeff}`` net conversion for one EFM column.

    Read from the exchange-flux rows: negative flux = uptake (reactant),
    positive = secretion (product).  Coefficients with the same KEGG id are
    summed so currency mets shared by several exchanges don't collide.
    """
    ex = df_normalized.loc[df_normalized.index.str.startswith("EX"), efm_col]
    kegg = {}
    for rid, v in ex.items():
        v = float(v)
        if abs(v) < 1e-9:
            continue
        metab = rid[2:]  # EX{ABBREV} — convention is load-bearing
        cid = metabolite_map.get(metab, metab)
        kegg[cid] = kegg.get(cid, 0.0) + v
    return {k: v for k, v in kegg.items() if abs(v) > 1e-9}


def net_conversion_string(efm_col, df_normalized, id_to_name=None):
    """Human-readable net conversion in working IDs (substrates -> products)."""
    ex = df_normalized.loc[df_normalized.index.str.startswith("EX"), efm_col]
    subs, prods = [], []
    for rid, v in ex.items():
        v = float(v)
        if abs(v) < 1e-9:
            continue
        abbrev = rid[2:]
        label = (id_to_name or {}).get(abbrev, abbrev)
        coeff = abs(v)
        term = label if abs(coeff - 1) < 1e-6 else f"{coeff:g} {label}"
        (subs if v < 0 else prods).append(term)
    return f"{' + '.join(subs) or '0'} -> {' + '.join(prods) or '0'}"


def _canonical_key(kegg_dict):
    return tuple(sorted((k, round(v, 6)) for k, v in kegg_dict.items()))


def compute_thermodynamics(df_normalized, df_metabolites, model, pseudo_mets,
                           which="dGm", cc=None):
    """Compute per-EFM and deduplicated net-conversion thermodynamics.

    Parameters
    ----------
    which : {"dGm", "dG0prime"}
        Which ΔG drives the Ω measure.

    Returns
    -------
    (per_efm_df, net_df)
        ``per_efm_df`` has one row per EFM; ``net_df`` collapses identical net
        conversions, noting EFM multiplicity and Ω range.
    """
    cc = cc or get_cc()
    metabolite_map = dict(zip(df_metabolites["ID"], df_metabolites["KEGG ID"]))
    # Display net conversions in the short working IDs (paper notation), not the
    # verbose compound names.
    id_to_name = dict(zip(df_metabolites["ID"], df_metabolites["ID"]))

    counts, _excluded = count_metabolic_reactions(df_normalized, model, pseudo_mets)
    flux_sums = flux_sum_per_efm(df_normalized, model, pseudo_mets)

    rows = []
    dg_cache = {}
    for col in df_normalized.columns:
        kegg = net_conversion_dict(col, df_normalized, metabolite_map)
        key = _canonical_key(kegg)
        if key not in dg_cache:
            dg_cache[key] = _equilibrator_dg(cc, kegg)
        dg = dg_cache[key]

        n_rxn = float(counts.get(col, float("nan")))
        flux_sum = float(flux_sums.get(col, float("nan")))
        dg_drive = dg["dGm value"] if which == "dGm" else dg["dG0prime value"]
        omega = (-dg_drive / flux_sum) if flux_sum and pd.notna(dg_drive) else None
        eq_bound = (-dg_drive / n_rxn) if n_rxn and pd.notna(dg_drive) else None

        rows.append({
            "EFM": col,
            "net_conversion": net_conversion_string(col, df_normalized, id_to_name),
            "_key": key,
            "dG0prime value": dg["dG0prime value"],
            "dG0prime error": dg["dG0prime error"],
            "dGm value": dg["dGm value"],
            "dGm error": dg["dGm error"],
            "thermo_error": dg["error"],
            "n_metabolic_reactions": n_rxn,
            "flux_sum": flux_sum,
            "Omega": omega,
            "equal_flux_bound": eq_bound,
        })

    per_efm = pd.DataFrame(rows)
    net_df = _deduplicate(per_efm)
    return per_efm.drop(columns=["_key"]), net_df


def _equilibrator_dg(cc, kegg_dict):
    """Query eQuilibrator for ΔG°′ and ΔGm′ of a net conversion."""
    out = {"dG0prime value": None, "dG0prime error": None,
           "dGm value": None, "dGm error": None, "error": None}
    if not kegg_dict:
        out["error"] = "empty net conversion"
        return out
    try:
        rxn = cc.parse_reaction_formula(dict_to_kegg_reaction(kegg_dict))
        dg0 = cc.standard_dg_prime(rxn)
        dgm = cc.physiological_dg_prime(rxn)
        out["dG0prime value"] = float(dg0.value.magnitude)
        out["dG0prime error"] = float(dg0.error.magnitude)
        out["dGm value"] = float(dgm.value.magnitude)
        out["dGm error"] = float(dgm.error.magnitude)
    except Exception as e:  # missing compound, parse failure, offline, ...
        out["error"] = str(e)
    return out


def _deduplicate(per_efm):
    """Collapse EFMs with identical net conversions into one row each."""
    if per_efm.empty:
        return per_efm
    out = []
    for key, grp in per_efm.groupby("_key", sort=False):
        first = grp.iloc[0]
        omegas = grp["Omega"].dropna().tolist()
        bounds = grp["equal_flux_bound"].dropna().tolist()
        out.append({
            "net_conversion": first["net_conversion"],
            "EFM_multiplicity": len(grp),
            "EFMs": ", ".join(grp["EFM"].tolist()),
            "dG0prime value": first["dG0prime value"],
            "dG0prime error": first["dG0prime error"],
            "dGm value": first["dGm value"],
            "dGm error": first["dGm error"],
            "Omega_min": min(omegas) if omegas else None,
            "Omega_max": max(omegas) if omegas else None,
            "equal_flux_bound_min": min(bounds) if bounds else None,
            "equal_flux_bound_max": max(bounds) if bounds else None,
            "thermo_error": first["thermo_error"],
        })
    return pd.DataFrame(out)
