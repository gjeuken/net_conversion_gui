"""Thin CLI to run the full pipeline on an Excel workbook.

    python -m netconv_app.cli netconv_app/examples/Example1_EMP_lactate.xlsx

Use ``--no-thermo`` to skip eQuilibrator (works fully offline, no first-run
cache download).  ``--which {dGm,dG0prime}`` selects the ΔG driving Ω.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from pipeline import efm, io, run

pd.set_option("display.width", 160)
pd.set_option("display.max_columns", 40)


def _default_config(df_reactions):
    """Fallback config if the workbook has no Config sheet."""
    ex = [r for r in df_reactions["ID"] if str(r).startswith("EX")]
    subs = [r for r in ex if "GLC" in r] or ex[:1]
    return {
        "MODEL_NAME": "model",
        "SUBSTRATES": subs,
        "PRODUCTS": [r for r in ex if r not in subs],
        "CARBON_PRODUCTS": [],
        "ENERGY_PRODUCT": next((r for r in ex if "ATP" in r), ex[0]),
        "REV_ALLOWED": [r for r in ex if any(
            k in r for k in ("ADP", "Pi", "H2O", "EXH"))],
        "PSEUDO_METS": set(),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("workbook", help="Excel file with Metabolites/Reactions sheets")
    ap.add_argument("--no-thermo", action="store_true",
                    help="skip eQuilibrator ΔG / Ω computation")
    ap.add_argument("--which", choices=["dGm", "dG0prime"], default="dGm",
                    help="which ΔG drives the Ω measure")
    args = ap.parse_args(argv)

    if not efm.java_available():
        print("ERROR: Java (JVM) not found on PATH — efmtool needs it.\n"
              "Install a JRE (e.g. `sudo apt install default-jre`) and retry.",
              file=sys.stderr)
        return 2

    df_met, df_rxn = io.load_excel(args.workbook)
    config = io.load_config(args.workbook) or _default_config(df_rxn)
    print(f"Loaded {len(df_met)} metabolites, {len(df_rxn)} reactions")
    print(f"Config: {config}\n")

    res = run.run_all(df_met, df_rxn, config,
                      with_thermo=not args.no_thermo, which=args.which)

    print("=" * 70)
    print("STAGE 2 — Balance check")
    print("=" * 70)
    cols = ["reaction_id", "atom_balanced", "charge_balanced", "is_exchange",
            "atom_imbalance", "notes"]
    print(res.df_balance[cols].to_string(index=False))
    print(f"\nBalance OK (non-exchange): {res.balance_ok}")
    if res.balance_errors:
        print("Imbalanced non-exchange reactions:")
        for e in res.balance_errors:
            print(f"  {e['reaction_id']}: {e['imbalance']} [{e['notes']}]")

    print("\n" + "=" * 70)
    print("STAGE 5 — FBA sanity checks")
    print("=" * 70)
    print(f"Removed exchanges: {res.removed_exchanges}")
    print(f"Max {config['ENERGY_PRODUCT']} with no substrate: "
          f"{res.atp_without_substrate:.6g}  "
          f"({'OK ~0' if abs(res.atp_without_substrate) < 1e-6 else 'WARNING: cycle!'})")
    if res.atp_cycle_fluxes:
        print(f"  offending fluxes: {res.atp_cycle_fluxes}")
    print("Max product yields:")
    for p, y in res.product_yields.items():
        print(f"  {p[2:] if p.startswith('EX') else p}: "
              f"{'infeasible' if y is None else f'{y:.4g}'}")

    print("\n" + "=" * 70)
    print("STAGE 8-10 — EFM enumeration")
    print("=" * 70)
    print(f"Blocked reactions removed: {len(res.blocked_reactions)}")
    print(f"Number of EFMs: {res.n_efms}")
    print(f"Rank-check failures: {len(res.rank_failures)}"
          + (f" {res.rank_failures}" if res.rank_failures else ""))

    print("\n" + "=" * 70)
    print("STAGE 11-12 — Net conversions & thermodynamics")
    print("=" * 70)
    if res.net_conversions is not None and not res.net_conversions.empty:
        nc = res.net_conversions.copy()
        for c in ["dG0prime value", "dG0prime error", "dGm value", "dGm error",
                  "Omega_min", "Omega_max", "equal_flux_bound_min",
                  "equal_flux_bound_max"]:
            if c in nc.columns:
                nc[c] = nc[c].map(lambda v: f"{v:.2f}" if isinstance(v, (int, float)) and pd.notna(v) else v)
        show = ["net_conversion", "EFM_multiplicity", "dG0prime value",
                "dGm value", "Omega_min", "Omega_max", "equal_flux_bound_min"]
        print(nc[show].to_string(index=False))
        if nc["thermo_error"].notna().any():
            print("\nThermo notes:")
            for _, r in nc.iterrows():
                if pd.notna(r["thermo_error"]):
                    print(f"  {r['net_conversion']}: {r['thermo_error']}")
    else:
        print("(thermo skipped)")
        for col in res.df_normalized.columns:
            from pipeline.thermo import net_conversion_string
            print(f"  {col}: {net_conversion_string(col, res.df_normalized, dict(zip(df_met['ID'], df_met['ID'])))}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
