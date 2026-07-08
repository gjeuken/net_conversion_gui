# Catabolic net conversions, EFMs & ΔG — Dash GUI

A locally-runnable [Dash](https://dash.plotly.com/) app that reproduces the
computational pipeline of the paper *"A computation tool for substrate,
energy-carrier and product stoichiometries and thermodynamic feasibility of
catabolic pathways"* (Bruggeman et al.) and the accompanying notebook
`educational_paper_pipeline.ipynb`.

Given a metabolic pathway you assemble in the workbench, the app:

1. Checks elemental + charge balance of every reaction (via KEGG).
2. Builds a COBRApy model and runs FBA sanity checks.
3. Enumerates **elementary flux modes (EFMs)** with efmtool.
4. Computes each EFM's **net conversion** (overall stoichiometry).
5. Computes **ΔG°′** and physiological **ΔGm′** of each net conversion via
   eQuilibrator.
6. Computes the pathway thermodynamic measure **Ω** (flux-weighted ΔG per unit
   flux).

It is a computational tool: every intermediate dataframe is surfaced so you can
stop and fix the input between stages.

## Requirements

* **Python 3.11+**
* **A Java runtime (JVM).** efmtool runs on the JVM. Without it the EFM step
  fails with a clear message. The **conda** install below bundles it
  (`openjdk`); otherwise install one system-wide, e.g.:
  ```
  sudo apt install default-jre      # Debian/Ubuntu
  brew install openjdk              # macOS
  ```
  Verify with `java -version`.
* **Outbound network** for KEGG (compound lookups during balance checks) and the
  eQuilibrator first-run cache download. See *Offline behaviour* below.

## Install

**With conda (recommended — bundles Java, no separate JRE install).** Run from
the repo root so the `pip -r` path in `environment.yml` resolves:

```
conda env create -f environment.yml
conda activate catabolic
```

**With venv** (needs Java installed separately, see Requirements):

```
python -m venv ~/envs/catabolic        # or use your own env
~/envs/catabolic/bin/pip install -r netconv_app/requirements.txt
```

## Run

From the **project root** (the directory containing `netconv_app/`):

```
~/envs/catabolic/bin/python -m netconv_app.app
```

Then open <http://127.0.0.1:8050>. The app opens pre-loaded with the
**EMP-glycolysis → lactate** example.

> **First-run eQuilibrator note.** The first time you click *Compute net
> conversions & ΔG*, `ComponentContribution()` downloads a large thermodynamic
> cache (minutes, needs network). This happens **once** for the machine; the
> cache is reused afterwards. The app instantiates it lazily — balance, model,
> FBA and EFM steps never wait for it.

## Using the app

The page is a five-tab stepper:

1. **① Workbench** — the heart of the tool. Build the two canonical tables
   (Metabolites, Reactions) by any mix of:
   * **Manual entry**: edit the tables directly, *+ Add … row* for new rows.
     Enter transport / ion-translocation / designed reactions as needed.
   * **Excel upload** to resume a saved session; **Download workbook** to save
     one (a reproducible `.xlsx` that travels with the paper's supplements).

   A **live balance panel** re-checks every reaction after each change and
   colours rows green / yellow / red. Exchange reactions (`EX…`) are expected to
   be unbalanced and are not flagged. Use the **one-click H⁺/H₂O fixes** to
   patch a flagged reaction. The **run-time selectors** (substrates, products,
   energy product, reversible exchanges, pseudo/currency metabolites) populate
   from the exchanges present — these are *your* modelling choices, never silent
   defaults.

2. **② Model & FBA** — blocked-reaction list, the ATP-without-substrate check
   (green if ≈ 0, red if an energy-generating cycle exists), and per-product max
   yields.

3. **③ EFMs** — runs efmtool as a **background task** (with Cancel and an EFM
   cap). Shows the EFM count, the normalised EFM table, and the rank-check
   status.

4. **④ Net conversions & ΔG** — one row per **distinct** net conversion, with
   EFM multiplicity, ΔG°′ ± err, ΔGm′ ± err, **Ω**, and the equal-flux upper
   bound −ΔG/n_R. Pick whether ΔGm′ or ΔG°′ drives Ω.

5. **⑤ Downloads** — COBRApy JSON, SBML, and a results workbook (EFM table + net
   conversions + per-EFM thermo).

### The Ω measure

For an EFM normalised so the net-conversion rate v_NC = 1:

```
Ω = −ΔG_CAT / Σᵢ (vᵢ / v_NC)
```

i.e. the net-conversion ΔG divided by the **sum of normalised metabolic
fluxes** (exchanges and the pseudo-transport reactions of currency metabolites
in `PSEUDO_METS` are excluded). Ω equals the MDF result with no concentration
bounds. The separate **equal-flux bound** column −ΔG/n_metabolic_reactions
coincides with Ω only when all fluxes are equal — they are not conflated.

## Command-line use

The full pipeline also runs headless, handy for the acceptance tests:

```
# Reproduces 1 GLC + 2 ADP + 2 Pi -> 2 LAC + 2 ATP (+ 2 H2O)
python -m netconv_app.cli netconv_app/examples/Example1_EMP_lactate.xlsx

# 24 EFMs; skip eQuilibrator (fully offline)
python -m netconv_app.cli netconv_app/examples/pan_glycolysis.xlsx --no-thermo
```

## Offline behaviour

* **EFM enumeration and FBA** work fully offline.
* **Balance checking** degrades gracefully: compounds that can't be fetched are
  treated as unknown (a warning, not a crash). Cached KEGG compound entries live
  in `netconv_app/.cache/kegg_cache.json`.
* **Thermodynamics** needs the eQuilibrator cache; skip it with `--no-thermo`
  in the CLI, or simply don't open tab ④.

## Project layout

```
netconv_app/
  app.py                 # Dash app: layout + callbacks
  cli.py                 # headless pipeline runner
  pipeline/
    io.py                # canonical dataframes, Excel + config save/load
    kegg.py              # KEGG REST: compound lookups (formula/charge/name)
    balance.py           # formula/charge parsing, atom+charge balance
    model.py             # COBRApy model, exchanges, FBA, pruning
    efm.py               # S matrix, efmtool, rank check, normalisation
    thermo.py            # eQuilibrator ΔG, Ω measure
    cache.py             # on-disk KEGG cache
    run.py               # stage orchestration (shared by app + CLI)
  examples/
    Example1_EMP_lactate.xlsx        # default
    Example1_EMPglycolysis.xlsx      # ethanol/acetate/formate net
    pan_glycolysis.xlsx              # 24-EFM stretch test
    _generate_examples.py            # regenerates the above
  assets/style.css
  requirements.txt
```

## Reproducing the examples

```
python netconv_app/examples/_generate_examples.py
```
