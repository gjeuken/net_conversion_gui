# Catabolic net conversions, EFMs & ΔG — a Dash GUI

A locally-runnable [Dash](https://dash.plotly.com/) web app that reproduces the
computational pipeline of the paper *"A computation tool for substrate,
energy-carrier and product stoichiometries and thermodynamic feasibility of
catabolic pathways"* (Bruggeman et al.) and its companion notebook
`educational_paper_pipeline.ipynb`.

You assemble a metabolic pathway in an interactive **workbench**, and the app
walks it through the full analysis, surfacing every intermediate result:

1. **Balance** — elemental + charge balance of every reaction (compound
   formulas/charges fetched from KEGG).
2. **Model & FBA** — builds a [COBRApy](https://opencobra.github.io/cobrapy/)
   model and runs sanity checks (no ATP without substrate; per-product yields).
3. **EFMs** — enumerates the **elementary flux modes** with
   [efmtool](https://csb.ethz.ch/tools/software/efmtool.html).
4. **Net conversions** — the overall stoichiometry of each EFM.
5. **Thermodynamics** — **ΔG°′** and physiological **ΔGm′** of each net
   conversion via [eQuilibrator](https://equilibrator.weizmann.ac.il/), plus the
   pathway measure **Ω** (flux-weighted ΔG per unit flux).

It is an *educational* tool: the intermediate dataframes are the point, so you
can stop and fix the input between any two stages.

---

## Quick start

```bash
# 1. Create an environment and install dependencies
python -m venv ~/envs/catabolic
~/envs/catabolic/bin/pip install -r netconv_app/requirements.txt

# 2. Run the app (from the repo root — the directory containing netconv_app/)
~/envs/catabolic/bin/python -m netconv_app.app
```

Open <http://127.0.0.1:8050>. The app opens pre-loaded with the
**EMP-glycolysis → lactate** example, so it is non-empty on first run.

### Requirements

* **Python 3.11+**
* **A Java runtime (JRE).** efmtool runs on the JVM. Without it the EFM step
  fails with a clear message; the footer of the app shows whether a JVM was
  detected. Install one with e.g.:
  ```bash
  sudo apt install default-jre      # Debian/Ubuntu
  brew install openjdk              # macOS
  ```
  Verify with `java -version`.
* **Outbound network** for KEGG lookups and the eQuilibrator first-run cache
  download. Most steps degrade gracefully offline — see [Offline
  behaviour](#offline-behaviour).

> **First-run eQuilibrator note.** The first time you click *Compute net
> conversions & ΔG*, `ComponentContribution()` downloads a large thermodynamic
> cache (minutes, needs network). This happens **once** per machine and is
> reused afterwards. The app instantiates it lazily — balance, model, FBA and
> EFM steps never wait on it.

---

## Using the app

The page is a five-tab stepper. You generally move left to right, but every tab
re-reads the current state, so you can go back to the workbench, edit, and
re-run any downstream stage.

### ① Workbench — assemble & balance

The heart of the tool. You build **two canonical tables** and the whole
pipeline reads only these:

* **Reactions** — `ID`, `Name`, `Reaction stoichiometry` (a text string in your
  working IDs using any arrow: `->`, `<=>`, `→`), `Reversibility` (1 = reversible,
  0 = irreversible).
* **Metabolites** — `ID` (working abbreviation), `Name`, `KEGG ID` (e.g.
  `C00031`), and an optional `Chemical formula`. Balance checking normally pulls
  each compound's formula/charge from KEGG, but a metabolite with **no KEGG ID**
  (a designed or non-KEGG compound) can't be looked up — supply its
  `Chemical formula` (e.g. `C4H8O2`) and the balance check uses that directly.
  When present, an explicit formula overrides the KEGG lookup. (Charge defaults
  to 0 for formula-only metabolites; add an optional `Charge` column if you need
  a non-zero one.)

You can populate them by any mix of:

* **Manual entry** — edit the tables directly; *+ Add row* appends a blank row.
  This is the escape hatch for everything KEGG can't give: transport,
  membrane/ion translocation, PTS-style import, and designed reactions.
* **Excel upload** — resume a saved session.
* **Download workbook (.xlsx)** — save the current state as a reproducible
  artifact (it round-trips through upload, and carries your run-time selectors
  in an optional `Config` sheet).

A **live balance panel** re-checks every reaction after each edit and colours
each row, showing the atom/charge imbalance string (e.g. `+1 H, −1 charge`)
right next to it. **Exchange reactions** — named `EX{METABOLITE_ID}` — are
expected to be unbalanced and are *not* flagged. Use the **one-click H⁺ / H₂O
fixes** to patch a flagged reaction from either side.

> The `EX` naming convention is **load-bearing**: the thermodynamics step strips
> the `EX` prefix to map an exchange back to its KEGG compound.

**Run-time selectors** (substrates, products, carbon products, energy product,
freely-reversible exchanges, and pseudo/currency metabolites) populate from the
exchanges actually present. Currency metabolites (ATP, ADP, Pi, H₂O, H⁺) are
pre-suggested — but which metabolites are external, and reaction
reversibility/direction, remain **your** modelling choices. The app never makes
those scientific calls silently.

### ② Model & FBA

Builds the COBRApy model and reports:

* the **blocked-reaction** list,
* the **ATP-without-substrate** check (green if ≈ 0; red means an
  energy-generating cycle exists, and the offending fluxes are surfaced),
* per-product **max yields** on the substrate.

### ③ EFMs

Runs efmtool as a **background task** (with a **Cancel** button and an EFM cap,
since networks with bypasses explode combinatorially — the paper notes
ecolicore → ~272M EFMs). Reports the EFM count, the normalised EFM table
(reactions × EFMs, per unit substrate uptake), and the **rank-check** status
(`n_active − rank = 1`).

### ④ Net conversions & ΔG

One row per **distinct** net conversion (EFMs sharing a net conversion are
deduplicated, noting the multiplicity), with the overall equation, ΔG°′ ± err,
ΔGm′ ± err, **Ω**, and the equal-flux upper bound −ΔG/n_R. A selector lets you
choose whether **ΔGm′** or **ΔG°′** drives Ω.

### ⑤ Downloads

COBRApy **JSON**, **SBML**, and a **results workbook** (EFM table + net
conversions + per-EFM thermodynamics).

---

## The Ω measure

For an EFM normalised so the net-conversion rate *v*<sub>NC</sub> = 1:

```
Ω = −ΔG_CAT / Σᵢ (vᵢ / v_NC)
```

i.e. the net-conversion ΔG divided by the **sum of normalised metabolic
fluxes** — *not* the reaction count. Exchanges and the pseudo-transport
reactions of currency metabolites (`PSEUDO_METS`) are excluded from the sum,
exactly as the notebook excludes them. Ω equals the MDF result with no
concentration bounds.

The separate **equal-flux bound** column, −ΔG / n_metabolic_reactions,
coincides with Ω only when all fluxes are equal — the two are reported
separately and never conflated.

---

## Command-line use

The full pipeline also runs headless — handy for the acceptance tests and for
scripting:

```bash
# Reproduces the headline result: 1 GLC + 2 ADP + 2 Pi -> 2 LAC + 2 ATP (+ 2 H2O)
~/envs/catabolic/bin/python -m netconv_app.cli \
    netconv_app/examples/Example1_EMP_lactate.xlsx

# 24 EFMs; skip eQuilibrator so it runs fully offline
~/envs/catabolic/bin/python -m netconv_app.cli \
    netconv_app/examples/pan_glycolysis.xlsx --no-thermo
```

Flags: `--no-thermo` skips eQuilibrator (and its first-run download);
`--which {dGm,dG0prime}` selects which ΔG drives Ω (default `dGm`).

---

## Examples

| File | What it is |
| --- | --- |
| `netconv_app/examples/Example1_EMP_lactate.xlsx` | EMP glycolysis → lactate (the default; reproduces the paper's headline net conversion) |
| `netconv_app/examples/Example1_EMPglycolysis.xlsx` | EMP glycolysis → ethanol / acetate / formate |
| `netconv_app/examples/pan_glycolysis.xlsx` | Pan-glycolysis (2 imports × EMP/ED × GAPDH/GAPN × 3 fermentation routes) — the 24-EFM stretch test |
| `netconv_app/examples/BDO_production_pathways.xlsx` | 1,4-butanediol production pathways — demonstrates a non-KEGG metabolite (`H4BUAL`) balanced via an explicit `Chemical formula` |

Regenerate them with:

```bash
~/envs/catabolic/bin/python netconv_app/examples/_generate_examples.py
```

### Acceptance tests

* **EMP-glycolysis → lactate** yields the single net conversion
  **1 GLC + 2 ADP + 2 Pi → 2 LAC + 2 ATP** (the app also shows the chemically
  explicit + 2 H₂O from ATP synthesis) with a sensibly negative ΔGm′.
* **Pan-glycolysis** enumerates exactly **24 EFMs**, collapsing to a smaller set
  of distinct net conversions.

---

## Offline behaviour

* **EFM enumeration and FBA** work fully offline.
* **Balance checking** degrades gracefully: compounds that can't be fetched are
  treated as unknown (a warning, not a crash). Cached KEGG compound/reaction
  entries live in `netconv_app/.cache/kegg_cache.json`.
* **Thermodynamics** needs the eQuilibrator cache; skip it with `--no-thermo`
  in the CLI, or simply don't open tab ④.

KEGG reactions frequently omit H⁺/H₂O, aren't balanced as written, or contain
generic R-groups / polymer coefficients (`n`, `n+1`). The first two are caught
by the live balance loop; R-group/polymer rows are flagged explicitly and
routed to manual entry rather than guessed.

---

## Project layout

```
netconv_app/
  app.py                 # Dash app: layout + callbacks (serves 127.0.0.1:8050)
  cli.py                 # headless pipeline runner
  pipeline/
    io.py                # canonical dataframes, Excel + config save/load
    kegg.py              # KEGG REST: reaction + compound fetch (formula/charge/name)
    balance.py           # formula/charge parsing, atom + charge balance
    model.py             # COBRApy model, exchanges, FBA, pruning
    efm.py               # S matrix, efmtool, rank check, normalisation
    thermo.py            # eQuilibrator ΔG, Ω measure
    cache.py             # on-disk KEGG cache
    run.py               # stage orchestration (shared by app + CLI)
  examples/              # example workbooks + their generator
  assets/style.css
  requirements.txt
CLAUDE.md                # the original build brief / specification
educational_paper_pipeline.ipynb   # the reference notebook
```

---

## Notes & conventions

* ΔG sign conventions and units (**kJ/mol**) are stated explicitly in the UI
  headers, reusing the paper's notation (ΔG°′, ΔGm′, Ω, ΔG_CAT).
* eQuilibrator prints a harmless `AttributeError` from its `__del__` at
  interpreter shutdown (a known cleanup-ordering quirk); it appears after
  results are already returned and does not affect any run.
