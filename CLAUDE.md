# Build brief: Dash GUI for catabolic net conversions, EFMs, and ΔG

Paste this as the opening prompt in Claude Code (run it from an empty project
directory). It is a self-contained spec; ask me to clarify anything ambiguous
before writing large amounts of code.

---

## 1. Goal

Build a locally-runnable **Dash** web app that reproduces the computational
pipeline in the accompanying Jupyter notebook (`educational_paper_pipeline.ipynb`)
and the methods of the accompanying paper (*"A computation tool for substrate,
energy-carrier and product stoichiometries and thermodynamic feasibility of
catabolic pathways"*, Bruggeman et al.).

Given a user-supplied metabolic pathway, the app must:

0. Let the user assemble the network in a workbench — fetching reaction
   stoichiometries and IDs from KEGG, adding non-KEGG reactions manually, and
   iterating against a live balance check until the chemistry is consistent.
1. Check elemental + charge balance of every reaction (via KEGG).
2. Build a COBRApy model and run FBA sanity checks.
3. Enumerate **elementary flux modes (EFMs)** with efmtool.
4. For each EFM, compute the **net conversion** (overall stoichiometry).
5. Compute the standard Gibbs free energy change **ΔG°′** and physiological
   **ΔGm′** of each net conversion via eQuilibrator.
6. Compute the pathway thermodynamic measure **Ω** (flux-weighted ΔG per unit
   flux) defined in the paper.

This is an educational tool accompanying the paper. Prioritise clarity,
inspectable intermediate results, and honest error messages over feature
breadth. It must run for a single user on a laptop with no cloud services.

## 2. Stack

- Python 3.11+
- `dash` (use the built-in `background=True` callbacks with a `DiskcacheManager`
  for the long-running EFM step — do NOT block the main thread)
- `dash-bootstrap-components` for layout (or plain Dash if you prefer; keep it
  simple and uncluttered)
- `cobra` (COBRApy), `efmtool`, `equilibrator-api`, `numpy`, `pandas`,
  `requests`, `openpyxl`, `diskcache`
- A JVM is required by efmtool — see §8.

## 3. Project layout

```
netconv_app/
  app.py                 # Dash app, layout, callbacks, server entry point
  pipeline/
    __init__.py
    io.py                # any input source -> the two canonical dataframes; Excel save/load
    kegg.py              # fetch reaction EQUATION + compound formula/charge/name from KEGG REST
    balance.py           # formula/charge parsing, atom+charge balance checks
    model.py             # build COBRApy model, configure exchanges, FBA checks, prune blocked
    efm.py               # S matrix + reversibility, efmtool enumeration, rank validation, normalisation
    thermo.py            # eQuilibrator wrapper: net conversion ΔG°′/ΔGm′, Ω measure
    cache.py             # KEGG + eQuilibrator caching helpers
  assets/                # css
  examples/
    Example1_EMPglycolysis.xlsx     # ask user to provide; see §5 for schema
    pan_glycolysis.xlsx
  requirements.txt
  README.md
```

Port the notebook's functions into `pipeline/` largely as-is (they are correct
and tested) rather than rewriting from scratch. The functions to carry over:
`get_compound_info`, `normalize_arrow`, `parse_equation`, `parse_formula`,
`multiply_formula`, `map_custom_to_kegg_equation`, `check_balance`,
`analyze_custom_equations`, `analyze_stoichiometry_matrix`,
`dict_to_kegg_reaction`, `clean_dataframe_whitespace`, plus the model-building,
exchange-configuration, blocked-reaction-pruning, S-matrix/reversibility,
efmtool call, rank-check, and normalisation blocks. Keep each as a pure function
that takes data in and returns data out (no globals, no `print`; return
structured results the callbacks can render).

The one genuinely new module is `kegg.py`: a *reaction* fetcher (the notebook
only fetches compounds). Everything it produces lands in the same two canonical
dataframes (§5) that the rest of the pipeline already consumes, so nothing
downstream of input changes.

## 4. Pipeline (wire these to the UI as discrete, inspectable stages)

Mirror the notebook stages so a user can stop and fix the input between any two:

1. **Assemble the network** (the workbench, §5): the user builds the two
   canonical dataframes by any mix of KEGG reaction fetch, manual entry, and
   Excel upload, iterating against a live balance check until consistent.
2. **Balance check** → return a dataframe with per-reaction atom/charge balance,
   the imbalance summary string, and a clear pass/fail flag. Exchange reactions
   are expected to be unbalanced; flag only non-exchange imbalances as errors.
3. **Build model** (COBRApy): add metabolites, add reactions, set irreversible
   reactions' lower bound to 0.
4. **Configure exchanges**: from the user's SUBSTRATES / PRODUCTS / REV_ALLOWED,
   restrict or remove exchange reactions not in those sets.
5. **FBA sanity checks**:
   - max ATP (ENERGY_PRODUCT) with all substrate uptake blocked → must be ~0;
     a nonzero value means an energy-generating cycle. Surface the offending
     fluxes.
   - max individual yield of each CARBON_PRODUCT on the substrate.
6. **Prune blocked reactions** (`cobra.flux_analysis.find_blocked_reactions`).
7. **Build S + reversibility array**, drop all-zero rows.
8. **EFM enumeration** (efmtool) — *background callback*, with a progress/elapsed
   indicator. efmtool time scales super-exponentially; warn the user and allow
   cancellation. Report number of EFMs.
9. **Validate** EFMs (rank check: `n_active - rank == 1`); flag any failures.
10. **Normalise** each EFM by |substrate uptake flux|.
11. **Net conversions**: read exchange-flux rows per EFM → overall stoichiometry,
    expressed as yields on the normalised substrate.
12. **Thermodynamics** per EFM: ΔG°′ and ΔGm′ via eQuilibrator (value + error),
    and **Ω** (§7). Many EFMs share a net conversion — deduplicate net
    conversions for the thermo display and note the EFM multiplicity.

## 5. Input — the workbench (primary feature)

The user does not need a finished spreadsheet to start. Input is an **iterative
workbench**: they accumulate reactions from any mix of sources into two
canonical dataframes, and a **live balance check** tells them what is still
inconsistent and when the network is ready for efmtool. This loop is the heart
of the app — build it well.

### The two canonical dataframes (the single source of truth)

Every input source writes into these; the whole downstream pipeline reads only
these. Keep the exact schema — it matches the notebook and the paper's
supplementary files:

- **Metabolites**: columns `ID`, `Name`, `KEGG ID`. `ID` is the working
  abbreviation; `KEGG ID` like `C00031`.
- **Reactions**: columns `ID`, `Name`, `Reaction stoichiometry`, `Reversibility`
  (1 = reversible, 0 = irreversible). Stoichiometry is a text string in the
  working IDs using any common arrow (`->`, `<=>`, `→`).
- Exchange reactions are named `EX{METABOLITE_ID}` exactly — the thermo step
  strips the `EX` prefix to map back to KEGG, so this convention is load-bearing.

### Three input sources, all feeding the same dataframes

1. **KEGG reaction fetch (primary).** User pastes/loads a list of KEGG reaction
   IDs (e.g. `R00200`, `R01786`) — or, as a convenience, a KEGG **module** ID
   (e.g. `M00001` for EMP) to pull a whole reaction set in one go. For each
   reaction, `kegg.py` fetches the `EQUATION` line from `rest.kegg.jp/get/rn:Rxxxxx`,
   parses it (reuse `parse_equation`/`normalize_arrow`), and for every compound
   ID auto-fills the Metabolites table — the KEGG ID *is* the compound ID, and
   names come from `rest.kegg.jp/get/cpd:Cxxxxx`. This removes both manual
   stoichiometry entry and the ID↔KEGG mapping in one step.
2. **Manual entry / edit.** Add or edit a row directly: type an ID, a
   stoichiometry string, set reversibility. This is the escape hatch for
   everything KEGG cannot give (see §8): transport, membrane/ion translocation,
   PTS-style import, and designed reactions (the paper's xylose co-consumption
   and Calvin-shunt networks have no KEGG IDs).
3. **Excel upload.** Load a previously saved workbook to resume. Excel is
   **save/load**, not the only door in.

Always offer **download to `.xlsx`** of the current dataframes so a session is a
reproducible artifact that travels with the paper's supplements.

### The live balance loop (the convergence signal)

After *any* change from any source, re-run the balance check
(`analyze_custom_equations` / `analyze_stoichiometry_matrix`) and show a
per-reaction status with the `atom_imbalance` string right next to each row
("+1 H, −1 charge"). Provide one-click fixes that append H⁺ / H₂O to either side
of a flagged reaction. Exchange reactions are expected to be unbalanced — flag
only non-exchange imbalances as errors. This single feedback loop is what makes
curation convergent: the user always knows what is left to fix and when the
network is clean enough to model.

### Run-time selectors (the notebook's USER INPUTS cell)

These become form controls, populated from the metabolites/exchanges actually
present: `MODEL_NAME`, `SUBSTRATES`, `PRODUCTS`, `CARBON_PRODUCTS`,
`ENERGY_PRODUCT`, `REV_ALLOWED`, `PSEUDO_METS`. Pre-flag known currency
metabolites (ATP `C00002`, ADP `C00008`, Pi `C00009`, H₂O `C00001`, H⁺ `C00080`)
to auto-suggest `PSEUDO_METS` and external-metabolite candidates — but these
remain user-confirmable modelling choices, never silent defaults.

### Scientific judgments the workbench must NOT make for the user

The balance loop converges for *chemistry only*. Two choices stay explicitly with
the user, exactly as the paper stresses:

- **Reversibility and direction.** KEGG writes every reaction `<=>` in an
  arbitrary orientation. For reactions the user marks irreversible, they must
  also confirm/flip the physiological forward direction — do not auto-derive it.
  You may *offer* an eQuilibrator-ΔG-based suggestion of likely-irreversible
  reactions, clearly labelled as a suggestion.
- **Which metabolites are external.** The user designates the boundary; the app
  only supplies candidates.

Make the mechanics frictionless without implying the app has made these calls.

## 6. UI

Single page, stepwise, uncluttered. Suggested sections (tabs or a vertical
stepper):

- **Workbench (Input + Balance combined)**: the KEGG-fetch box (reaction or
  module IDs), an editable reactions/metabolites `DataTable` for manual entry,
  Excel upload/download, and the run-time selectors — all above a live
  per-reaction balance panel showing pass/fail and the `atom_imbalance` string
  with one-click H⁺/H₂O fixes. Keep Input and Balance on one screen: the user
  edits and sees consistency update together. Default to a loaded EMP-glycolysis
  → lactate example so the app is non-empty on first open.
- **Model & FBA**: blocked-reaction list, ATP-without-substrate result (green if
  ~0, red otherwise), per-product max yields.
- **EFMs**: count, the normalised EFM table (reactions × EFMs), rank-check status.
- **Net conversions & thermodynamics**: one row per distinct net conversion —
  the overall equation, EFM multiplicity, ΔG°′ ± err, ΔGm′ ± err, Ω, and the
  n_R approximation (−ΔG/n_R). Let the user pick which ΔG drives Ω (ΔG°′ vs ΔGm′).
- **Downloads**: COBRApy JSON, SBML, and a results CSV/XLSX (EFM table + net
  conversions + thermo). Matches the paper's "CSV export" expectation.

Keep ΔG sign conventions and units explicit in headers (kJ/mol). Reuse the
paper's notation (ΔG°′, ΔGm′, Ω, ΔG_CAT).

## 7. The Ω measure — implement the real definition

The paper defines, for an EFM with steady-state rates vᵢ normalised so the net
conversion rate v_NC = 1:

    Ω = −ΔG_CAT / Σᵢ (vᵢ / v_NC)
      = Σᵢ ωᵢ (−ΔGᵢ),   ωᵢ = vᵢ / Σⱼ vⱼ

i.e. ΔG_CAT divided by the **flux sum**, not the reaction count. Compute Ω from
the net-conversion ΔG (ΔG_CAT) and the sum of the normalised metabolic fluxes
(exclude exchange reactions and the pseudo-transport reactions of currency
metabolites in `PSEUDO_METS`, exactly as the notebook excludes them when counting
reactions). Ω equals the MDF result with no concentration bounds, so label it as
such in the UI.

Also display the notebook's `−ΔG / n_metabolic_reactions` as a separate column
labelled clearly as the **equal-flux upper bound** on Ω — they coincide only when
all fluxes are equal. Do not conflate the two.

(If you later want it: ΔG_CAT = Σᵢ (vᵢ/v_NC) ΔGᵢ, so Ω can also be assembled from
per-reaction ΔGᵢ if those are fetched — but the net-conversion ΔG from
eQuilibrator is sufficient and avoids per-reaction ΔG uncertainty, which is the
whole point of Ω.)

## 8. Environment gotchas — handle these explicitly

- **efmtool needs Java.** Detect the JVM at startup and print a clear message
  (with install hint) if missing, rather than letting JPype throw a cryptic
  error. Pin a working `efmtool` version in requirements and note the Java
  requirement in the README.
- **eQuilibrator first-run cache.** `ComponentContribution()` downloads a large
  cache on first instantiation (minutes, needs network). Instantiate it **once**,
  lazily, on first thermo request — not at import — and show a "downloading
  thermodynamic data, first run only" notice. Cache it for the app lifetime.
- **KEGG REST.** Keep the `@lru_cache` on compound lookups and add a small
  on-disk cache (JSON/sqlite keyed by ID) covering **both** compound *and*
  reaction entries, so repeat runs and re-fetches don't re-hit KEGG. Handle
  timeouts/HTTP errors gracefully per entry.
- **KEGG reaction parsing limits — surface, don't paper over.** Fetched
  reactions frequently (a) omit H⁺/H₂O, (b) aren't atom/charge balanced as
  written, or (c) contain generic R-groups or polymer coefficients (`n`, `n+1`,
  `(n)`) that have no numeric stoichiometry. For (a)/(b) the live balance loop is
  the safety net — expect more flags than with hand-curated input. For (c),
  detect the non-numeric/R-group case explicitly, do **not** guess, flag the row
  ("R-group / polymer coefficient — enter manually"), and route the user to
  manual entry. KEGG has no compartments or transport, so PTS-style import and
  ion-translocation reactions must always be added by hand.
- **Network.** The app needs outbound network for KEGG and the eQuilibrator
  cache download. State this in the README. If offline, balance-check and thermo
  steps should degrade gracefully (skip with a warning), but EFM enumeration and
  FBA still work fully offline.
- **Long EFM runs.** Use Dash background callbacks (DiskcacheManager) + a cancel
  button + an EFM-count cap with a warning, since networks with bypasses explode
  combinatorially (the paper notes ecolicore → ~272M EFMs).

## 9. Deployment

Local-first: `python app.py` serving `http://127.0.0.1:8050`. Include a
`requirements.txt` and a README with: Java install, `pip install -r`, first-run
eQuilibrator note, and how to load the example. A Dockerfile (bundling a JRE) is
a nice-to-have for reproducibility but optional.

## 10. Acceptance test

The app is correct when, loading the **EMP-glycolysis-to-lactate** example, it
reproduces the paper's headline result:

    1 GLC + 2 ADP + 2 Pi → 2 LAC + 2 ATP

as the single net conversion, with a sensible negative ΔGm′ for the net
conversion. As a stretch test, the **pan-glycolysis** network in the paper
(two imports × EMP/ED × GAPDH/GAPN × three fermentation routes) should enumerate
**24 EFMs** collapsing to a smaller set of distinct net conversions.

## 11. Build order

Build in this order so each layer is verifiable before the next:

1. **Compute core behind a thin CLI.** Port the notebook functions; reproduce the
   EMP→lactate acceptance test (§10) end-to-end from a static Excel *before* any
   UI or KEGG-fetch work. This de-risks efmtool/eQuilibrator first.
2. **KEGG reaction fetch** (`kegg.py`) as a standalone function: reaction ID →
   rows in the two dataframes, with the R-group/polymer/unbalanced cases handled
   per §8. Test it independently against a known reaction (e.g. `R00200`).
3. **The workbench UI**: editable tables + KEGG-fetch box + Excel save/load,
   wired to the live balance loop as the primary input mode.
4. **The rest of the pipeline UI**: model/FBA, EFM (background callback), net
   conversions, thermodynamics, downloads.

Surface every intermediate dataframe in the UI — this is a teaching tool; the
intermediate state is the point. Do the full editable two-table grid only as part
of the workbench (step 3); don't let it expand into a general spreadsheet clone.
