"""Dash GUI for catabolic net conversions, EFMs and ΔG.

Run locally::

    python -m netconv_app.app      # or: python netconv_app/app.py

Serves http://127.0.0.1:8050.  A JVM is required for the EFM step and outbound
network for KEGG balance checks and the eQuilibrator first-run cache download.
See the README.
"""

from __future__ import annotations

import base64
import io as _io
import os
import tempfile

import dash
import dash_bootstrap_components as dbc
import pandas as pd
from dash import (ALL, MATCH, Dash, DiskcacheManager, Input, Output, State, ctx,
                  dash_table, dcc, html, no_update)

import diskcache

from pipeline import efm, io, model as model_mod, run, thermo

# --------------------------------------------------------------------------- #
# App + background-callback manager (DiskcacheManager for the long EFM step)
# --------------------------------------------------------------------------- #
_CACHE_DIR = os.path.join(tempfile.gettempdir(), "netconv_dash_cache")
_dc = diskcache.Cache(_CACHE_DIR)
background_manager = DiskcacheManager(_dc)

app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    background_callback_manager=background_manager,
    suppress_callback_exceptions=True,
    title="Catabolic net conversions, EFMs & ΔG",
)
server = app.server

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), "examples")
DEFAULT_EXAMPLE = os.path.join(EXAMPLE_DIR, "Example1_EMP_lactate.xlsx")
EXAMPLES = {
    "EMP glycolysis → lactate (default)": "Example1_EMP_lactate.xlsx",
    "EMP glycolysis → ethanol/acetate/formate": "Example1_EMPglycolysis.xlsx",
    "Pan-glycolysis (24 EFMs)": "pan_glycolysis.xlsx",
    "1,4-butanediol production (non-KEGG metabolite)": "BDO_production_pathways.xlsx",
}

CURRENCY_KEGG = {"C00002": "ATP", "C00008": "ADP", "C00009": "Pi",
                 "C00001": "H2O", "C00080": "H+"}


# --------------------------------------------------------------------------- #
# Serialization helpers (dataframes <-> dcc.Store)
# --------------------------------------------------------------------------- #
def df_to_records(df):
    return df.to_dict("records")


def records_to_df(records, cols):
    df = pd.DataFrame(records or [])
    for c in cols:
        if c not in df.columns:
            df[c] = pd.Series(dtype=object)
    return df[cols]


def load_example(filename):
    path = os.path.join(EXAMPLE_DIR, filename)
    df_met, df_rxn = io.load_excel(path)
    config = io.load_config(path) or {}
    return df_met, df_rxn, config


# --------------------------------------------------------------------------- #
# Layout
# --------------------------------------------------------------------------- #
def _table(table_id, columns, editable=True, data=None):
    return dash_table.DataTable(
        id=table_id,
        columns=[{"name": c, "id": c} for c in columns],
        data=data or [],
        editable=editable,
        row_deletable=editable,
        page_size=50,
        style_table={"overflowX": "auto", "maxHeight": "420px",
                     "overflowY": "auto"},
        style_cell={"fontFamily": "monospace", "fontSize": "13px",
                    "textAlign": "left", "padding": "4px"},
        style_header={"fontWeight": "bold"},
    )


# --------------------------------------------------------------------------- #
# Manual reaction entry: dynamic substrate/product rows (dcc components, not
# DataTable cells — the DataTable's cell editor swallows Backspace/Delete/
# Home/arrow keys for its own cell-to-cell navigation, so composing a new
# reaction character-by-character in a table cell is effectively unfixable
# without replacing the editor).  Each row gets a stable, unique pattern-
# matching id ({"type": ..., "side": "subs"/"prods", "index": n}) so Dash can
# address an arbitrary number of rows without per-row callbacks.
# --------------------------------------------------------------------------- #
def _met_options(search_value, known_ids):
    """Dropdown options for a metabolite picker: existing ids filtered by
    ``search_value`` (substring match), plus a synthetic "new metabolite"
    entry when the typed text doesn't match an existing id (case-insensitive)
    — picking it uses the typed text itself as a brand-new working id."""
    known_ids = sorted(set(known_ids))
    known_lower = {i.lower() for i in known_ids}
    q = (search_value or "").strip()
    if not q:
        return [{"label": i, "value": i} for i in known_ids]
    matches = [i for i in known_ids if q.lower() in i.lower()]
    opts = [{"label": i, "value": i} for i in matches]
    if q.lower() not in known_lower:
        opts.insert(0, {"label": f"＋ New metabolite: {q}", "value": q})
    return opts


def _metab_row(side, idx, known_ids):
    return html.Div([
        dbc.Row([
            dbc.Col(dbc.Input(id={"type": "man-coeff", "side": side, "index": idx},
                              type="number", value=1, min=0, step="any", size="sm"),
                    width=3),
            dbc.Col(dcc.Dropdown(id={"type": "man-met", "side": side, "index": idx},
                                 options=_met_options(None, known_ids),
                                 placeholder="metabolite…", clearable=False,
                                 style={"fontSize": "13px"}),
                    width=7),
            dbc.Col(dbc.Button("×", id={"type": "man-rm", "side": side, "index": idx},
                               size="sm", color="danger", outline=True,
                               style={"padding": "0 8px"}),
                    width=2),
        ], className="g-1 mb-1"),
    ], id={"type": "man-row", "side": side, "index": idx})


def _new_met_panel(new_met_ids):
    """A follow-up mini-form for freshly-created placeholder metabolites: a
    KEGG id / chemical formula field per metabolite, so balance can resolve
    them right away instead of only via the Metabolites table below."""
    if not new_met_ids:
        return None
    rows = [
        dbc.Row([
            dbc.Col(html.Code(mid), width=3, className="pt-2"),
            dbc.Col(dbc.Input(id={"type": "new-met-kegg", "mid": mid},
                              placeholder="KEGG id (e.g. C00002)", size="sm"),
                    width=4),
            dbc.Col(dbc.Input(id={"type": "new-met-formula", "mid": mid},
                              placeholder="Chemical formula (e.g. C6H12O6)", size="sm"),
                    width=5),
        ], className="g-2 mb-1")
        for mid in new_met_ids
    ]
    return dbc.Card(dbc.CardBody([
        html.Div("New metabolite(s) — add a KEGG id or chemical formula so "
                 "they can balance (either is enough; both is fine too):",
                 className="small fw-bold mb-2"),
        *rows,
        dbc.Button("Save metabolite details", id="btn-save-new-mets",
                   size="sm", color="primary", className="mt-1"),
    ]), color="info", outline=True, className="mt-2")


def _next_row_index(children):
    if not children:
        return 0
    idxs = [c["props"]["id"]["index"] for c in children]
    return max(idxs) + 1


def workbench_tab():
    # Seed the tables with the default example directly in the layout so they
    # are never empty on first paint.  (A DataTable mounted with no data can
    # render blank until a reflow — e.g. switching tabs — forces a redraw.)
    df_met0, df_rxn0, _cfg0 = load_example("Example1_EMP_lactate.xlsx")
    met0, rxn0 = df_to_records(df_met0), df_to_records(df_rxn0)
    # Seed the selector options *and* values from the default example too, so the
    # run-time selectors (incl. currency metabolites) are populated on first
    # paint.  A Dropdown value set while its options are still empty gets dropped,
    # which is why the currency metabolites showed up blank on load.
    ex0 = [{"label": r, "value": r} for r in df_rxn0["ID"].dropna().astype(str)
           if r.startswith("EX")]
    metopt0 = [{"label": m, "value": m} for m in df_met0["ID"].dropna().astype(str)]
    pseudo0 = sorted(_cfg0.get("PSEUDO_METS") or [])
    met_ids0 = df_met0["ID"].dropna().astype(str).tolist()
    return dbc.Container([
        dbc.Row([
            dbc.Col([
                html.H5("1 · Load / save"),
                dbc.InputGroup([
                    dbc.Select(
                        id="example-select",
                        options=[{"label": k, "value": v} for k, v in EXAMPLES.items()],
                        value="Example1_EMP_lactate.xlsx"),
                    dbc.Button("Load example", id="btn-load-example", color="secondary"),
                ], className="mb-2"),
                dcc.Upload(
                    id="upload-excel",
                    children=html.Div(["⬆ Drag/drop or ", html.A("select an .xlsx")]),
                    style={"border": "1px dashed #aaa", "borderRadius": "6px",
                           "padding": "10px", "textAlign": "center"},
                    className="mb-2"),
                dbc.Button("⬇ Download workbook (.xlsx)", id="btn-download-xlsx",
                           color="secondary", outline=True, className="mb-2 w-100"),
                dcc.Download(id="download-xlsx"),
                dcc.Loading(html.Div(id="io-status",
                                     style={"maxHeight": "120px", "overflowY": "auto",
                                            "fontSize": "12px", "fontFamily": "monospace"})),
            ], md=6),
        ], className="mb-3"),

        html.Hr(),
        html.H5("Add a reaction manually"),
        html.Div("Pick substrates and products from existing metabolites, or "
                 "type a working id that doesn't exist yet and choose "
                 "\"＋ New metabolite: …\" to create it as a placeholder (fill "
                 "in its KEGG id or chemical formula afterwards so it can "
                 "balance). Add rows for extra reactants/products as needed.",
                 className="text-muted small mb-2"),
        dbc.Row([
            dbc.Col(dbc.Input(id="man-rxn-id", placeholder="Reaction id (e.g. PTS_GLC)"), md=3),
            dbc.Col(dbc.Input(id="man-rxn-name", placeholder="Name (optional)"), md=3),
            dbc.Col(dbc.Select(id="man-rxn-rev", value="0",
                               options=[{"label": "irreversible", "value": "0"},
                                        {"label": "reversible", "value": "1"}]),
                    md=2),
            dbc.Col(dbc.Button("Add reaction", id="btn-add-manual-rxn",
                               color="primary", className="w-100"), md=4),
        ], className="g-2 mb-2"),
        dbc.Row([
            dbc.Col([
                html.Div("Substrates", className="fw-bold small mb-1"),
                html.Div(id="man-subs-rows", children=[_metab_row("subs", 0, met_ids0)]),
                dbc.Button("+ Add metabolite", id="btn-man-add-subs", size="sm",
                           color="light"),
            ], md=6),
            dbc.Col([
                html.Div("Products", className="fw-bold small mb-1"),
                html.Div(id="man-prods-rows", children=[_metab_row("prods", 0, met_ids0)]),
                dbc.Button("+ Add metabolite", id="btn-man-add-prods", size="sm",
                           color="light"),
            ], md=6),
        ], className="mb-2"),
        html.Div(id="man-messages", style={"fontSize": "12px", "marginTop": "0.5rem"}),
        html.Div(id="man-new-met-panel"),

        html.Hr(),
        dbc.Row([
            dbc.Col([
                html.H5("Reactions"),
                dbc.Button("+ Add reaction row", id="btn-add-rxn", size="sm",
                           color="light", className="mb-1"),
                _table("table-reactions", io.REACTION_COLS, data=rxn0),
            ], md=7),
            dbc.Col([
                html.H5("Metabolites"),
                html.Div(id="met-missing-warning"),
                dbc.Button("+ Add metabolite row", id="btn-add-met", size="sm",
                           color="light", className="mb-1"),
                _table("table-metabolites", io.METABOLITE_COLS, data=met0),
            ], md=5),
        ], className="mb-3"),

        html.Hr(),
        html.H5("2 · Live balance check"),
        html.Div("Exchange reactions (EX…) are expected to be unbalanced and are "
                 "not flagged as errors.", className="text-muted small mb-2"),
        dcc.Loading(html.Div(id="balance-panel")),
        dbc.Card(dbc.CardBody([
            html.H6("One-click fixes", className="mb-2"),
            dbc.InputGroup([
                dbc.InputGroupText("Reaction"),
                dbc.Col(dcc.Dropdown(id="fix-reaction", options=[],
                                     placeholder="an imbalanced reaction…"),
                        className="p-0"),
            ], className="mb-2"),
            dbc.ButtonGroup([
                dbc.Button("+ H⁺ left", id="fix-h-left", size="sm", outline=True, color="primary"),
                dbc.Button("+ H⁺ right", id="fix-h-right", size="sm", outline=True, color="primary"),
                dbc.Button("+ H₂O left", id="fix-w-left", size="sm", outline=True, color="info"),
                dbc.Button("+ H₂O right", id="fix-w-right", size="sm", outline=True, color="info"),
            ]),
        ]), className="mb-3"),

        html.Hr(),
        html.H5("3 · Run-time selectors"),
        html.Div("Scientific choices (reversibility/direction, which metabolites "
                 "are external) stay with you — these are populated from the "
                 "exchanges present. Every exchange you want to keep must appear "
                 "in Substrates, Products, or Freely-reversible below — any "
                 "exchange in none of the three is deleted from the model.",
                 className="text-muted small mb-2"),
        dbc.Row([
            dbc.Col([dbc.Label("Model name"),
                     dbc.Input(id="cfg-model-name", value="model")], md=4),
            dbc.Col([dbc.Label("Substrate(s)"),
                     html.Div("Exchanges with uptake enabled — your carbon "
                              "source(s); put the primary one first, since only "
                              "the first entry sets EFM normalisation and the "
                              "yield basis, and all of them are blocked together "
                              "in the energy-cycle check.",
                              className="text-muted small mb-1"),
                     dcc.Dropdown(id="cfg-substrates", multi=True, options=ex0,
                                  value=_cfg0.get("SUBSTRATES", []))], md=4),
            dbc.Col([dbc.Label("Energy product (ΔG cycle check)"),
                     html.Div("The exchange (usually ATP) maximised with all "
                              "substrate uptake blocked, to catch "
                              "energy-generating cycles — not a co-substrate.",
                              className="text-muted small mb-1"),
                     dcc.Dropdown(id="cfg-energy", options=ex0,
                                  value=_cfg0.get("ENERGY_PRODUCT"))], md=4),
        ], className="mb-2"),
        dbc.Row([
            dbc.Col([dbc.Label("Products (secretable)"),
                     html.Div("Exchanges with uptake blocked but secretion "
                              "allowed — your pathway's outputs (e.g. lactate, "
                              "CO₂, any net-produced ATP/ADP/Pi).",
                              className="text-muted small mb-1"),
                     dcc.Dropdown(id="cfg-products", multi=True, options=ex0,
                                  value=_cfg0.get("PRODUCTS", []))], md=4),
            dbc.Col([dbc.Label("Carbon products (yield check)"),
                     html.Div("A subset of Products: the carbon-containing end "
                              "product(s) to report max theoretical FBA yield "
                              "for — not currency carriers like ATP/ADP/Pi.",
                              className="text-muted small mb-1"),
                     dcc.Dropdown(id="cfg-carbon", multi=True, options=ex0,
                                  value=_cfg0.get("CARBON_PRODUCTS", []))], md=4),
            dbc.Col([dbc.Label("Freely-reversible exchanges"),
                     html.Div("Exchanges kept reversible (uptake and secretion "
                              "both allowed) — cofactor/currency exchanges "
                              "(ADP, Pi, H₂O, H⁺, …) whose net direction isn't "
                              "fixed a priori; the EFMs decide the sign.",
                              className="text-muted small mb-1"),
                     dcc.Dropdown(id="cfg-rev", multi=True, options=ex0,
                                  value=_cfg0.get("REV_ALLOWED", []))], md=4),
        ], className="mb-2"),
        dbc.Row([
            dbc.Col([dbc.Label("Currency / unbalanced metabolites — excluded from Ω flux sum"),
                     html.Div("Metabolites net-produced/-consumed alongside the "
                              "main carbon conversion (ATP, ADP, Pi, H₂O, H⁺, …).",
                              className="text-muted small mb-1"),
                     dcc.Dropdown(id="cfg-pseudo", multi=True, options=metopt0,
                                  value=pseudo0)], md=12),
        ], className="mb-3"),
    ], fluid=True)


def sanity_tab():
    return dbc.Container([
        html.P("Before the (potentially expensive) EFM enumeration, these checks "
               "catch the most common modelling mistakes — reactions pointing the "
               "wrong way, products that can't actually be made, and dead-end "
               "reactions. Each card explains what it means and what to do if it "
               "fails.", className="text-muted mt-3"),
        dbc.Button("Run sanity checks", id="btn-fba", color="primary",
                   className="mb-3"),
        dcc.Loading(html.Div(id="fba-output")),
    ], fluid=True)


def efm_tab():
    return dbc.Container([
        dbc.Alert("efmtool time scales super-exponentially with network size. "
                  "Networks with bypasses can explode (the paper notes "
                  "ecolicore → ~272M EFMs). Use the cap and Cancel if needed.",
                  color="warning", className="my-2"),
        dbc.Row([
            dbc.Col(dbc.Button("Enumerate EFMs", id="btn-efm", color="primary"), width="auto"),
            dbc.Col(dbc.Button("Cancel", id="btn-efm-cancel", color="danger",
                               outline=True), width="auto"),
            dbc.Col([dbc.InputGroupText("EFM cap"),
                     ], width="auto"),
            dbc.Col(dbc.Input(id="efm-cap", type="number", value=100000,
                              style={"width": "140px"}), width="auto"),
        ], className="g-2 align-items-center mb-2"),
        html.Div(id="efm-progress", className="text-muted"),
        dcc.Loading(html.Div(id="efm-output")),
    ], fluid=True)


def thermo_tab():
    return dbc.Container([
        dbc.Row([
            dbc.Col([dbc.Label("ΔG driving Ω"),
                     dbc.Select(id="thermo-which",
                                options=[{"label": "ΔGm′ (physiological)", "value": "dGm"},
                                         {"label": "ΔG°′ (standard)", "value": "dG0prime"}],
                                value="dGm")], md=3),
            dbc.Col(dbc.Button("Compute net conversions & ΔG", id="btn-thermo",
                               color="primary", className="mt-4"), width="auto"),
        ], className="mb-2"),
        dbc.Alert("eQuilibrator downloads a large thermodynamic cache on first "
                  "use only (minutes, needs network). Ω here equals the MDF with "
                  "no concentration bounds.", color="info", className="small"),
        dcc.Loading(html.Div(id="thermo-output")),
    ], fluid=True)


def downloads_tab():
    return dbc.Container([
        html.H5("Downloads", className="my-3"),
        dbc.ButtonGroup([
            dbc.Button("COBRApy JSON", id="btn-dl-json", color="secondary"),
            dbc.Button("SBML", id="btn-dl-sbml", color="secondary"),
            dbc.Button("Results (XLSX)", id="btn-dl-results", color="secondary"),
        ]),
        dcc.Download(id="download-json"),
        dcc.Download(id="download-sbml"),
        dcc.Download(id="download-results"),
        html.Div(id="dl-status", className="mt-2 text-muted"),
    ], fluid=True)


def _startup_banner():
    if not efm.java_available():
        return ("⚠ Java (JVM) not detected — the EFM step will fail. "
                "Install a JRE (e.g. `sudo apt install default-jre`).")
    return "JVM detected ✓  — ready for EFM enumeration."


app.layout = dbc.Container([
    dcc.Store(id="store-metabolites"),
    dcc.Store(id="store-reactions"),
    dcc.Store(id="store-efm"),       # serialized normalized EFM table + meta
    dcc.Store(id="store-thermo"),    # serialized thermo results

    html.H3("Catabolic net conversions, EFMs & ΔG", className="mt-3"),
    html.P("A computational tool: assemble a pathway, check balance, enumerate "
           "elementary flux modes, and compute net-conversion ΔG°′/ΔGm′ and the "
           "Ω measure (kJ/mol).", className="text-muted"),
    dbc.Tabs([
        dbc.Tab(workbench_tab(), label="① Workbench", tab_id="tab-workbench"),
        dbc.Tab(sanity_tab(), label="② Sanity checks", tab_id="tab-fba"),
        dbc.Tab(efm_tab(), label="③ EFMs", tab_id="tab-efm"),
        dbc.Tab(thermo_tab(), label="④ Net conversions & ΔG", tab_id="tab-thermo"),
        dbc.Tab(downloads_tab(), label="⑤ Downloads", tab_id="tab-downloads"),
    ], id="main-tabs", active_tab="tab-workbench"),
    html.Footer(html.Small(_startup_banner()), className="text-muted my-3"),
], fluid=True)


# --------------------------------------------------------------------------- #
# Workbench callbacks: load, KEGG fetch, edits → tables + balance
# --------------------------------------------------------------------------- #
@app.callback(
    Output("store-metabolites", "data"),
    Output("store-reactions", "data"),
    Output("table-metabolites", "data"),
    Output("table-reactions", "data"),
    Output("cfg-model-name", "value"),
    Output("io-status", "children"),
    Input("btn-load-example", "n_clicks"),
    Input("upload-excel", "contents"),
    Input("btn-add-rxn", "n_clicks"),
    Input("btn-add-met", "n_clicks"),
    State("example-select", "value"),
    State("table-metabolites", "data"),
    State("table-reactions", "data"),
    State("cfg-model-name", "value"),
    prevent_initial_call=False,
)
def update_tables(_n_ex, upload, _n_arxn, _n_amet, example_file,
                  met_data, rxn_data, model_name):
    trig = ctx.triggered_id
    status = no_update

    # First load (no trigger) → default example.
    if trig is None:
        df_met, df_rxn, cfg = load_example("Example1_EMP_lactate.xlsx")
        return (df_to_records(df_met), df_to_records(df_rxn),
                df_to_records(df_met), df_to_records(df_rxn),
                cfg.get("MODEL_NAME", "model"), "")

    if trig == "btn-load-example":
        df_met, df_rxn, cfg = load_example(example_file or "Example1_EMP_lactate.xlsx")
        return (df_to_records(df_met), df_to_records(df_rxn),
                df_to_records(df_met), df_to_records(df_rxn),
                cfg.get("MODEL_NAME", model_name), f"Loaded {example_file}")

    if trig == "upload-excel" and upload:
        _header, b64 = upload.split(",", 1)
        data = base64.b64decode(b64)
        try:
            df_met, df_rxn = io.load_excel(data)
            cfg = io.load_config(data) or {}
        except Exception as e:
            return (no_update, no_update, no_update, no_update, no_update,
                    f"⚠ Could not read workbook: {e}")
        return (df_to_records(df_met), df_to_records(df_rxn),
                df_to_records(df_met), df_to_records(df_rxn),
                cfg.get("MODEL_NAME", model_name), "Workbook uploaded ✓")

    # Work from whatever is currently in the tables.
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)

    if trig == "btn-add-rxn":
        df_rxn = pd.concat([df_rxn, pd.DataFrame([{c: "" for c in io.REACTION_COLS}])],
                           ignore_index=True)
    elif trig == "btn-add-met":
        df_met = pd.concat([df_met, pd.DataFrame([{c: "" for c in io.METABOLITE_COLS}])],
                           ignore_index=True)

    return (df_to_records(df_met), df_to_records(df_rxn),
            df_to_records(df_met), df_to_records(df_rxn), no_update, status)


# Sync edited DataTables back into the stores (so balance + downstream see edits).
@app.callback(
    Output("store-metabolites", "data", allow_duplicate=True),
    Input("table-metabolites", "data"),
    prevent_initial_call=True,
)
def sync_metabolites(data):
    return data


@app.callback(
    Output("store-reactions", "data", allow_duplicate=True),
    Input("table-reactions", "data"),
    prevent_initial_call=True,
)
def sync_reactions(data):
    return data


# --------------------------------------------------------------------------- #
# Manual reaction entry: add/remove substrate & product rows, live metabolite
# search, and the final "Add reaction" submit.
# --------------------------------------------------------------------------- #
@app.callback(
    Output("man-subs-rows", "children"),
    Input("btn-man-add-subs", "n_clicks"),
    State("man-subs-rows", "children"),
    State("store-metabolites", "data"),
    prevent_initial_call=True,
)
def add_subs_row(_n, children, met_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    known_ids = df_met["ID"].dropna().astype(str).tolist()
    children = children or []
    return children + [_metab_row("subs", _next_row_index(children), known_ids)]


@app.callback(
    Output("man-prods-rows", "children"),
    Input("btn-man-add-prods", "n_clicks"),
    State("man-prods-rows", "children"),
    State("store-metabolites", "data"),
    prevent_initial_call=True,
)
def add_prods_row(_n, children, met_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    known_ids = df_met["ID"].dropna().astype(str).tolist()
    children = children or []
    return children + [_metab_row("prods", _next_row_index(children), known_ids)]


@app.callback(
    Output("man-subs-rows", "children", allow_duplicate=True),
    Input({"type": "man-rm", "side": "subs", "index": ALL}, "n_clicks"),
    State("man-subs-rows", "children"),
    prevent_initial_call=True,
)
def remove_subs_row(_clicks, children):
    if not children or len(children) <= 1 or not any(_clicks):
        return no_update
    trig = ctx.triggered_id
    return [c for c in children if c["props"]["id"]["index"] != trig["index"]]


@app.callback(
    Output("man-prods-rows", "children", allow_duplicate=True),
    Input({"type": "man-rm", "side": "prods", "index": ALL}, "n_clicks"),
    State("man-prods-rows", "children"),
    prevent_initial_call=True,
)
def remove_prods_row(_clicks, children):
    if not children or len(children) <= 1 or not any(_clicks):
        return no_update
    trig = ctx.triggered_id
    return [c for c in children if c["props"]["id"]["index"] != trig["index"]]


@app.callback(
    Output({"type": "man-met", "side": MATCH, "index": MATCH}, "options"),
    Input({"type": "man-met", "side": MATCH, "index": MATCH}, "search_value"),
    Input({"type": "man-met", "side": MATCH, "index": MATCH}, "value"),
    State("store-metabolites", "data"),
)
def update_man_met_options(search_value, current_value, met_data):
    # Re-runs both when the user types (search_value) and right after they
    # pick something (value) — a freshly-picked "new metabolite" doesn't
    # exist in store-metabolites yet (it's only created on submit), so
    # without this the very next options refresh (e.g. the search box
    # closing) would drop it from the list and the dropdown would lose its
    # selection silently.
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    known_ids = df_met["ID"].dropna().astype(str).tolist()
    opts = _met_options(search_value, known_ids)
    if (current_value and current_value not in known_ids
            and not any(o["value"] == current_value for o in opts)):
        opts.append({"label": f"＋ New metabolite: {current_value}",
                     "value": current_value})
    return opts


def _side_terms(coeffs, mets):
    terms = []
    for c, m in zip(coeffs or [], mets or []):
        if not m:
            continue
        try:
            c = float(c) if c not in (None, "") else 1.0
        except (TypeError, ValueError):
            c = 1.0
        terms.append((c, str(m).strip()))
    return terms


def _terms_to_str(terms):
    parts = []
    for c, m in terms:
        parts.append(m if abs(c - 1.0) < 1e-9 else f"{c:g} {m}")
    return " + ".join(parts) if parts else "0"


@app.callback(
    Output("store-metabolites", "data", allow_duplicate=True),
    Output("store-reactions", "data", allow_duplicate=True),
    Output("table-metabolites", "data", allow_duplicate=True),
    Output("table-reactions", "data", allow_duplicate=True),
    Output("man-rxn-id", "value"),
    Output("man-rxn-name", "value"),
    Output("man-subs-rows", "children", allow_duplicate=True),
    Output("man-prods-rows", "children", allow_duplicate=True),
    Output("man-messages", "children"),
    Output("man-new-met-panel", "children"),
    Input("btn-add-manual-rxn", "n_clicks"),
    State("man-rxn-id", "value"),
    State("man-rxn-name", "value"),
    State("man-rxn-rev", "value"),
    State({"type": "man-coeff", "side": "subs", "index": ALL}, "value"),
    State({"type": "man-met", "side": "subs", "index": ALL}, "value"),
    State({"type": "man-coeff", "side": "prods", "index": ALL}, "value"),
    State({"type": "man-met", "side": "prods", "index": ALL}, "value"),
    State("store-metabolites", "data"),
    State("store-reactions", "data"),
    prevent_initial_call=True,
)
def add_manual_reaction(_n, rid, name, rev, subs_c, subs_m, prods_c, prods_m,
                        met_data, rxn_data):
    rid = (rid or "").strip()
    if not rid:
        return (no_update,) * 4 + (no_update, no_update, no_update, no_update,
                                    dbc.Alert("Enter a reaction id.", color="warning"),
                                    no_update)

    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if (df_rxn["ID"].astype(str).str.strip() == rid).any():
        return (no_update,) * 4 + (no_update, no_update, no_update, no_update,
                                    dbc.Alert(f"⚠ Reaction '{rid}' already exists.",
                                              color="warning"),
                                    no_update)

    subs_terms = _side_terms(subs_c, subs_m)
    prods_terms = _side_terms(prods_c, prods_m)
    if not subs_terms or not prods_terms:
        return (no_update,) * 4 + (no_update, no_update, no_update, no_update,
                                    dbc.Alert("Pick at least one substrate and one "
                                              "product.", color="warning"),
                                    no_update)

    eq = f"{_terms_to_str(subs_terms)} <=> {_terms_to_str(prods_terms)}"
    df_rxn = pd.concat([df_rxn, pd.DataFrame([{
        "ID": rid, "Name": (name or "").strip() or rid,
        "Reaction stoichiometry": eq,
        "Reversibility": int(rev) if rev in ("0", "1") else 0,
    }])], ignore_index=True)

    known = {str(x).strip() for x in df_met["ID"].dropna()}
    new_mets = []
    for _c, mid in subs_terms + prods_terms:
        if mid and mid not in known:
            new_mets.append({"ID": mid, "Name": mid, "KEGG ID": "",
                             "Chemical formula": ""})
            known.add(mid)
    if new_mets:
        df_met = pd.concat([df_met, pd.DataFrame(new_mets)], ignore_index=True)

    parts = [f"✓ Added reaction {rid}: {eq}"]
    if new_mets:
        parts.append("new metabolite placeholder(s): "
                     + ", ".join(m["ID"] for m in new_mets)
                     + " — fill in a KEGG id or chemical formula below so "
                       "they balance.")

    fresh_known = df_met["ID"].dropna().astype(str).tolist()
    return (df_to_records(df_met), df_to_records(df_rxn),
            df_to_records(df_met), df_to_records(df_rxn),
            "", "",
            [_metab_row("subs", 0, fresh_known)],
            [_metab_row("prods", 0, fresh_known)],
            dbc.Alert(" — ".join(parts), color="success"),
            _new_met_panel([m["ID"] for m in new_mets]))


@app.callback(
    Output("store-metabolites", "data", allow_duplicate=True),
    Output("table-metabolites", "data", allow_duplicate=True),
    Output("man-new-met-panel", "children", allow_duplicate=True),
    Input("btn-save-new-mets", "n_clicks"),
    State({"type": "new-met-kegg", "mid": ALL}, "value"),
    State({"type": "new-met-kegg", "mid": ALL}, "id"),
    State({"type": "new-met-formula", "mid": ALL}, "value"),
    State("store-metabolites", "data"),
    prevent_initial_call=True,
)
def save_new_met_details(_n, kegg_vals, kegg_ids, formula_vals, met_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    for id_dict, kegg_val, formula_val in zip(kegg_ids, kegg_vals, formula_vals):
        mid = id_dict["mid"]
        mask = df_met["ID"].astype(str) == mid
        if kegg_val and kegg_val.strip():
            df_met.loc[mask, "KEGG ID"] = kegg_val.strip()
        if formula_val and formula_val.strip():
            df_met.loc[mask, "Chemical formula"] = formula_val.strip()
    return df_to_records(df_met), df_to_records(df_met), None


def _is_blank(v):
    return v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() == ""


@app.callback(
    Output("met-missing-warning", "children"),
    Input("store-metabolites", "data"),
)
def check_missing_met_info(met_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_met = df_met[df_met["ID"].astype(str).str.strip() != ""]
    missing = [str(r["ID"]) for _, r in df_met.iterrows()
              if _is_blank(r.get("KEGG ID")) and _is_blank(r.get("Chemical formula"))]
    if not missing:
        return None
    return dbc.Alert([
        html.B(f"{len(missing)} metabolite(s) "),
        "have neither a KEGG id nor a chemical formula — balance can't be "
        "checked for these: ",
        html.Code(", ".join(missing)),
    ], color="warning", className="small py-2 mb-1")


# Live balance panel + populate selectors from the current tables.
@app.callback(
    Output("balance-panel", "children"),
    Output("fix-reaction", "options"),
    Output("fix-reaction", "value"),
    Output("cfg-substrates", "options"),
    Output("cfg-products", "options"),
    Output("cfg-carbon", "options"),
    Output("cfg-energy", "options"),
    Output("cfg-rev", "options"),
    Output("cfg-pseudo", "options"),
    Input("store-metabolites", "data"),
    Input("store-reactions", "data"),
    State("fix-reaction", "value"),
)
def update_balance(met_data, rxn_data, fix_value):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if df_rxn.empty:
        return (html.Div("No reactions yet."), [], None, [], [], [], [], [], [])

    df_bal, ok, errors = run.run_balance(df_met, df_rxn)

    def status_cell(row):
        if row["atom_balanced"] is True and row["charge_balanced"] is True:
            return "✓ balanced"
        if row.get("is_exchange"):
            return "exchange (unbalanced ok)"
        if row["atom_balanced"] is None:
            return f"⚠ {row['notes']}"
        return f"✗ {row['atom_imbalance'] or row['notes']}"

    disp = pd.DataFrame({
        "Reaction": df_bal["reaction_id"],
        "Status": df_bal.apply(status_cell, axis=1),
        "KEGG equation": df_bal["equation_kegg"],
    })
    table = dash_table.DataTable(
        data=disp.to_dict("records"),
        columns=[{"name": c, "id": c} for c in disp.columns],
        style_cell={"fontFamily": "monospace", "fontSize": "12px",
                    "textAlign": "left", "padding": "4px"},
        style_data_conditional=[
            {"if": {"filter_query": '{Status} contains "✗"'},
             "backgroundColor": "#f8d7da"},
            {"if": {"filter_query": '{Status} contains "⚠"'},
             "backgroundColor": "#fff3cd"},
            {"if": {"filter_query": '{Status} contains "✓"'},
             "backgroundColor": "#d4edda"},
        ],
        page_size=60,
        style_table={"overflowX": "auto", "maxHeight": "360px", "overflowY": "auto"},
    )
    if ok:
        banner = dbc.Alert("✓ All non-exchange reactions are atom- and "
                           "charge-balanced — ready to model.", color="success")
    else:
        items = [html.Li(f"{e['reaction_id']}: {e['imbalance'] or e['notes']}")
                 for e in errors]
        banner = dbc.Alert([html.B(f"{len(errors)} imbalanced reaction(s):"),
                            html.Ul(items)], color="danger")

    # selector option sets
    exchanges = [r for r in df_rxn["ID"].dropna().astype(str) if r.startswith("EX")]
    ex_opts = [{"label": r, "value": r} for r in exchanges]
    met_ids = [str(m) for m in df_met["ID"].dropna()]
    met_opts = [{"label": m, "value": m} for m in met_ids]
    flagged = [{"label": e["reaction_id"], "value": e["reaction_id"]} for e in errors]

    # Keep the current pick if it's still imbalanced; otherwise advance to the
    # first remaining flagged reaction (so the fix buttons always act on a real,
    # currently-selected reaction — never a stale one).
    flagged_ids = [e["reaction_id"] for e in errors]
    fix_val = fix_value if fix_value in flagged_ids else (
        flagged_ids[0] if flagged_ids else None)

    return (html.Div([banner, table]), flagged, fix_val,
            ex_opts, ex_opts, ex_opts, ex_opts, ex_opts, met_opts)


# Default selector *values* whenever a workbook is (re)loaded.
@app.callback(
    Output("cfg-substrates", "value"),
    Output("cfg-products", "value"),
    Output("cfg-carbon", "value"),
    Output("cfg-energy", "value"),
    Output("cfg-rev", "value"),
    Output("cfg-pseudo", "value"),
    Input("btn-load-example", "n_clicks"),
    Input("upload-excel", "contents"),
    State("example-select", "value"),
    State("store-metabolites", "data"),
    prevent_initial_call=False,
)
def default_selectors(_n, upload, example_file, met_data):
    trig = ctx.triggered_id
    cfg = {}
    # Read the metabolite table from the SAME source we're loading the config
    # from — not from the store, which lags a paint behind on first load and
    # would leave the currency-metabolite suggestion empty.
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    if trig == "upload-excel" and upload:
        try:
            _h, b64 = upload.split(",", 1)
            data = base64.b64decode(b64)
            cfg = io.load_config(data) or {}
            df_met, _r = io.load_excel(data)
        except Exception:
            cfg = {}
    else:
        df_met, _r, cfg = load_example(example_file or "Example1_EMP_lactate.xlsx")

    # Suggest currency pseudo-mets present in the table.
    suggested_pseudo = cfg.get("PSEUDO_METS")
    if not suggested_pseudo and not df_met.empty:
        suggested_pseudo = sorted({
            str(r["ID"]) for _, r in df_met.iterrows()
            if str(r.get("KEGG ID")) in CURRENCY_KEGG})
    return (cfg.get("SUBSTRATES", []), cfg.get("PRODUCTS", []),
            cfg.get("CARBON_PRODUCTS", []), cfg.get("ENERGY_PRODUCT"),
            cfg.get("REV_ALLOWED", []), sorted(suggested_pseudo or []))


# One-click H+/H2O fixes on a flagged reaction.
@app.callback(
    Output("store-reactions", "data", allow_duplicate=True),
    Output("table-reactions", "data", allow_duplicate=True),
    Input("fix-h-left", "n_clicks"),
    Input("fix-h-right", "n_clicks"),
    Input("fix-w-left", "n_clicks"),
    Input("fix-w-right", "n_clicks"),
    State("fix-reaction", "value"),
    State("store-reactions", "data"),
    State("store-metabolites", "data"),
    prevent_initial_call=True,
)
def apply_fix(_hl, _hr, _wl, _wr, rxn_id, rxn_data, met_data):
    if not rxn_id:
        return no_update, no_update
    trig = ctx.triggered_id
    # Resolve the working ID for H+ / H2O from the metabolite table (KEGG-keyed),
    # preferring the intracellular form (IDs not ending in "ex").
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    kegg_to_id = {}
    for _, r in df_met.iterrows():
        k, mid = str(r["KEGG ID"]), str(r["ID"])
        if k not in kegg_to_id or kegg_to_id[k].endswith("ex"):
            kegg_to_id[k] = mid
    species = (kegg_to_id.get("C00080", "H") if "h-" in trig
               else kegg_to_id.get("C00001", "H2O"))
    side = "left" if trig.endswith("left") else "right"

    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    mask = df_rxn["ID"].astype(str) == str(rxn_id)
    if not mask.any():
        return no_update, no_update
    idx = df_rxn.index[mask][0]
    eq = str(df_rxn.loc[idx, "Reaction stoichiometry"])
    eq_norm = eq
    for arrow in ["<=>", "->", "=>", "→", "<->"]:
        if arrow in eq:
            lhs, rhs = eq.split(arrow, 1)
            if side == "left":
                lhs = f"{lhs.strip()} + {species}"
            else:
                rhs = f"{rhs.strip()} + {species}"
            eq_norm = f"{lhs.strip()} {arrow} {rhs.strip()}"
            break
    df_rxn.loc[idx, "Reaction stoichiometry"] = eq_norm
    recs = df_to_records(df_rxn)
    return recs, recs


# Download current workbook.
@app.callback(
    Output("download-xlsx", "data"),
    Input("btn-download-xlsx", "n_clicks"),
    State("store-metabolites", "data"),
    State("store-reactions", "data"),
    State("cfg-model-name", "value"),
    State("cfg-substrates", "value"),
    State("cfg-products", "value"),
    State("cfg-carbon", "value"),
    State("cfg-energy", "value"),
    State("cfg-rev", "value"),
    State("cfg-pseudo", "value"),
    prevent_initial_call=True,
)
def download_workbook(_n, met_data, rxn_data, name, subs, prods, carbon, energy,
                      rev, pseudo):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    config = _collect_config(name, subs, prods, carbon, energy, rev, pseudo)
    data = io.save_excel(df_met, df_rxn, config=config)
    return dcc.send_bytes(data, f"{name or 'workbook'}.xlsx")


def _collect_config(name, subs, prods, carbon, energy, rev, pseudo):
    return {
        "MODEL_NAME": name or "model",
        "SUBSTRATES": subs or [],
        "PRODUCTS": prods or [],
        "CARBON_PRODUCTS": carbon or [],
        "ENERGY_PRODUCT": energy,
        "REV_ALLOWED": rev or [],
        "PSEUDO_METS": set(pseudo or []),
    }


# --------------------------------------------------------------------------- #
# Model & FBA
# --------------------------------------------------------------------------- #
# A sanity-check card: coloured badge + plain-language "what it means" +
# the numeric result + an interpretation of what the result implies.
_STATUS = {
    "pass": ("success", "✓ Pass"),
    "warn": ("warning", "⚠ Check this"),
    "fail": ("danger", "✗ Problem"),
    "info": ("secondary", "ℹ For info"),
}


def _ex_name(rid):
    return rid[2:] if isinstance(rid, str) and rid.startswith("EX") else rid


def _check_card(status, title, what, result, implication):
    color, badge = _STATUS[status]
    return dbc.Card(dbc.CardBody([
        html.Div([dbc.Badge(badge, color=color, className="me-2"),
                  html.Span(title, className="fw-bold")], className="mb-2"),
        html.Div(what, className="text-muted small mb-2"),
        html.Div(result, className="mb-2"),
        html.Div(implication),
    ]), className="mb-3", color=color, outline=True)


@app.callback(
    Output("fba-output", "children"),
    Input("btn-fba", "n_clicks"),
    State("store-metabolites", "data"),
    State("store-reactions", "data"),
    State("cfg-model-name", "value"),
    State("cfg-substrates", "value"),
    State("cfg-products", "value"),
    State("cfg-carbon", "value"),
    State("cfg-energy", "value"),
    State("cfg-rev", "value"),
    State("cfg-pseudo", "value"),
    prevent_initial_call=True,
)
def run_fba(_n, met_data, rxn_data, name, subs, prods, carbon, energy, rev, pseudo):
    cfg = _collect_config(name, subs, prods, carbon, energy, rev, pseudo)
    if not cfg["SUBSTRATES"] or not cfg["ENERGY_PRODUCT"]:
        return dbc.Alert("Set at least one substrate and the energy product in "
                         "the Workbench selectors.", color="warning")
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    try:
        _model, checks = run.build_and_check(df_met, df_rxn, cfg)
    except Exception as e:
        return dbc.Alert(f"Model build/FBA failed: {e}", color="danger")

    ename = _ex_name(cfg["ENERGY_PRODUCT"])
    sname = _ex_name(cfg["SUBSTRATES"][0])

    # --- Check 1: no energy-generating cycle --------------------------------
    atp = checks["atp_without_substrate"]
    atp_ok = abs(atp) < 1e-6 if atp == atp else False  # nan-safe
    if atp_ok:
        imp1 = html.Span("No energy-generating cycle — the reaction directions "
                         "are thermodynamically consistent.", className="text-success")
    else:
        fluxes = checks["atp_cycle_fluxes"] or {}
        flux_lines = "\n".join(f"{rid:<12} {val:+.3g}"
                               for rid, val in sorted(fluxes.items()))
        imp1 = html.Div([
            html.Span("This is a modelling error: the network can create energy "
                      "from nothing. It almost always means a reaction's "
                      "reversibility or direction is wrong. The reactions "
                      "carrying flux in this impossible cycle:",
                      className="text-danger"),
            html.Pre(flux_lines, className="mt-2 small"),
        ])
    card1 = _check_card(
        "pass" if atp_ok else "fail",
        "No free energy from nothing",
        f"With every substrate uptake switched off, the network should not be "
        f"able to make the energy carrier ({ename}). If it can, some reaction is "
        f"letting it run a perpetual-motion energy cycle.",
        html.Span([f"Most {ename} the model can make with no substrate: ",
                   html.B(f"{atp:.3g}"), "  (should be ≈ 0)"]),
        imp1)

    # --- Check 2: products are actually reachable ---------------------------
    yields = checks["product_yields"]
    infeasible = [p for p, y in yields.items() if y is None]
    yrows = [{"Product": _ex_name(p),
              "Max amount per 1 " + sname: "not reachable" if y is None else f"{y:.3g}"}
             for p, y in yields.items()]
    ytable = dash_table.DataTable(
        data=yrows,
        columns=[{"name": c, "id": c}
                 for c in ["Product", "Max amount per 1 " + sname]],
        style_cell={"textAlign": "left", "fontFamily": "monospace",
                    "fontSize": "12px", "padding": "4px"},
        style_header={"fontWeight": "bold"},
        style_table={"overflowX": "auto"}) if yrows else None
    if not yields:
        status2 = "info"
        imp2 = html.Span("No carbon products were selected. Add them in the "
                         "Workbench selectors if you want to check reachability.",
                         className="text-muted")
    elif infeasible:
        status2 = "warn"
        names = ", ".join(_ex_name(p) for p in infeasible)
        imp2 = html.Span(
            f"{names} cannot be produced from {sname}. Either a reaction is "
            f"missing or points the wrong way, or it isn't really a product of "
            f"this network. It will not appear in any EFM.",
            className="text-warning")
    else:
        status2 = "pass"
        imp2 = html.Span(f"Every selected product can be made from {sname}.",
                         className="text-success")
    card2 = _check_card(
        status2, "Products can be made from the substrate",
        f"Each product you flagged should be reachable from {sname}. The number "
        f"is the most product the model can make per unit of {sname} (an FBA "
        f"optimum — an upper bound, not a specific pathway).",
        ytable or html.Span("—"),
        imp2)

    # --- Check 3: network connectivity (blocked reactions) ------------------
    blocked = checks["blocked_reactions"]
    removed = checks["removed_exchanges"]

    # A blocked reaction is fatal (not merely a prunable dead end) if it is one
    # of the run-time selectors from the Workbench — substrate uptake, product
    # formation, or the energy product: if any of these can never carry flux,
    # that metabolite is silently absent from every EFM (or, for the substrate/
    # energy product, later stages error out entirely).
    fatal_ids = ({cfg["ENERGY_PRODUCT"]} | set(cfg["SUBSTRATES"])
                 | set(cfg["PRODUCTS"]) | set(cfg["CARBON_PRODUCTS"]))
    fatal_ids.discard(None)
    fatal_blocked = [b for b in blocked if b in fatal_ids]

    if blocked:
        status3 = "fail" if fatal_blocked else "warn"
        if fatal_blocked:
            imp3 = html.Span(
                f"Fatal: {', '.join(fatal_blocked)} — configured in the "
                f"Workbench as a substrate, product, or energy product, but "
                f"blocked (cannot carry any flux). Do not enumerate EFMs until "
                f"this is fixed — add the missing connecting reaction(s), or "
                f"reconsider the selector.", className="text-danger fw-bold")
        else:
            imp3 = html.Span(
                "These are inactive as the network stands — usually a missing "
                "connecting reaction or exchange. If you expected one to carry "
                "flux, check that each of its metabolites is both produced and "
                "consumed somewhere. They are pruned automatically, so the EFM "
                "step is unaffected.", className="text-warning")
        blocked_result = html.Div([
            html.Span([html.B(f"{len(blocked)}"), " reaction(s) pruned: "]),
            html.Code(", ".join(blocked)),
        ])
    else:
        status3 = "pass"
        imp3 = html.Span("Every reaction can carry flux — the network is fully "
                         "connected.", className="text-success")
        blocked_result = html.Span([html.B("0"), " blocked reactions."])
    card3 = _check_card(
        status3, "Every reaction can carry flux",
        "Blocked reactions can never carry any flux at steady state (dead "
        "ends). They are removed before enumerating EFMs.",
        blocked_result,
        imp3)

    # --- Overall banner + boundary note -------------------------------------
    if fatal_blocked:
        banner = dbc.Alert(
            html.Span([
                "Do not enumerate EFMs: ", html.B(", ".join(fatal_blocked)),
                " — a configured substrate/product/energy-product exchange — "
                "is blocked and cannot carry flux. Fix this in the Workbench "
                "first (see the connectivity check below).",
            ]), color="danger")
    elif not atp_ok:
        banner = dbc.Alert("A problem was found that would corrupt the EFM "
                           "results — fix it in the Workbench before continuing.",
                           color="danger")
    elif infeasible:
        banner = dbc.Alert("No blocking problems, but check the warnings below "
                           "before continuing.", color="warning")
    else:
        banner = dbc.Alert("All checks passed — ready to enumerate EFMs (tab ③).",
                           color="success")

    boundary = html.P(
        ["Boundary applied: exchanges removed because they were not marked as a "
         "substrate, product, or reversible exchange — ",
         html.Code(", ".join(map(str, removed)) if removed else "none"), "."],
        className="text-muted small")

    return html.Div([banner, card1, card2, card3, boundary])


# --------------------------------------------------------------------------- #
# EFM enumeration (background callback)
# --------------------------------------------------------------------------- #
@app.callback(
    Output("efm-output", "children"),
    Output("store-efm", "data"),
    Input("btn-efm", "n_clicks"),
    State("store-metabolites", "data"),
    State("store-reactions", "data"),
    State("cfg-model-name", "value"),
    State("cfg-substrates", "value"),
    State("cfg-products", "value"),
    State("cfg-carbon", "value"),
    State("cfg-energy", "value"),
    State("cfg-rev", "value"),
    State("cfg-pseudo", "value"),
    State("efm-cap", "value"),
    background=True,
    cancel=[Input("btn-efm-cancel", "n_clicks")],
    running=[
        (Output("btn-efm", "disabled"), True, False),
        (Output("efm-progress", "children"),
         "⏳ Enumerating EFMs (JVM running)…", ""),
    ],
    prevent_initial_call=True,
)
def run_efm(_n, met_data, rxn_data, name, subs, prods, carbon, energy, rev,
            pseudo, cap):
    if not efm.java_available():
        return dbc.Alert("Java (JVM) not found — efmtool cannot run.",
                         color="danger"), no_update
    cfg = _collect_config(name, subs, prods, carbon, energy, rev, pseudo)
    if not cfg["SUBSTRATES"]:
        return dbc.Alert("Set a substrate first.", color="warning"), no_update

    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    try:
        model, _checks = run.build_and_check(df_met, df_rxn, cfg)
        df_efms, df_norm, rank_failures = run.enumerate_efms(model, cfg)
    except Exception as e:
        return dbc.Alert(f"EFM enumeration failed: {e}", color="danger"), no_update

    n = df_efms.shape[1]
    cap = int(cap) if cap else 100000
    capped = n > cap

    rank_msg = ("✓ all EFMs pass the rank check (n_active − rank = 1)"
                if not rank_failures
                else f"⚠ {len(rank_failures)} EFM(s) fail the rank check: {rank_failures}")
    rank_color = "success" if not rank_failures else "warning"

    show = df_norm.reset_index().rename(columns={"index": "reaction"})
    if capped:
        keep = ["reaction"] + list(df_norm.columns[:cap])
        show = show[keep]
    # Fixed per-cell widths + a shrink-to-fit container give one consistent
    # layout regardless of EFM count: the table is exactly as wide as its
    # columns (reaction label + one 70px column per EFM) and scrolls
    # horizontally once it outgrows the viewport.  We deliberately avoid
    # `fixed_columns` (freezing the reaction column) — its two-pane rendering
    # collapses the overall table width, making a many-EFM table render
    # *narrower* than a single-EFM one.
    efm_table = dash_table.DataTable(
        data=show.round(4).to_dict("records"),
        columns=[{"name": c, "id": c} for c in show.columns],
        css=[{"selector": ".dash-spreadsheet-inner table",
              "rule": "width: auto !important; min-width: 0 !important;"}],
        style_cell={"fontFamily": "monospace", "fontSize": "12px", "padding": "3px",
                    "minWidth": "70px", "width": "70px", "maxWidth": "70px",
                    "textAlign": "right"},
        style_cell_conditional=[
            {"if": {"column_id": "reaction"},
             "textAlign": "left", "minWidth": "160px", "width": "160px",
             "maxWidth": "160px"}],
        style_table={"overflowX": "auto", "overflowY": "auto",
                     "maxHeight": "480px", "width": "fit-content",
                     "maxWidth": "100%"},
        page_size=80,
        style_data_conditional=[
            {"if": {"row_index": "odd"}, "backgroundColor": "#f7f7f7"}],
    )

    children = html.Div([
        html.H4(f"{n} EFM{'s' if n != 1 else ''}"),
        dbc.Alert(rank_msg, color=rank_color),
        (dbc.Alert(f"Showing first {cap} of {n} (cap). Raise the cap to see all.",
                   color="warning") if capped else None),
        html.H6("Normalised EFM table (per unit substrate uptake)"),
        efm_table,
    ])
    store = {"normalized": df_norm.to_json(orient="split"),
             "n_efms": n, "config": _config_for_store(cfg)}
    return children, store


def _config_for_store(cfg):
    c = dict(cfg)
    c["PSEUDO_METS"] = sorted(cfg.get("PSEUDO_METS") or [])
    return c


# --------------------------------------------------------------------------- #
# Thermodynamics
# --------------------------------------------------------------------------- #
@app.callback(
    Output("thermo-output", "children"),
    Output("store-thermo", "data"),
    Input("btn-thermo", "n_clicks"),
    State("store-efm", "data"),
    State("store-metabolites", "data"),
    State("store-reactions", "data"),
    State("thermo-which", "value"),
    background=True,
    running=[(Output("btn-thermo", "disabled"), True, False)],
    prevent_initial_call=True,
)
def run_thermo(_n, efm_store, met_data, rxn_data, which):
    if not efm_store:
        return dbc.Alert("Enumerate EFMs first (tab ③).", color="warning"), no_update
    df_norm = pd.read_json(_io.StringIO(efm_store["normalized"]), orient="split")
    cfg = efm_store["config"]
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    try:
        model, _ = run.build_and_check(df_met, df_rxn,
                                       {**cfg, "PSEUDO_METS": set(cfg["PSEUDO_METS"])})
        per_efm, net_df = thermo.compute_thermodynamics(
            df_norm, df_met, model, cfg["PSEUDO_METS"], which=which)
    except Exception as e:
        return dbc.Alert(f"Thermodynamics failed: {e}", color="danger"), no_update

    disp = _format_net_table(net_df, which)
    table = dash_table.DataTable(
        data=disp.to_dict("records"),
        columns=[{"name": c, "id": c} for c in disp.columns],
        style_cell={"fontFamily": "monospace", "fontSize": "12px", "textAlign": "left",
                    "padding": "4px", "whiteSpace": "normal", "height": "auto"},
        style_header={"fontWeight": "bold"},
        style_table={"overflowX": "auto"},
    )
    notes = [html.Li(f"{r['net_conversion']}: {r['thermo_error']}")
             for _, r in net_df.iterrows() if pd.notna(r.get("thermo_error"))]
    details = _per_efm_details(per_efm, which)
    children = html.Div([
        html.H5(f"{len(net_df)} distinct net conversion(s) "
                f"from {efm_store['n_efms']} EFMs"),
        html.P("Ω = −ΔG_CAT / Σ(vᵢ/v_NC) over metabolic reactions (= MDF, no "
               "concentration bounds).", className="text-muted small"),
        table,
        (html.Div([
            html.H6("Per-EFM breakdown", className="mt-4"),
            html.P("Net conversions produced by more than one EFM: each EFM "
                   "shares the same net-conversion ΔG but has its own flux sum, "
                   "so Ω can differ per EFM.",
                   className="text-muted small"),
            *details,
        ]) if details else None),
        (html.Div([html.B("Thermo notes:"), html.Ul(notes)]) if notes else None),
    ])
    store = {"net": net_df.to_json(orient="split"),
             "per_efm": per_efm.to_json(orient="split")}
    return children, store


def _per_efm_details(per_efm, which):
    """Collapsible per-EFM tables for net conversions produced by >1 EFM.

    ΔG is a property of the net conversion (shared by all its EFMs), so it is
    shown once in the section header; the per-EFM table carries the quantities
    that actually differ between EFMs — the flux sum, n_R and Ω.
    """
    if per_efm.empty:
        return []

    def fmt(v):
        return f"{v:.2f}" if isinstance(v, (int, float)) and pd.notna(v) else "—"

    dg_label = "ΔGm′" if which == "dGm" else "ΔG°′"
    blocks = []
    for nc, grp in per_efm.groupby("net_conversion", sort=False):
        if len(grp) < 2:
            continue
        first = grp.iloc[0]
        dg_val = first["dGm value"] if which == "dGm" else first["dG0prime value"]
        dg_err = first["dGm error"] if which == "dGm" else first["dG0prime error"]
        sub = pd.DataFrame({
            "EFM": grp["EFM"],
            "n_R": [fmt(v) for v in grp["n_metabolic_reactions"]],
            "Σ(vᵢ/v_NC)": [fmt(v) for v in grp["flux_sum"]],
            f"Ω (from {dg_label}, kJ/mol)": [fmt(v) for v in grp["Omega"]],
        })
        sub_table = dash_table.DataTable(
            data=sub.to_dict("records"),
            columns=[{"name": c, "id": c} for c in sub.columns],
            style_cell={"fontFamily": "monospace", "fontSize": "12px",
                        "textAlign": "left", "padding": "4px"},
            style_header={"fontWeight": "bold"},
            style_table={"overflowX": "auto"},
        )
        blocks.append(html.Details([
            html.Summary(
                f"{nc}  —  {len(grp)} EFMs, "
                f"{dg_label} = {fmt(dg_val)} ± {fmt(dg_err)} kJ/mol"),
            html.Div(sub_table, className="my-2"),
        ], open=True))
    return blocks


def _format_net_table(net_df, which):
    df = net_df.copy()

    def fmt(v):
        return f"{v:.2f}" if isinstance(v, (int, float)) and pd.notna(v) else "—"

    def rng(lo, hi):
        if pd.isna(lo) and pd.isna(hi):
            return "—"
        return fmt(lo) if abs((lo or 0) - (hi or 0)) < 1e-9 else f"{fmt(lo)} … {fmt(hi)}"

    out = pd.DataFrame({
        "Net conversion": df["net_conversion"],
        "EFMs": df["EFM_multiplicity"],
        "ΔG°′ (kJ/mol)": [f"{fmt(v)} ± {fmt(e)}" for v, e in
                          zip(df["dG0prime value"], df["dG0prime error"])],
        "ΔGm′ (kJ/mol)": [f"{fmt(v)} ± {fmt(e)}" for v, e in
                          zip(df["dGm value"], df["dGm error"])],
        f"Ω (from {'ΔGm′' if which == 'dGm' else 'ΔG°′'}, kJ/mol)":
            [rng(lo, hi) for lo, hi in zip(df["Omega_min"], df["Omega_max"])],
    })
    return out


# --------------------------------------------------------------------------- #
# Downloads: model JSON/SBML, results XLSX
# --------------------------------------------------------------------------- #
def _build_model_for_download(met_data, rxn_data, name):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    model, _ = model_mod.build_model(df_met, df_rxn, name or "model")
    model_mod.sanitize_model(model)
    return model


@app.callback(
    Output("download-json", "data"),
    Input("btn-dl-json", "n_clicks"),
    State("store-metabolites", "data"),
    State("store-reactions", "data"),
    State("cfg-model-name", "value"),
    prevent_initial_call=True,
)
def dl_json(_n, met_data, rxn_data, name):
    import cobra
    model = _build_model_for_download(met_data, rxn_data, name)
    with tempfile.NamedTemporaryFile("w+", suffix=".json", delete=False) as fh:
        cobra.io.save_json_model(model, fh.name)
        path = fh.name
    with open(path) as fh:
        text = fh.read()
    os.unlink(path)
    return dict(content=text, filename=f"{name or 'model'}.json")


@app.callback(
    Output("download-sbml", "data"),
    Input("btn-dl-sbml", "n_clicks"),
    State("store-metabolites", "data"),
    State("store-reactions", "data"),
    State("cfg-model-name", "value"),
    prevent_initial_call=True,
)
def dl_sbml(_n, met_data, rxn_data, name):
    import cobra
    model = _build_model_for_download(met_data, rxn_data, name)
    with tempfile.NamedTemporaryFile("w+", suffix=".xml", delete=False) as fh:
        cobra.io.write_sbml_model(model, fh.name)
        path = fh.name
    with open(path) as fh:
        text = fh.read()
    os.unlink(path)
    return dict(content=text, filename=f"{name or 'model'}.xml")


@app.callback(
    Output("download-results", "data"),
    Output("dl-status", "children"),
    Input("btn-dl-results", "n_clicks"),
    State("store-efm", "data"),
    State("store-thermo", "data"),
    State("cfg-model-name", "value"),
    prevent_initial_call=True,
)
def dl_results(_n, efm_store, thermo_store, name):
    if not efm_store:
        return no_update, "Run EFM enumeration first."
    buf = _io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xl:
        df_norm = pd.read_json(_io.StringIO(efm_store["normalized"]), orient="split")
        df_norm.reset_index().rename(columns={"index": "reaction"}).to_excel(
            xl, sheet_name="EFM_table", index=False)
        if thermo_store:
            pd.read_json(_io.StringIO(thermo_store["net"]), orient="split").to_excel(
                xl, sheet_name="Net_conversions", index=False)
            pd.read_json(_io.StringIO(thermo_store["per_efm"]), orient="split").to_excel(
                xl, sheet_name="Per_EFM_thermo", index=False)
    return (dcc.send_bytes(buf.getvalue(), f"{name or 'results'}_results.xlsx"),
            "Results downloaded ✓")


if __name__ == "__main__":
    app.run(debug=False, host="127.0.0.1", port=8050)
