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
from dash import (Dash, DiskcacheManager, Input, Output, State, ctx, dash_table,
                  dcc, html, no_update)

import diskcache

from .pipeline import efm, io, model as model_mod, run, thermo

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


def workbench_tab():
    # Seed the tables with the default example directly in the layout so they
    # are never empty on first paint.  (A DataTable mounted with no data can
    # render blank until a reflow — e.g. switching tabs — forces a redraw.)
    df_met0, df_rxn0, _cfg0 = load_example("Example1_EMP_lactate.xlsx")
    met0, rxn0 = df_to_records(df_met0), df_to_records(df_rxn0)
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
        dbc.Row([
            dbc.Col([
                html.H5("Reactions"),
                dbc.Button("+ Add reaction row", id="btn-add-rxn", size="sm",
                           color="light", className="mb-1"),
                _table("table-reactions", io.REACTION_COLS, data=rxn0),
            ], md=7),
            dbc.Col([
                html.H5("Metabolites"),
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
                dbc.Select(id="fix-reaction", options=[]),
            ], className="mb-2"),
            dbc.ButtonGroup([
                dbc.Button("+ H⁺ left", id="fix-h-left", size="sm", outline=True, color="primary"),
                dbc.Button("+ H⁺ right", id="fix-h-right", size="sm", outline=True, color="primary"),
                dbc.Button("+ H₂O left", id="fix-w-left", size="sm", outline=True, color="info"),
                dbc.Button("+ H₂O right", id="fix-w-right", size="sm", outline=True, color="info"),
            ]),
        ]), className="mb-3"),

        html.Hr(),
        html.H5("4 · Run-time selectors"),
        html.Div("Scientific choices (reversibility/direction, which metabolites "
                 "are external) stay with you — these are populated from the "
                 "exchanges present.", className="text-muted small mb-2"),
        dbc.Row([
            dbc.Col([dbc.Label("Model name"),
                     dbc.Input(id="cfg-model-name", value="model")], md=4),
            dbc.Col([dbc.Label("Substrate(s)"),
                     dcc.Dropdown(id="cfg-substrates", multi=True)], md=4),
            dbc.Col([dbc.Label("Energy product (ΔG cycle check)"),
                     dcc.Dropdown(id="cfg-energy")], md=4),
        ], className="mb-2"),
        dbc.Row([
            dbc.Col([dbc.Label("Products (secretable)"),
                     dcc.Dropdown(id="cfg-products", multi=True)], md=4),
            dbc.Col([dbc.Label("Carbon products (yield check)"),
                     dcc.Dropdown(id="cfg-carbon", multi=True)], md=4),
            dbc.Col([dbc.Label("Freely-reversible exchanges"),
                     dcc.Dropdown(id="cfg-rev", multi=True)], md=4),
        ], className="mb-2"),
        dbc.Row([
            dbc.Col([dbc.Label("Pseudo (currency) metabolites — excluded from Ω flux sum"),
                     dcc.Dropdown(id="cfg-pseudo", multi=True)], md=12),
        ], className="mb-3"),
    ], fluid=True)


def fba_tab():
    return dbc.Container([
        dbc.Button("Build model & run FBA checks", id="btn-fba", color="primary",
                   className="my-3"),
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
        dbc.Tab(fba_tab(), label="② Model & FBA", tab_id="tab-fba"),
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


# Live balance panel + populate selectors from the current tables.
@app.callback(
    Output("balance-panel", "children"),
    Output("fix-reaction", "options"),
    Output("cfg-substrates", "options"),
    Output("cfg-products", "options"),
    Output("cfg-carbon", "options"),
    Output("cfg-energy", "options"),
    Output("cfg-rev", "options"),
    Output("cfg-pseudo", "options"),
    Input("store-metabolites", "data"),
    Input("store-reactions", "data"),
)
def update_balance(met_data, rxn_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if df_rxn.empty:
        return (html.Div("No reactions yet."), [], [], [], [], [], [], [])

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

    return (html.Div([banner, table]), flagged,
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
    if trig == "upload-excel" and upload:
        try:
            _h, b64 = upload.split(",", 1)
            cfg = io.load_config(base64.b64decode(b64)) or {}
        except Exception:
            cfg = {}
    else:
        _m, _r, cfg = load_example(example_file or "Example1_EMP_lactate.xlsx")

    # Suggest currency pseudo-mets present in the table.
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
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

    atp = checks["atp_without_substrate"]
    atp_ok = abs(atp) < 1e-6 if atp == atp else False  # nan-safe
    atp_card = dbc.Alert(
        [html.B(f"Max {cfg['ENERGY_PRODUCT']} with substrate blocked: {atp:.4g}"),
         html.Br(),
         "≈ 0 ✓ no energy-generating cycle" if atp_ok
         else "⚠ nonzero — energy-generating cycle present!"],
        color="success" if atp_ok else "danger")
    cycle = (html.Pre(str(checks["atp_cycle_fluxes"]))
             if checks["atp_cycle_fluxes"] else None)

    yields = checks["product_yields"]
    yrows = [{"Product": p[2:] if p.startswith("EX") else p,
              "Max yield / substrate": "infeasible" if y is None else f"{y:.4g}"}
             for p, y in yields.items()]
    ytable = dash_table.DataTable(
        data=yrows, columns=[{"name": c, "id": c} for c in ["Product", "Max yield / substrate"]],
        style_cell={"textAlign": "left", "fontFamily": "monospace"})

    removed = checks["removed_exchanges"]
    blocked = checks["blocked_reactions"]
    return html.Div([
        atp_card, cycle,
        html.H6("Max individual product yields"), ytable,
        html.Hr(),
        html.P(f"Exchanges removed (not in substrate/product/rev sets): "
               f"{removed or 'none'}"),
        html.P(f"Blocked reactions pruned ({len(blocked)}): {blocked or 'none'}"),
        dbc.Alert("Model built and pruned. Proceed to the EFMs tab.", color="info"),
    ])


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
               "concentration bounds). The equal-flux bound −ΔG/n_R coincides "
               "with Ω only when all fluxes are equal.", className="text-muted small"),
        table,
        (html.Div([
            html.H6("Per-EFM breakdown", className="mt-4"),
            html.P("Net conversions produced by more than one EFM: each EFM "
                   "shares the same net-conversion ΔG but has its own flux sum, "
                   "so Ω and the equal-flux bound differ per EFM.",
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
    that actually differ between EFMs — the flux sum, n_R, Ω and the bound.
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
            "Equal-flux bound −ΔG/n_R": [fmt(v) for v in grp["equal_flux_bound"]],
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
        "Equal-flux bound −ΔG/n_R": [rng(lo, hi) for lo, hi in
                                     zip(df["equal_flux_bound_min"],
                                         df["equal_flux_bound_max"])],
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
