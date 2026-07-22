"""KEGG workbook builder.

A small, single-purpose companion to the main analysis app: assemble the two
canonical dataframes (Metabolites, Reactions) **from KEGG** — reaction ids or a
whole module — check balance live, patch the usual KEGG omissions (H+/H2O), and
download the ``.xlsx`` the analysis app (or the paper's supplements) consumes.

It depends only on the shared :mod:`pipeline` package, needs **no Java** (there
is no efmtool/EFM step here), and serves ``http://127.0.0.1:8051`` so it can run
alongside the analysis app.

Run from the repo root:  ``python -m kegg_builder_app.app``
"""

from __future__ import annotations

import base64

import dash_bootstrap_components as dbc
import pandas as pd
from dash import (ALL, MATCH, Dash, Input, Output, State, ctx, dash_table, dcc,
                  html, no_update)

from pipeline import bigg, exchanges, idmap, io, kegg, run
from pipeline import model as model_mod

app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY],
           title="KEGG workbook builder")
server = app.server


# --------------------------------------------------------------------------- #
# Serialization helpers (dataframes <-> DataTable records)
# --------------------------------------------------------------------------- #
def df_to_records(df):
    return df.to_dict("records")


def records_to_df(records, cols):
    df = pd.DataFrame(records or [])
    for c in cols:
        if c not in df.columns:
            df[c] = pd.Series(dtype=object)
    return df[cols]


def _table(table_id, columns, data=None):
    return dash_table.DataTable(
        id=table_id,
        columns=[{"name": c, "id": c} for c in columns],
        data=data or [],
        editable=True,
        row_deletable=True,
        page_size=100,
        style_table={"overflowX": "auto", "maxHeight": "360px", "overflowY": "auto"},
        style_cell={"fontFamily": "monospace", "fontSize": "13px",
                    "textAlign": "left", "padding": "4px"},
        style_header={"fontWeight": "bold"},
    )


# --------------------------------------------------------------------------- #
# Manual reaction entry: dynamic substrate/product rows (dcc components, not
# DataTable cells — the DataTable's cell editor swallows Backspace/Delete/
# Home/arrow keys for its own cell-to-cell navigation, so composing a new
# reaction character-by-character in a table cell is effectively unfixable
# without replacing the editor). Each row gets a stable, unique pattern-
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


# --------------------------------------------------------------------------- #
# Layout
# --------------------------------------------------------------------------- #
app.layout = dbc.Container([
    html.H3("KEGG workbook builder", className="mt-3"),
    html.P("Assemble the Metabolites/Reactions workbook from KEGG, check its "
           "balance, and download the .xlsx for the analysis app. No Java "
           "needed — this app does not run efmtool.", className="text-muted"),

    dbc.Card(dbc.CardBody([
        html.H5("1 · Fetch from KEGG"),
        html.P("Paste KEGG reaction ids (e.g. R00200, R01786) and/or module ids "
               "(e.g. M00001 for EMP glycolysis), separated by spaces, commas or "
               "new lines. Fetching appends to the tables below.",
               className="text-muted small"),
        dbc.Textarea(id="kegg-input", placeholder="R00200, R01786  M00001",
                     style={"height": "70px", "fontFamily": "monospace"},
                     className="mb-2"),
        dbc.Button("⬇ Fetch from KEGG", id="btn-fetch", color="primary"),
        dcc.Loading(html.Div(id="fetch-messages",
                             style={"maxHeight": "160px", "overflowY": "auto",
                                    "fontSize": "12px", "fontFamily": "monospace",
                                    "marginTop": "0.5rem"})),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("1a · Human-readable ids (optional)"),
        html.P("Translate KEGG compound ids to BiGG ids (e.g. C00031 → glc__D) "
               "and KEGG reaction ids to BiGG reaction ids (e.g. R00200 → PYK). "
               "The KEGG ID / KEGG Reaction ID columns are kept, so balance and "
               "thermodynamics still resolve everything correctly. Anything "
               "with no clean BiGG id falls back to a slugified KEGG name and "
               "is flagged for you to hand-edit. First use downloads two small "
               "BiGG files (needs network).",
               className="text-muted small"),
        dbc.Checklist(
            id="chk-use-mnx",
            options=[{
                "label": ("Also try MetaNetX cross-references for ids BiGG "
                          "doesn't have — finds many more matches (BiGG's own "
                          "files link a KEGG id for only ~10% of its "
                          "reactions), but is a large one-time download "
                          "(~700MB compounds + ~80MB reactions) and can take "
                          "several minutes on first use. Already have "
                          "chem_xref.tsv / reac_xref.tsv from metanetx.org? "
                          "Drop them in .cache/ (next to the pipeline "
                          "package) to skip the download entirely."),
                "value": "mnx",
            }],
            value=[], switch=True, className="small text-muted mb-2",
        ),
        dbc.Button("Translate KEGG → BiGG ids", id="btn-bigg",
                   color="secondary", outline=True),
        dcc.Loading(html.Div(id="bigg-messages",
                             style={"maxHeight": "220px", "overflowY": "auto",
                                    "fontSize": "12px", "marginTop": "0.5rem"})),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("1b · Custom metabolite ids (optional)"),
        html.P("Set your own working ids. Click Load to fill the table from the "
               "current metabolites, edit the New ID column, then Apply — the "
               "reaction stoichiometries and EX/transport ids are rewritten to "
               "match, and the KEGG ID column is kept. Do this before adding "
               "exchanges so those pick up your ids. Press Enter after editing "
               "a cell to save it — clicking away without Enter does not save.",
               className="text-muted small"),
        dbc.Button("Load current metabolites", id="btn-load-rename",
                   color="secondary", outline=True, size="sm", className="mb-2"),
        dash_table.DataTable(
            id="tbl-rename",
            columns=[{"name": "KEGG ID", "id": "KEGG ID", "editable": False},
                     {"name": "Current ID", "id": "Current ID", "editable": False},
                     {"name": "New ID", "id": "New ID", "editable": True}],
            data=[], editable=True, page_size=100,
            style_table={"overflowX": "auto", "maxHeight": "300px",
                         "overflowY": "auto"},
            style_cell={"fontFamily": "monospace", "fontSize": "12px",
                        "textAlign": "left", "padding": "4px"},
            style_header={"fontWeight": "bold"},
            style_data_conditional=[
                {"if": {"column_id": "New ID"}, "backgroundColor": "#f0f7ff"}],
        ),
        dbc.Button("Apply custom ids", id="btn-apply-rename",
                   color="primary", className="mt-2 w-100"),
        html.Div(id="rename-messages",
                 style={"fontSize": "12px", "marginTop": "0.5rem"}),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("1c · Custom reaction ids (optional)"),
        html.P("Same idea for reactions — handy after KEGG gives you bare "
               "Rxxxxx ids. Click Load to fill the table from the current "
               "reactions, edit the New ID column, then Apply. Reaction ids "
               "aren't referenced anywhere else in the tables (only "
               "metabolite ids are, inside stoichiometries), so this never "
               "touches your equations. Press Enter after editing a cell to "
               "save it — clicking away without Enter does not save.",
               className="text-muted small"),
        dbc.Button("Load current reactions", id="btn-load-rxn-rename",
                   color="secondary", outline=True, size="sm", className="mb-2"),
        dash_table.DataTable(
            id="tbl-rxn-rename",
            columns=[{"name": "KEGG Reaction ID", "id": "KEGG Reaction ID",
                     "editable": False},
                     {"name": "Current ID", "id": "Current ID", "editable": False},
                     {"name": "Stoichiometry", "id": "Stoichiometry",
                     "editable": False},
                     {"name": "New ID", "id": "New ID", "editable": True}],
            data=[], editable=True, page_size=100,
            style_table={"overflowX": "auto", "maxHeight": "300px",
                         "overflowY": "auto"},
            style_cell={"fontFamily": "monospace", "fontSize": "12px",
                        "textAlign": "left", "padding": "4px"},
            style_header={"fontWeight": "bold"},
            style_data_conditional=[
                {"if": {"column_id": "New ID"}, "backgroundColor": "#f0f7ff"}],
        ),
        dbc.Button("Apply custom reaction ids", id="btn-apply-rxn-rename",
                   color="primary", className="mt-2 w-100"),
        html.Div(id="rxn-rename-messages",
                 style={"fontSize": "12px", "marginTop": "0.5rem"}),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("1d · Add a reaction manually (non-KEGG)"),
        html.P("For transport, PTS import, ion translocation or designed "
               "reactions KEGG has no entry for. Pick substrates and products "
               "from existing metabolites, or type a working id that doesn't "
               "exist yet and choose \"＋ New metabolite: …\" to create it as "
               "a placeholder (fill in its KEGG id or chemical formula "
               "afterwards so it can balance). Add rows for extra "
               "reactants/products as needed.",
               className="text-muted small"),
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
                html.Div(id="man-subs-rows", children=[_metab_row("subs", 0, [])]),
                dbc.Button("+ Add metabolite", id="btn-man-add-subs", size="sm",
                           color="light"),
            ], md=6),
            dbc.Col([
                html.Div("Products", className="fw-bold small mb-1"),
                html.Div(id="man-prods-rows", children=[_metab_row("prods", 0, [])]),
                dbc.Button("+ Add metabolite", id="btn-man-add-prods", size="sm",
                           color="light"),
            ], md=6),
        ], className="mb-2"),
        html.Div(id="manual-messages",
                 style={"fontSize": "12px", "marginTop": "0.5rem"}),
        html.Div(id="man-new-met-panel"),
    ]), className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.H5("Reactions"),
            dbc.Button("+ Add reaction row", id="btn-add-rxn", size="sm",
                       color="light", className="mb-1"),
            html.Div("Press Enter after editing a cell to save it — clicking "
                     "away without Enter does not save.",
                     className="text-muted small mb-1"),
            _table("tbl-rxn", io.REACTION_COLS),
        ], md=7),
        dbc.Col([
            html.H5("Metabolites"),
            html.Div(id="met-missing-warning"),
            dbc.Button("+ Add metabolite row", id="btn-add-met", size="sm",
                       color="light", className="mb-1"),
            html.Div("Press Enter after editing a cell to save it — clicking "
                     "away without Enter does not save.",
                     className="text-muted small mb-1"),
            _table("tbl-met", io.METABOLITE_COLS),
        ], md=5),
    ], className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("Remove reactions"),
        html.P("Trim reactions you don't want — handy after loading a KEGG "
               "module that pulls in more than you need. Optionally drop any "
               "metabolites left unused afterwards.", className="text-muted small"),
        dcc.Dropdown(id="remove-rxns", multi=True, className="mb-2",
                     placeholder="Reactions to remove…"),
        dbc.ButtonGroup([
            dbc.Button("Remove selected reactions", id="btn-remove-rxns",
                       color="danger", outline=True, size="sm"),
            dbc.Button("Remove unused metabolites", id="btn-prune-mets",
                       color="secondary", outline=True, size="sm"),
        ]),
        html.Div(id="remove-messages",
                 style={"fontSize": "12px", "marginTop": "0.5rem"}),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("Reaction reversibility"),
        html.P("New reactions default to irreversible — reversibility is a "
               "deliberate modelling choice (CLAUDE.md §5), never a silent "
               "default. Pick the reactions that are actually reversible and "
               "mark them here; check the Reversibility column in the "
               "Reactions table above to see the current state.",
               className="text-muted small"),
        dcc.Dropdown(id="reversible-rxns", multi=True, className="mb-2",
                     placeholder="Reactions to mark…"),
        dbc.ButtonGroup([
            dbc.Button("Mark selected reversible", id="btn-mark-reversible",
                       color="primary", outline=True, size="sm"),
            dbc.Button("Mark selected irreversible", id="btn-mark-irreversible",
                       color="secondary", outline=True, size="sm"),
        ]),
        html.Div(id="reversibility-messages",
                 style={"fontSize": "12px", "marginTop": "0.5rem"}),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("1e · Exchange & transport reactions"),
        html.P("Pick the metabolites that cross the system boundary. For each, "
               "the app adds an extracellular counterpart, a transport reaction "
               "(X ⇌ Xex) and an exchange reaction (EXX). Direction stays a "
               "modelling choice — everything is created reversible; set "
               "substrates/products later in the analysis app.",
               className="text-muted small"),
        dcc.Dropdown(id="exchange-mets", multi=True, className="mb-2",
                     placeholder="Select boundary metabolites…"),
        dbc.Button("Add exchange + transport reactions", id="btn-exchange",
                   color="secondary", outline=True),
        html.Div(id="exchange-messages",
                 style={"fontSize": "12px", "fontFamily": "monospace",
                        "marginTop": "0.5rem"}),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("1f · Network connectivity check"),
        html.P("Build a model straight from the tables above (every reaction "
               "reversible except those marked irreversible — no substrate/"
               "product/energy-product choices yet, those come later in the "
               "analysis app) and find reactions that can never carry any "
               "flux at steady state. Run this after adding exchange/"
               "transport reactions: a blocked exchange usually means a "
               "missing transport or connecting reaction for that boundary "
               "metabolite.", className="text-muted small"),
        dbc.Button("Check for blocked reactions", id="btn-connectivity",
                   color="secondary", outline=True),
        dcc.Loading(html.Div(id="connectivity-output", className="mt-2")),
    ]), className="mb-3"),

    html.Hr(),
    html.H5("2 · Live balance check"),
    html.P("KEGG reactions often omit H⁺/H₂O or aren't balanced as written. "
           "Exchange reactions (EX…) are expected to be unbalanced and are not "
           "flagged. Use the one-click fixes to patch a flagged reaction.",
           className="text-muted small"),
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
    html.H5("3 · Save / load"),
    dbc.Row([
        dbc.Col(dcc.Upload(
            id="upload-xlsx",
            children=html.Div(["⬆ Drag/drop or ", html.A("select an .xlsx")]),
            style={"border": "1px dashed #aaa", "borderRadius": "6px",
                   "padding": "10px", "textAlign": "center"}), md=6),
        dbc.Col([
            dbc.Input(id="wb-name", value="kegg_workbook", className="mb-2"),
            dbc.Button("⬇ Download workbook (.xlsx)", id="btn-download",
                       color="secondary", outline=True, className="w-100"),
            dcc.Download(id="download-xlsx"),
        ], md=6),
    ], className="mb-4"),
], fluid=True)


# --------------------------------------------------------------------------- #
# Callbacks
# --------------------------------------------------------------------------- #
@app.callback(
    Output("tbl-met", "data"),
    Output("tbl-rxn", "data"),
    Output("fetch-messages", "children"),
    Input("btn-fetch", "n_clicks"),
    Input("upload-xlsx", "contents"),
    Input("btn-add-rxn", "n_clicks"),
    Input("btn-add-met", "n_clicks"),
    State("kegg-input", "value"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def update_tables(_nf, upload, _nar, _nam, kegg_text, met_data, rxn_data):
    trig = ctx.triggered_id
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)

    if trig == "btn-fetch":
        ids = [t for t in _split_ids(kegg_text)]
        if not ids:
            return no_update, no_update, "Enter at least one KEGG id."
        df_met, df_rxn, messages = kegg.fetch_reactions(
            ids, existing_metabolites=df_met, existing_reactions=df_rxn)
        msg = html.Ul([html.Li(m) for m in messages]) if messages else "Nothing fetched."
        return df_to_records(df_met), df_to_records(df_rxn), msg

    if trig == "upload-xlsx" and upload:
        _header, b64 = upload.split(",", 1)
        try:
            df_met, df_rxn = io.load_excel(base64.b64decode(b64))
        except Exception as e:
            return no_update, no_update, f"⚠ Could not read workbook: {e}"
        return (df_to_records(df_met), df_to_records(df_rxn),
                "Workbook uploaded ✓")

    if trig == "btn-add-rxn":
        df_rxn = pd.concat([df_rxn, pd.DataFrame([{c: "" for c in io.REACTION_COLS}])],
                           ignore_index=True)
    elif trig == "btn-add-met":
        df_met = pd.concat([df_met, pd.DataFrame([{c: "" for c in io.METABOLITE_COLS}])],
                           ignore_index=True)
    return df_to_records(df_met), df_to_records(df_rxn), no_update


def _split_ids(text):
    if not text:
        return []
    out = []
    for tok in str(text).replace(",", " ").split():
        tok = tok.strip()
        if tok:
            out.append(tok)
    return out


_BIGG_STATUS_LABEL = {
    "bigg-ambiguous": "ambiguous — picked shortest BiGG id",
    "name-fallback": "no BiGG id — used slugified KEGG name",
    "kegg-fallback": "no KEGG id — left unchanged",
    "collision": "id clashed — suffixed to keep it unique",
    "merged": "duplicate of an existing metabolite — merged",
}


@app.callback(
    Output("tbl-met", "data", allow_duplicate=True),
    Output("tbl-rxn", "data", allow_duplicate=True),
    Output("bigg-messages", "children"),
    Input("btn-bigg", "n_clicks"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    State("chk-use-mnx", "value"),
    prevent_initial_call=True,
)
def translate_bigg(_n, met_data, rxn_data, use_mnx_value):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if df_met.empty or df_met["ID"].dropna().empty:
        return no_update, no_update, "Fetch some reactions first."
    use_mnx = bool(use_mnx_value)
    try:
        m2, r2, met_report = bigg.translate_to_bigg(df_met, df_rxn, use_mnx=use_mnx)
        r2, rxn_report = bigg.translate_reactions_to_bigg(r2, use_mnx=use_mnx)
    except bigg.BiggError as e:
        return no_update, no_update, dbc.Alert(
            f"⚠ {e} — needs network on first use.", color="warning")

    def _summary_bits(report, label):
        n = len(report)
        clean = sum(1 for x in report if x["status"] == "bigg")
        flags = [x for x in report if x["status"] != "bigg"]
        text = f"{label}: {clean}/{n} clean BiGG ids"
        if flags:
            text += f" ({len(flags)} need a look)"
        return text, flags

    met_text, met_flags = _summary_bits(met_report, "Metabolites")
    rxn_text, rxn_flags = _summary_bits(rxn_report, "Reactions")
    summary = dbc.Alert(f"{met_text}. {rxn_text}.",
                        color="success" if not (met_flags or rxn_flags) else "info")

    def _detail(flags, label):
        if not flags:
            return None
        return html.Div([
            html.B(label),
            html.Ul([
                html.Li([html.Code(f"{x['kegg'] or x['old']} → {x['new']}"),
                         f"  ({_BIGG_STATUS_LABEL.get(x['status'], x['status'])})"])
                for x in flags]),
        ])

    detail_children = [d for d in
                       [_detail(met_flags, "Metabolites:"), _detail(rxn_flags, "Reactions:")]
                       if d]
    return (df_to_records(m2), df_to_records(r2),
            html.Div([summary] + detail_children))


@app.callback(
    Output("man-subs-rows", "children"),
    Input("btn-man-add-subs", "n_clicks"),
    State("man-subs-rows", "children"),
    State("tbl-met", "data"),
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
    State("tbl-met", "data"),
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
    State("tbl-met", "data"),
)
def update_man_met_options(search_value, current_value, met_data):
    # Re-runs both when the user types (search_value) and right after they
    # pick something (value) — a freshly-picked "new metabolite" doesn't
    # exist in tbl-met yet (it's only created on submit), so without this
    # the very next options refresh (e.g. the search box closing) would drop
    # it from the list and the dropdown would lose its selection silently.
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    known_ids = df_met["ID"].dropna().astype(str).tolist()
    opts = _met_options(search_value, known_ids)
    if (current_value and current_value not in known_ids
            and not any(o["value"] == current_value for o in opts)):
        opts.append({"label": f"＋ New metabolite: {current_value}",
                     "value": current_value})
    return opts


@app.callback(
    Output("tbl-met", "data", allow_duplicate=True),
    Output("tbl-rxn", "data", allow_duplicate=True),
    Output("man-rxn-id", "value"),
    Output("man-rxn-name", "value"),
    Output("man-subs-rows", "children", allow_duplicate=True),
    Output("man-prods-rows", "children", allow_duplicate=True),
    Output("manual-messages", "children"),
    Output("man-new-met-panel", "children"),
    Input("btn-add-manual-rxn", "n_clicks"),
    State("man-rxn-id", "value"),
    State("man-rxn-name", "value"),
    State("man-rxn-rev", "value"),
    State({"type": "man-coeff", "side": "subs", "index": ALL}, "value"),
    State({"type": "man-met", "side": "subs", "index": ALL}, "value"),
    State({"type": "man-coeff", "side": "prods", "index": ALL}, "value"),
    State({"type": "man-met", "side": "prods", "index": ALL}, "value"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def add_manual_reaction(_n, rid, name, rev, subs_c, subs_m, prods_c, prods_m,
                        met_data, rxn_data):
    rid = (rid or "").strip()
    if not rid:
        return (no_update,) * 6 + (dbc.Alert("Enter a reaction id.", color="warning"),
                                    no_update)

    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if (df_rxn["ID"].astype(str).str.strip() == rid).any():
        return (no_update,) * 6 + (dbc.Alert(f"⚠ Reaction '{rid}' already exists.",
                                              color="warning"),
                                    no_update)

    subs_terms = _side_terms(subs_c, subs_m)
    prods_terms = _side_terms(prods_c, prods_m)
    if not subs_terms or not prods_terms:
        return (no_update,) * 6 + (dbc.Alert("Pick at least one substrate and one "
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
            "", "",
            [_metab_row("subs", 0, fresh_known)],
            [_metab_row("prods", 0, fresh_known)],
            dbc.Alert(" — ".join(parts), color="success"),
            _new_met_panel([m["ID"] for m in new_mets]))


@app.callback(
    Output("tbl-met", "data", allow_duplicate=True),
    Output("man-new-met-panel", "children", allow_duplicate=True),
    Input("btn-save-new-mets", "n_clicks"),
    State({"type": "new-met-kegg", "mid": ALL}, "value"),
    State({"type": "new-met-kegg", "mid": ALL}, "id"),
    State({"type": "new-met-formula", "mid": ALL}, "value"),
    State("tbl-met", "data"),
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
    return df_to_records(df_met), None


def _is_blank(v):
    return v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() == ""


@app.callback(
    Output("met-missing-warning", "children"),
    Input("tbl-met", "data"),
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


@app.callback(
    Output("tbl-rename", "data"),
    Input("btn-load-rename", "n_clicks"),
    State("tbl-met", "data"),
    prevent_initial_call=True,
)
def populate_rename(_n, met_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    rows = []
    for _, r in df_met.iterrows():
        mid = str(r.get("ID", "")).strip()
        if not mid:
            continue
        rows.append({"KEGG ID": r.get("KEGG ID"), "Current ID": mid,
                     "New ID": mid})
    return rows


@app.callback(
    Output("tbl-met", "data", allow_duplicate=True),
    Output("tbl-rxn", "data", allow_duplicate=True),
    Output("rename-messages", "children"),
    Input("btn-apply-rename", "n_clicks"),
    State("tbl-rename", "data"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def apply_custom_ids(_n, rename_rows, met_data, rxn_data):
    if not rename_rows:
        return no_update, no_update, "Click 'Load current metabolites' first."
    id_map = {}
    for row in rename_rows:
        cur = str(row.get("Current ID", "")).strip()
        new = str(row.get("New ID", "")).strip()
        if cur and new and new != cur:
            id_map[cur] = new
    if not id_map:
        return no_update, no_update, "No id changes to apply."

    # Validate: the resulting ids must stay unique (no two metabolites collide).
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    current_ids = [str(x).strip() for x in df_met["ID"].dropna()]
    resulting = [id_map.get(i, i) for i in current_ids]
    dupes = sorted({x for x in resulting if resulting.count(x) > 1})
    if dupes:
        return no_update, no_update, dbc.Alert(
            f"⚠ These ids would collide — pick distinct New IDs: "
            f"{', '.join(dupes)}", color="danger")

    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    m2, r2 = idmap.apply_id_map(df_met, df_rxn, id_map)
    changes = ", ".join(f"{k}→{v}" for k, v in id_map.items())
    return (df_to_records(m2), df_to_records(r2),
            dbc.Alert(f"Renamed {len(id_map)} metabolite(s): {changes}",
                      color="success"))


@app.callback(
    Output("tbl-rxn-rename", "data"),
    Input("btn-load-rxn-rename", "n_clicks"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def populate_rxn_rename(_n, rxn_data):
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    rows = []
    for _, r in df_rxn.iterrows():
        rid = str(r.get("ID", "")).strip()
        if not rid:
            continue
        rows.append({"KEGG Reaction ID": r.get("KEGG Reaction ID"),
                     "Current ID": rid,
                     "Stoichiometry": r.get("Reaction stoichiometry"),
                     "New ID": rid})
    return rows


@app.callback(
    Output("tbl-rxn", "data", allow_duplicate=True),
    Output("rxn-rename-messages", "children"),
    Input("btn-apply-rxn-rename", "n_clicks"),
    State("tbl-rxn-rename", "data"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def apply_rxn_custom_ids(_n, rename_rows, rxn_data):
    if not rename_rows:
        return no_update, "Click 'Load current reactions' first."
    id_map = {}
    for row in rename_rows:
        cur = str(row.get("Current ID", "")).strip()
        new = str(row.get("New ID", "")).strip()
        if cur and new and new != cur:
            id_map[cur] = new
    if not id_map:
        return no_update, "No id changes to apply."

    # Validate: the resulting ids must stay unique (no two reactions collide).
    # Unlike metabolites, reaction ids aren't referenced anywhere else in the
    # tables, so this is a plain rename — no equation rewriting needed.
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    current_ids = [str(x).strip() for x in df_rxn["ID"].dropna()]
    resulting = [id_map.get(i, i) for i in current_ids]
    dupes = sorted({x for x in resulting if resulting.count(x) > 1})
    if dupes:
        return no_update, dbc.Alert(
            f"⚠ These ids would collide — pick distinct New IDs: "
            f"{', '.join(dupes)}", color="danger")

    df_rxn["ID"] = df_rxn["ID"].map(lambda x: id_map.get(str(x).strip(), x))
    changes = ", ".join(f"{k}→{v}" for k, v in id_map.items())
    return (df_to_records(df_rxn),
            dbc.Alert(f"Renamed {len(id_map)} reaction(s): {changes}",
                      color="success"))


_ARROWS = {"<=>", "->", "=>", "<->", "→", "<-", "<="}


def _referenced_ids(df_rxn):
    """Working metabolite ids that appear in any reaction stoichiometry."""
    used = set()
    for eq in df_rxn["Reaction stoichiometry"].dropna().astype(str):
        for tok in eq.replace("+", " ").split():
            if tok in _ARROWS:
                continue
            try:
                float(tok)  # a stoichiometric coefficient, not a metabolite
                continue
            except ValueError:
                pass
            used.add(tok)
    return used


@app.callback(
    Output("remove-rxns", "options"),
    Input("tbl-rxn", "data"),
)
def update_remove_options(rxn_data):
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    return [{"label": r, "value": r}
            for r in df_rxn["ID"].dropna().astype(str) if r.strip()]


@app.callback(
    Output("tbl-rxn", "data", allow_duplicate=True),
    Output("remove-messages", "children"),
    Input("btn-remove-rxns", "n_clicks"),
    State("remove-rxns", "value"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def remove_reactions(_n, selected, rxn_data):
    if not selected:
        return no_update, "Select at least one reaction to remove."
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    drop = set(selected)
    kept = df_rxn[~df_rxn["ID"].astype(str).isin(drop)]
    n = len(df_rxn) - len(kept)
    return (df_to_records(kept),
            dbc.Alert(f"Removed {n} reaction(s): {', '.join(sorted(drop))}. "
                      f"Metabolites are kept — use 'Remove unused metabolites' "
                      f"to prune any orphans.", color="success"))


@app.callback(
    Output("tbl-met", "data", allow_duplicate=True),
    Output("remove-messages", "children", allow_duplicate=True),
    Input("btn-prune-mets", "n_clicks"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def prune_metabolites(_n, met_data, rxn_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    used = _referenced_ids(df_rxn)
    mask_keep = df_met["ID"].astype(str).str.strip().isin(used)
    removed = [str(x) for x in df_met.loc[~mask_keep, "ID"].dropna()]
    if not removed:
        return no_update, "No unused metabolites — nothing to remove."
    return (df_to_records(df_met[mask_keep]),
            dbc.Alert(f"Removed {len(removed)} unused metabolite(s): "
                      f"{', '.join(removed)}.", color="success"))


@app.callback(
    Output("reversible-rxns", "options"),
    Input("tbl-rxn", "data"),
)
def update_reversible_options(rxn_data):
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    return [{"label": r, "value": r}
            for r in df_rxn["ID"].dropna().astype(str) if r.strip()]


@app.callback(
    Output("tbl-rxn", "data", allow_duplicate=True),
    Output("reversibility-messages", "children"),
    Input("btn-mark-reversible", "n_clicks"),
    Input("btn-mark-irreversible", "n_clicks"),
    State("reversible-rxns", "value"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def mark_reversibility(_n_rev, _n_irrev, selected, rxn_data):
    if not selected:
        return no_update, "Select at least one reaction first."
    trig = ctx.triggered_id
    new_value = 1 if trig == "btn-mark-reversible" else 0
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    mask = df_rxn["ID"].astype(str).isin(set(selected))
    df_rxn.loc[mask, "Reversibility"] = new_value
    label = "reversible" if new_value else "irreversible"
    return (df_to_records(df_rxn),
            dbc.Alert(f"Marked {mask.sum()} reaction(s) {label}: "
                      f"{', '.join(sorted(selected))}.", color="success"))


@app.callback(
    Output("exchange-mets", "options"),
    Input("tbl-met", "data"),
    Input("tbl-rxn", "data"),
)
def update_exchange_options(met_data, rxn_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    have_ex = {str(r)[2:] for r in df_rxn["ID"].dropna().astype(str)
               if str(r).startswith("EX")}
    opts = []
    for mid in df_met["ID"].dropna().astype(str):
        mid = mid.strip()
        if not mid or mid.endswith("ex") or mid in have_ex:
            continue
        opts.append({"label": mid, "value": mid})
    return opts


@app.callback(
    Output("tbl-met", "data", allow_duplicate=True),
    Output("tbl-rxn", "data", allow_duplicate=True),
    Output("exchange-messages", "children"),
    Input("btn-exchange", "n_clicks"),
    State("exchange-mets", "value"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def add_exchanges(_n, selected, met_data, rxn_data):
    if not selected:
        return no_update, no_update, "Select at least one boundary metabolite."
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    m2, r2, messages = exchanges.add_exchange_transport(df_met, df_rxn, selected)
    msg = html.Ul([html.Li(m) for m in messages]) if messages else "Nothing added."
    return df_to_records(m2), df_to_records(r2), msg


@app.callback(
    Output("connectivity-output", "children"),
    Input("btn-connectivity", "n_clicks"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def check_connectivity(_n, met_data, rxn_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if df_rxn.empty or df_rxn["ID"].dropna().empty:
        return dbc.Alert("No reactions yet — fetch some from KEGG above.",
                         color="warning")
    try:
        blocked = model_mod.find_blocked_reactions(df_met, df_rxn,
                                                    "connectivity_check")
    except Exception as e:
        return dbc.Alert(f"⚠ Could not build the model: {e}", color="danger")

    if not blocked:
        return dbc.Alert(
            "✓ Every reaction can carry flux — no dead ends given the "
            "reactions and exchanges added so far.", color="success")

    blocked_ex = sorted(r for r in blocked if str(r).startswith("EX"))
    blocked_other = sorted(r for r in blocked if not str(r).startswith("EX"))
    parts = [html.B(f"{len(blocked)} reaction(s) can never carry flux "
                    f"(default bounds — no substrate/product/energy-product "
                    f"choices applied yet):")]
    if blocked_ex:
        parts.append(html.Div([
            "Blocked exchanges — usually a missing transport reaction, or "
            "this boundary metabolite isn't reachable from anything else "
            "in the network yet: ",
            html.Code(", ".join(blocked_ex)),
        ], className="mt-1"))
    if blocked_other:
        parts.append(html.Div([
            "Other blocked reactions — check that each of their metabolites "
            "is both produced and consumed somewhere (a missing exchange "
            "for a byproduct is a common cause): ",
            html.Code(", ".join(blocked_other)),
        ], className="mt-1"))
    return dbc.Alert(parts, color="warning")


@app.callback(
    Output("balance-panel", "children"),
    Output("fix-reaction", "options"),
    Output("fix-reaction", "value"),
    Input("tbl-met", "data"),
    Input("tbl-rxn", "data"),
    State("fix-reaction", "value"),
)
def update_balance(met_data, rxn_data, fix_value):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if df_rxn.empty or df_rxn["ID"].dropna().empty:
        return html.Div("No reactions yet — fetch some from KEGG above."), [], None

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
        "Equation (working ids)": df_bal["equation_custom"],
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
        page_size=100,
        style_table={"overflowX": "auto", "maxHeight": "360px", "overflowY": "auto"},
    )
    if ok:
        banner = dbc.Alert("✓ All non-exchange reactions are atom- and "
                           "charge-balanced — ready to download.", color="success")
    else:
        items = [html.Li(f"{e['reaction_id']}: {e['imbalance'] or e['notes']}")
                 for e in errors]
        banner = dbc.Alert([html.B(f"{len(errors)} imbalanced reaction(s):"),
                            html.Ul(items)], color="danger")

    flagged = [{"label": e["reaction_id"], "value": e["reaction_id"]} for e in errors]
    # Keep the current pick if it's still imbalanced; otherwise advance to the
    # first remaining flagged reaction, so the fix buttons never act on a stale
    # selection.
    flagged_ids = [e["reaction_id"] for e in errors]
    fix_val = fix_value if fix_value in flagged_ids else (
        flagged_ids[0] if flagged_ids else None)
    return html.Div([banner, table]), flagged, fix_val


@app.callback(
    Output("tbl-rxn", "data", allow_duplicate=True),
    Input("fix-h-left", "n_clicks"),
    Input("fix-h-right", "n_clicks"),
    Input("fix-w-left", "n_clicks"),
    Input("fix-w-right", "n_clicks"),
    State("fix-reaction", "value"),
    State("tbl-rxn", "data"),
    State("tbl-met", "data"),
    prevent_initial_call=True,
)
def apply_fix(_hl, _hr, _wl, _wr, rxn_id, rxn_data, met_data):
    if not rxn_id:
        return no_update
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
        return no_update
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
    return df_to_records(df_rxn)


@app.callback(
    Output("download-xlsx", "data"),
    Input("btn-download", "n_clicks"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    State("wb-name", "value"),
    prevent_initial_call=True,
)
def download_workbook(_n, met_data, rxn_data, name):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    data = io.save_excel(df_met, df_rxn)
    return dcc.send_bytes(data, f"{name or 'kegg_workbook'}.xlsx")


if __name__ == "__main__":
    app.run(debug=False, host="127.0.0.1", port=8051)
