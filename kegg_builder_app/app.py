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
from dash import (Dash, Input, Output, State, ctx, dash_table, dcc, html,
                  no_update)

from pipeline import balance, bigg, exchanges, idmap, io, kegg, run
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
        html.H5("1a · Add a reaction manually (non-KEGG)"),
        html.P("For transport, PTS import, ion translocation or designed "
               "reactions KEGG has no entry for. Write the stoichiometry in "
               "working ids with any arrow (-> or <=>). Any metabolite not yet "
               "in the table is added as a placeholder for you to complete "
               "(KEGG id and/or chemical formula, so it can balance).",
               className="text-muted small"),
        dbc.Row([
            dbc.Col(dbc.Input(id="man-id", placeholder="Reaction id (e.g. PTS_GLC)"),
                    md=3),
            dbc.Col(dbc.Input(id="man-name", placeholder="Name (optional)"), md=3),
            dbc.Col(dbc.Input(id="man-eq",
                              placeholder="e.g. GLC + PEP -> G6P + PYR"), md=4),
            dbc.Col(dbc.Select(id="man-rev", value="1",
                               options=[{"label": "reversible", "value": "1"},
                                        {"label": "irreversible", "value": "0"}]),
                    md=1),
            dbc.Col(dbc.Button("Add", id="btn-add-manual", color="primary"), md=1),
        ], className="g-2"),
        html.Div(id="manual-messages",
                 style={"fontSize": "12px", "marginTop": "0.5rem"}),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("1b · Human-readable ids (optional)"),
        html.P("Translate KEGG compound ids to BiGG ids (e.g. C00031 → glc__D). "
               "The KEGG ID column is kept, so balance and thermodynamics still "
               "resolve compounds correctly. Compounds with no clean BiGG id "
               "fall back to a slugified KEGG name and are flagged for you to "
               "hand-edit. First use downloads a small BiGG file (needs network).",
               className="text-muted small"),
        dbc.Button("Translate KEGG → BiGG ids", id="btn-bigg",
                   color="secondary", outline=True),
        dcc.Loading(html.Div(id="bigg-messages",
                             style={"maxHeight": "220px", "overflowY": "auto",
                                    "fontSize": "12px", "marginTop": "0.5rem"})),
    ]), className="mb-3"),

    dbc.Card(dbc.CardBody([
        html.H5("1c · Custom metabolite ids (optional)"),
        html.P("Set your own working ids. Click Load to fill the table from the "
               "current metabolites, edit the New ID column, then Apply — the "
               "reaction stoichiometries and EX/transport ids are rewritten to "
               "match, and the KEGG ID column is kept. Do this before adding "
               "exchanges so those pick up your ids.", className="text-muted small"),
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

    dbc.Row([
        dbc.Col([
            html.H5("Reactions"),
            dbc.Button("+ Add reaction row", id="btn-add-rxn", size="sm",
                       color="light", className="mb-1"),
            _table("tbl-rxn", io.REACTION_COLS),
        ], md=7),
        dbc.Col([
            html.H5("Metabolites"),
            dbc.Button("+ Add metabolite row", id="btn-add-met", size="sm",
                       color="light", className="mb-1"),
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
        html.H5("1d · Exchange & transport reactions"),
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
        html.H5("1e · Network connectivity check"),
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
    prevent_initial_call=True,
)
def translate_bigg(_n, met_data, rxn_data):
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if df_met.empty or df_met["ID"].dropna().empty:
        return no_update, no_update, "Fetch some reactions first."
    try:
        m2, r2, report = bigg.translate_to_bigg(df_met, df_rxn)
    except bigg.BiggError as e:
        return no_update, no_update, dbc.Alert(
            f"⚠ {e} — needs network on first use.", color="warning")

    n = len(report)
    clean = sum(1 for x in report if x["status"] == "bigg")
    flags = [x for x in report if x["status"] != "bigg"]
    summary = dbc.Alert(
        f"Translated {clean}/{n} metabolites to clean BiGG ids"
        + (f"; {len(flags)} need a look." if flags else "."),
        color="success" if not flags else "info")
    detail = html.Ul([
        html.Li([html.Code(f"{x['kegg'] or x['old']} → {x['new']}"),
                 f"  ({_BIGG_STATUS_LABEL.get(x['status'], x['status'])})"])
        for x in flags]) if flags else None
    return (df_to_records(m2), df_to_records(r2),
            html.Div([summary, detail] if detail else [summary]))


@app.callback(
    Output("tbl-met", "data", allow_duplicate=True),
    Output("tbl-rxn", "data", allow_duplicate=True),
    Output("manual-messages", "children"),
    Input("btn-add-manual", "n_clicks"),
    State("man-id", "value"),
    State("man-name", "value"),
    State("man-eq", "value"),
    State("man-rev", "value"),
    State("tbl-met", "data"),
    State("tbl-rxn", "data"),
    prevent_initial_call=True,
)
def add_manual_reaction(_n, rid, name, eq, rev, met_data, rxn_data):
    rid = (rid or "").strip()
    eq = (eq or "").strip()
    if not rid or not eq:
        return no_update, no_update, "Enter a reaction id and stoichiometry."
    df_met = records_to_df(met_data, io.METABOLITE_COLS)
    df_rxn = records_to_df(rxn_data, io.REACTION_COLS)
    if (df_rxn["ID"].astype(str).str.strip() == rid).any():
        return no_update, no_update, dbc.Alert(
            f"⚠ Reaction '{rid}' already exists.", color="warning")
    try:
        subs, prods = balance.parse_equation(eq)
    except Exception as e:
        return no_update, no_update, dbc.Alert(
            f"⚠ Could not parse the stoichiometry: {e}", color="danger")

    df_rxn = pd.concat([df_rxn, pd.DataFrame([{
        "ID": rid, "Name": (name or "").strip() or rid,
        "Reaction stoichiometry": balance.normalize_arrow(eq),
        "Reversibility": int(rev) if rev in ("0", "1") else 1,
    }])], ignore_index=True)

    known = {str(x).strip() for x in df_met["ID"].dropna()}
    new_mets = []
    for _c, cid in subs + prods:
        cid = str(cid).strip()
        if cid and cid not in known:
            new_mets.append({"ID": cid, "Name": cid, "KEGG ID": "",
                             "Chemical formula": ""})
            known.add(cid)
    if new_mets:
        df_met = pd.concat([df_met, pd.DataFrame(new_mets)], ignore_index=True)

    parts = [f"✓ Added reaction {rid}"]
    if new_mets:
        parts.append("new metabolite placeholder(s): "
                     + ", ".join(m["ID"] for m in new_mets)
                     + " — add a KEGG id or chemical formula so they balance.")
    return (df_to_records(df_met), df_to_records(df_rxn),
            dbc.Alert(" — ".join(parts), color="success"))


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
