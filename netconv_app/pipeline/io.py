"""Input/output: Excel save/load and the canonical-dataframe helpers.

The two canonical dataframes (CLAUDE.md §5) are the single source of truth:

* **Metabolites** — columns ``ID``, ``Name``, ``KEGG ID``
* **Reactions**   — columns ``ID``, ``Name``, ``Reaction stoichiometry``,
  ``Reversibility`` (1 = reversible, 0 = irreversible)
"""

from __future__ import annotations

import io as _io

import pandas as pd

METABOLITE_COLS = ["ID", "Name", "KEGG ID"]
REACTION_COLS = ["ID", "Name", "Reaction stoichiometry", "Reversibility"]

# Run-time selectors (the notebook's USER INPUTS cell).  Stored in an optional
# "Config" sheet so an example workbook carries its own modelling choices.
CONFIG_KEYS = ["MODEL_NAME", "SUBSTRATES", "PRODUCTS", "CARBON_PRODUCTS",
               "ENERGY_PRODUCT", "REV_ALLOWED", "PSEUDO_METS"]
_LIST_KEYS = {"SUBSTRATES", "PRODUCTS", "CARBON_PRODUCTS", "REV_ALLOWED",
              "PSEUDO_METS"}


def clean_dataframe_whitespace(df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
    """Strip newline/tab characters from every string cell."""
    target = df if inplace else df.copy()
    target = target.apply(
        lambda col: col.map(
            lambda x: x.replace("\n", "").replace("\t", "") if isinstance(x, str) else x
        )
    )
    return target


def empty_tables():
    """Return a pair of empty canonical dataframes with the right columns."""
    return (pd.DataFrame(columns=METABOLITE_COLS),
            pd.DataFrame(columns=REACTION_COLS))


def _coerce(df: pd.DataFrame, cols) -> pd.DataFrame:
    """Ensure all required columns exist and are in canonical order."""
    df = df.copy()
    for c in cols:
        if c not in df.columns:
            df[c] = pd.Series(dtype=object)
    return df[cols]


def load_excel(source):
    """Load Metabolites + Reactions sheets from an Excel workbook.

    ``source`` may be a path or a file-like / bytes object (Dash uploads).
    """
    if isinstance(source, (bytes, bytearray)):
        source = _io.BytesIO(source)
    df_met = pd.read_excel(source, sheet_name="Metabolites")
    if isinstance(source, _io.BytesIO):
        source.seek(0)
    df_rxn = pd.read_excel(source, sheet_name="Reactions")

    df_met = clean_dataframe_whitespace(_coerce(df_met, METABOLITE_COLS))
    df_rxn = clean_dataframe_whitespace(_coerce(df_rxn, REACTION_COLS))
    # Reversibility defaults to reversible if blank.
    df_rxn["Reversibility"] = (
        pd.to_numeric(df_rxn["Reversibility"], errors="coerce").fillna(1).astype(int)
    )
    return df_met, df_rxn


def save_excel(df_metabolites, df_reactions, path=None, config=None) -> bytes:
    """Write the two canonical dataframes (and optional config) to an .xlsx.

    If ``path`` is given the file is written there; the raw bytes are always
    returned (handy for Dash download callbacks).
    """
    buf = _io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xl:
        _coerce(df_metabolites, METABOLITE_COLS).to_excel(
            xl, sheet_name="Metabolites", index=False)
        _coerce(df_reactions, REACTION_COLS).to_excel(
            xl, sheet_name="Reactions", index=False)
        if config:
            _config_to_df(config).to_excel(xl, sheet_name="Config", index=False)
    data = buf.getvalue()
    if path is not None:
        with open(path, "wb") as fh:
            fh.write(data)
    return data


def _config_to_df(config) -> pd.DataFrame:
    rows = []
    for k in CONFIG_KEYS:
        if k not in config:
            continue
        v = config[k]
        if isinstance(v, (list, tuple, set)):
            v = ", ".join(str(x) for x in v)
        rows.append({"Key": k, "Value": v})
    return pd.DataFrame(rows, columns=["Key", "Value"])


def load_config(source):
    """Read the optional Config sheet into a dict, or return ``None``.

    List-valued keys are split on commas; ``PSEUDO_METS`` becomes a set.
    """
    if isinstance(source, (bytes, bytearray)):
        source = _io.BytesIO(source)
    try:
        df = pd.read_excel(source, sheet_name="Config")
    except (ValueError, KeyError):
        return None
    if "Key" not in df.columns or "Value" not in df.columns:
        return None
    config = {}
    for _, row in df.iterrows():
        key = str(row["Key"]).strip()
        val = row["Value"]
        if key in _LIST_KEYS:
            items = [x.strip() for x in str(val).split(",") if str(x).strip()]
            config[key] = set(items) if key == "PSEUDO_METS" else items
        else:
            config[key] = None if pd.isna(val) else str(val)
    return config


def dict_to_kegg_reaction(met_dict) -> str:
    """Convert a ``{kegg_id: signed_coeff}`` dict into a KEGG reaction string.

    Negative coefficients are reactants, positive are products — the format
    eQuilibrator's ``parse_reaction_formula`` expects.
    """
    reactants, products = [], []
    for met, coeff in met_dict.items():
        if coeff < 0:
            term = f"{-coeff:g} kegg:{met}" if abs(coeff) != 1 else f"kegg:{met}"
            reactants.append(term)
        elif coeff > 0:
            term = f"{coeff:g} kegg:{met}" if abs(coeff) != 1 else f"kegg:{met}"
            products.append(term)
    return f"{' + '.join(reactants)} = {' + '.join(products)}"
