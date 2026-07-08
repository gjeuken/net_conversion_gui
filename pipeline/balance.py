"""Formula/charge parsing and atom + charge balance checks.

Ported verbatim (logic-wise) from the notebook's balance cell.  These are pure
functions: equation strings + a metabolite mapping table in, balance results
out.  KEGG REST access lives in :mod:`pipeline.kegg`; the ``analyze_*``
functions below pull compound formula/charge through it (lazily imported to
avoid a circular dependency).
"""

from __future__ import annotations

import re
from collections import Counter

import pandas as pd

# -------- ELEMENT SYMBOL WHITELIST (up to Og) --------
ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
}

# Coefficient tokens KEGG uses for polymers / generic R-groups that have no
# numeric stoichiometry.  Surfaced (not guessed) per CLAUDE.md §8.
_NON_NUMERIC_COEFF = re.compile(r"^\(?n([+\-]\d+)?\)?$|^m$", re.IGNORECASE)


# -------- Equation parsing (robust to common arrows) --------
def normalize_arrow(eq: str) -> str:
    """Normalise the various reaction arrows to ``<=>``."""
    return re.sub(r"\s*(<[-=]+>|<=>|<->|⇌|↔|→|->|=>)\s*", " <=> ", eq)


def parse_equation(equation: str):
    """Parse a reaction equation into ``(substrates, products)``.

    Each side is a list of ``(coeff, compound_id)`` with float coefficients.
    Supports integer and fractional coefficients.
    """
    equation = normalize_arrow(equation)
    if "<=>" not in equation:
        raise ValueError(f"No recognized arrow in equation: {equation}")
    left, right = equation.split("<=>")

    def parse_side(side):
        parts = [p.strip() for p in side.strip().split("+") if p.strip()]
        parsed = []
        for p in parts:
            m = re.match(r"^\s*([\d\.]+)\s+(\S+)\s*$", p)
            if m:
                coeff = float(m.group(1))
                cid = m.group(2)
            else:
                coeff, cid = 1.0, p.strip()
            parsed.append((coeff, cid))
        return parsed

    return parse_side(left), parse_side(right)


def has_nonnumeric_coefficient(equation: str) -> bool:
    """True if any term uses an R-group/polymer coefficient (n, n+1, (n), m)."""
    try:
        equation = normalize_arrow(equation)
        sides = equation.split("<=>")
    except Exception:
        return False
    for side in sides:
        for p in side.split("+"):
            p = p.strip()
            if not p:
                continue
            tok = p.split()
            if len(tok) >= 2 and _NON_NUMERIC_COEFF.match(tok[0]):
                return True
    return False


# -------- Formula parsing (only real elements, supports hydrates) --------
def parse_formula(formula) -> Counter:
    """Parse a chemical formula string into element counts."""
    if not formula:
        return Counter()

    counts = Counter()
    parts = re.split(r"[·\.]\s*|\s{2,}", str(formula).strip())
    for part in parts:
        if not part:
            continue
        part = re.sub(r"[^A-Za-z0-9]", "", part)
        i = 0
        while i < len(part):
            if i + 1 < len(part) and part[i:i + 2] in ELEMENTS:
                sym = part[i:i + 2]
                i += 2
            elif part[i] in [s[:1] for s in ELEMENTS]:
                sym1 = part[i]
                i += 1
                if sym1 in ELEMENTS:
                    sym = sym1
                else:
                    continue
            else:
                i += 1
                continue

            j = i
            while j < len(part) and part[j].isdigit():
                j += 1
            num = int(part[i:j]) if j > i else 1
            counts[sym] += num
            i = j
    return counts


def multiply_formula(counts: Counter, coeff: float) -> Counter:
    return Counter({atom: n * coeff for atom, n in counts.items()})


# -------- User-supplied formulas/charges (for non-KEGG metabolites) --------
def formula_overrides(mapping_df: pd.DataFrame):
    """Pull explicit ``Chemical formula`` / ``Charge`` from the metabolite table.

    Metabolites without a KEGG ID can't be looked up, but the user may supply a
    ``Chemical formula`` column (and optionally ``Charge``) directly.  When
    present, that value wins over the KEGG lookup — it's the escape hatch for
    designed/non-KEGG compounds (CLAUDE.md §5, §8).

    Returns ``(formula_map, charge_map)`` keyed by the id each metabolite carries
    *after* mapping to KEGG — i.e. the KEGG ID when one exists, otherwise the
    working ID (which is exactly what ``map_custom_to_kegg_equation`` falls back
    to), so the keys line up with the equations passed to :func:`check_balance`.
    """
    formula_map, charge_map = {}, {}
    has_formula = "Chemical formula" in mapping_df.columns
    has_charge = "Charge" in mapping_df.columns
    if not (has_formula or has_charge):
        return formula_map, charge_map

    for _, row in mapping_df.iterrows():
        wid = row.get("ID")
        kid = row.get("KEGG ID")
        mapped = kid if isinstance(kid, str) and kid.strip() else wid
        if not (isinstance(mapped, str) and mapped.strip()):
            continue
        mapped = mapped.strip()
        if has_formula:
            f = row.get("Chemical formula")
            if isinstance(f, str) and f.strip():
                formula_map[mapped] = f.strip()
        if has_charge:
            c = row.get("Charge")
            if pd.notna(c) and str(c).strip() != "":
                try:
                    charge_map[mapped] = int(float(c))
                except (ValueError, TypeError):
                    pass
    return formula_map, charge_map


# -------- Clean mapping: custom IDs -> KEGG IDs by token --------
def _clean_kegg(kid, fallback):
    """A KEGG id, or ``fallback`` when it's missing/blank (NaN, None, "")."""
    if kid is None or (isinstance(kid, float) and pd.isna(kid)):
        return str(fallback)
    if isinstance(kid, str) and not kid.strip():
        return str(fallback)
    return str(kid)


def map_custom_to_kegg_equation(equation: str, mapping_df: pd.DataFrame):
    custom2kegg = dict(zip(mapping_df["ID"], mapping_df["KEGG ID"]))
    subs, prods = parse_equation(equation)

    def map_side(pairs):
        mapped = []
        for coeff, cid in pairs:
            kid = _clean_kegg(custom2kegg.get(cid, cid), cid)
            mapped.append((coeff, kid))
        return mapped

    subs_k = map_side(subs)
    prods_k = map_side(prods)

    def side_to_str(pairs):
        parts = []
        for coeff, cid in pairs:
            if abs(coeff - 1.0) < 1e-12:
                parts.append(cid)
            else:
                parts.append(f"{coeff:g} {cid}")
        return " + ".join(parts) if parts else "0"

    eq_kegg = f"{side_to_str(subs_k)} <=> {side_to_str(prods_k)}"
    return eq_kegg, subs_k, prods_k


# -------- Balance check with signed atom deltas --------
def check_balance(substrates, products, formula_map, charge_map, tol=1e-6):
    left_atoms, right_atoms = Counter(), Counter()
    left_charge, right_charge = 0.0, 0.0

    for coeff, cid in substrates:
        f = formula_map.get(cid)
        if f:
            left_atoms += multiply_formula(parse_formula(f), coeff)
        left_charge += charge_map.get(cid, 0) * coeff

    for coeff, cid in products:
        f = formula_map.get(cid)
        if f:
            right_atoms += multiply_formula(parse_formula(f), coeff)
        right_charge += charge_map.get(cid, 0) * coeff

    all_atoms = set(left_atoms) | set(right_atoms)
    delta = {
        a: right_atoms.get(a, 0.0) - left_atoms.get(a, 0.0) for a in all_atoms
    }

    atoms_bal = all(abs(d) < tol for d in delta.values())
    charge_bal = abs(left_charge - right_charge) < tol

    imbalance_summary = ", ".join(
        f"{'+' if d > 0 else ''}{d:g} {a}"
        for a, d in sorted(delta.items())
        if abs(d) >= tol
    ) or None
    charge_delta = right_charge - left_charge
    if abs(charge_delta) >= tol:
        chg = f"{'+' if charge_delta > 0 else ''}{charge_delta:g} charge"
        imbalance_summary = (
            f"{imbalance_summary}, {chg}" if imbalance_summary else chg
        )

    return (atoms_bal, charge_bal, left_atoms, right_atoms, left_charge,
            right_charge, delta, imbalance_summary)


# -------- Pipeline: custom eq -> KEGG eq -> balance --------
def analyze_custom_equations(equations, mapping_df, reaction_ids=None):
    """Map custom equations to KEGG ids and balance-check each.

    ``reaction_ids`` (optional, parallel to ``equations``) is carried through so
    callers can flag exchange reactions, which are expected to be unbalanced.
    """
    from .kegg import get_compound_info  # lazy: avoid circular import

    f_over, c_over = formula_overrides(mapping_df)

    rows = []
    for idx, eq in enumerate(equations):
        rid = reaction_ids[idx] if reaction_ids is not None else None
        if eq is None or (isinstance(eq, float) and pd.isna(eq)) or str(eq).strip() == "":
            rows.append(_blank_row(rid, eq, "Empty equation"))
            continue
        if has_nonnumeric_coefficient(str(eq)):
            rows.append(_blank_row(
                rid, eq, "R-group / polymer coefficient — enter manually"))
            continue
        try:
            eq_kegg, subs_k, prods_k = map_custom_to_kegg_equation(str(eq), mapping_df)
        except Exception as e:
            rows.append(_blank_row(rid, eq, f"Parse error: {e}"))
            continue

        cids = [cid for _, cid in subs_k + prods_k]
        uniq = list(dict.fromkeys(cids))
        # KEGG lookup first, then let any user-supplied formula/charge win.
        compound_info = {cid: get_compound_info(cid) for cid in uniq}
        formula_map = {cid: f_over.get(cid, v[0]) for cid, v in compound_info.items()}
        charge_map = {cid: c_over.get(cid, v[1]) for cid, v in compound_info.items()}

        (atoms_bal, charge_bal, _la, _ra, left_charge, right_charge, delta,
         imbalance_summary) = check_balance(subs_k, prods_k, formula_map, charge_map)

        is_exchange = bool(rid) and str(rid).startswith("EX")
        notes = []
        if not atoms_bal:
            notes.append("Atoms not balanced")
        if not charge_bal:
            notes.append("Charge not balanced")
        if not notes:
            notes.append("OK")
        elif is_exchange:
            notes.append("(exchange — expected)")

        rows.append({
            "reaction_id": rid,
            "equation_custom": eq,
            "equation_kegg": eq_kegg,
            "atom_balanced": atoms_bal,
            "charge_balanced": charge_bal,
            "is_exchange": is_exchange,
            "left_charge": left_charge,
            "right_charge": right_charge,
            "atom_delta": delta,
            "atom_imbalance": imbalance_summary,
            "notes": "; ".join(notes),
        })

    return pd.DataFrame(rows)


def _blank_row(rid, eq, note):
    return {
        "reaction_id": rid,
        "equation_custom": eq,
        "equation_kegg": None,
        "atom_balanced": None,
        "charge_balanced": None,
        "is_exchange": bool(rid) and str(rid).startswith("EX"),
        "left_charge": None,
        "right_charge": None,
        "atom_delta": None,
        "atom_imbalance": None,
        "notes": note,
    }


def analyze_stoichiometry_matrix(stoich_df, mapping_df, formula_map=None,
                                 charge_map=None, tol=1e-6, only_unbalanced=False):
    """Balance-check reactions given as a metabolite x reaction stoich matrix.

    Rows are metabolites with an ``EX``-prefixed index (abbrev = ``index[2:]``);
    columns are reactions; values are signed coefficients.
    """
    from .kegg import get_compound_info  # lazy

    f_over, c_over = formula_overrides(mapping_df)

    rows = []
    abbrev2kegg = dict(zip(mapping_df["ID"], mapping_df["KEGG ID"]))

    for rxn in stoich_df.columns:
        coeffs = stoich_df[rxn].dropna()
        substrates, products = [], []

        for met_full, coeff in coeffs.items():
            if abs(coeff) < tol:
                continue
            abbrev = met_full[2:]
            kegg_id = _clean_kegg(abbrev2kegg.get(abbrev, abbrev), abbrev)
            if coeff < 0:
                substrates.append((-coeff, kegg_id))
            else:
                products.append((coeff, kegg_id))

        def side_to_str(pairs):
            parts = []
            for coeff, cid in pairs:
                if abs(coeff - 1.0) < tol:
                    parts.append(cid)
                else:
                    parts.append(f"{coeff:g} {cid}")
            return " + ".join(parts) if parts else "0"

        eq_kegg = f"{side_to_str(substrates)} <=> {side_to_str(products)}"

        if formula_map is None or charge_map is None:
            cids = [cid for _, cid in substrates + products]
            uniq = list(dict.fromkeys(cids))
            compound_info = {cid: get_compound_info(cid) for cid in uniq}
            f_map = {cid: f_over.get(cid, v[0]) for cid, v in compound_info.items()}
            c_map = {cid: c_over.get(cid, v[1]) for cid, v in compound_info.items()}
        else:
            f_map, c_map = formula_map, charge_map

        (atoms_bal, charge_bal, _la, _ra, left_charge, right_charge, delta,
         imbalance_summary) = check_balance(substrates, products, f_map, c_map, tol=tol)

        if only_unbalanced and atoms_bal and charge_bal:
            continue

        notes = []
        if not atoms_bal:
            notes.append("Atoms not balanced")
        if not charge_bal:
            notes.append("Charge not balanced")
        if not notes:
            notes.append("OK")

        rows.append({
            "reaction_id": rxn,
            "equation_kegg": eq_kegg,
            "atom_balanced": atoms_bal,
            "charge_balanced": charge_bal,
            "left_charge": left_charge,
            "right_charge": right_charge,
            "atom_delta": delta,
            "atom_imbalance": imbalance_summary,
            "notes": "; ".join(notes),
        })

    return pd.DataFrame(rows)
