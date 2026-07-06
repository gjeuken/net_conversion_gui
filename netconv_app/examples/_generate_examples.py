"""Generate the bundled example workbooks.

Run once from the project root:

    python netconv_app/examples/_generate_examples.py

Produces, in this folder:
  * Example1_EMP_lactate.xlsx     — default; reproduces the paper's headline
                                     1 GLC + 2 ADP + 2 Pi -> 2 LAC + 2 ATP
  * Example1_EMPglycolysis.xlsx   — the notebook's ethanol/acetate/formate net
  * pan_glycolysis.xlsx           — copied from the repo root (+ Config sheet)
"""

import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))  # project root

from netconv_app.pipeline import io  # noqa: E402

# Master KEGG mapping for every metabolite used below (working ID -> KEGG id).
MET_KEGG = {
    "GLCex": ("Glucose, extracellular", "C00031"),
    "GLC": ("Glucose, intracellular", "C00031"),
    "G6P": ("Glucose-6-phosphate", "C00092"),
    "F6P": ("Fructose-6-phosphate", "C00085"),
    "FBP": ("Fructose-1,6-bisphosphate", "C00354"),
    "G3P": ("Glyceraldehyde-3-phosphate", "C00118"),
    "DHAP": ("Dihydroxyacetonephosphate", "C00111"),
    "BPG": ("1,3-bisphosphoglycerate", "C00236"),
    "PG3": ("3-phosphoglycerate", "C00197"),
    "PG2": ("2-phosphoglycerate", "C00631"),
    "PEP": ("phosphoenolpyruvate", "C00074"),
    "PYR": ("pyruvate", "C00022"),
    "LAC": ("lactate", "C00256"),
    "LACex": ("Lactate, extracellular", "C00256"),
    "FOR": ("Formate", "C00058"),
    "ACCOA": ("acetyl-CoA", "C00024"),
    "ACTP": ("acetyl-phosphate", "C00227"),
    "AC": ("acetate", "C00033"),
    "ACALD": ("acetaldehyde", "C00084"),
    "ETOH": ("ethanol", "C00469"),
    "ETOHex": ("Ethanol, extracellular", "C00469"),
    "ACex": ("Acetate, extracellular", "C00033"),
    "FORex": ("Formate, extracellular", "C00058"),
    "ATP": ("adenosine triphosphate", "C00002"),
    "ADP": ("adenosine diphosphate", "C00008"),
    "Pi": ("phosphate", "C00009"),
    "NAD": ("nicotinamide adenosine dinucleotide", "C00003"),
    "NADH": ("nicotinamide adenosine dinucleotide (reduced)", "C00004"),
    "H2O": ("water", "C00001"),
    "H": ("proton", "C00080"),
    "COA": ("coenzyme A", "C00010"),
    "ATPex": ("ATP, extracellular", "C00002"),
    "ADPex": ("ADP, extracellular", "C00008"),
    "Piex": ("Phosphate, extracellular", "C00009"),
    "H2Oex": ("water, extracellular", "C00001"),
    "Hex": ("Proton, extracellular", "C00080"),
}

PSEUDO_METS = ["H2O", "H2Oex", "H", "Hex", "ATP", "ATPex", "ADP", "ADPex",
               "Pi", "Piex"]


def met_df(ids):
    rows = [{"ID": i, "Name": MET_KEGG[i][0], "KEGG ID": MET_KEGG[i][1]}
            for i in ids]
    return pd.DataFrame(rows, columns=io.METABOLITE_COLS)


def rxn_df(rows):
    return pd.DataFrame(
        [{"ID": r[0], "Name": r[1], "Reaction stoichiometry": r[2],
          "Reversibility": r[3]} for r in rows],
        columns=io.REACTION_COLS)


# Shared upper-glycolysis backbone (HXK import -> pyruvate).
EMP_CORE = [
    ("HXK", "Hexokinase", "GLC + ATP <=> G6P + ADP", 0),
    ("PGI", "Glucose-6-phosphate isomerase", "G6P <=> F6P", 0),
    ("PFK", "Phosphofructokinase", "F6P + ATP <=> FBP + ADP", 0),
    ("ALDO", "Aldolase", "FBP <=> G3P + DHAP", 0),
    ("TPI", "Triose phosphate isomerase", "DHAP <=> G3P", 0),
    ("GAPDH", "Glyceraldehyde-3-phosphate dehydrogenase",
     "G3P + Pi + NAD <=> BPG + NADH + H", 0),
    ("PGK", "Phosphoglycerate kinase", "BPG + ADP <=> PG3 + ATP", 0),
    ("PGM", "Phosphoglycerate mutase", "PG3 <=> PG2", 0),
    ("ENO", "Enolase", "PG2 <=> PEP + H2O", 0),
    ("PYK", "Pyruvate kinase", "PEP + ADP <=> PYR + ATP", 0),
]

CURRENCY_EXCHANGES = [
    ("TATP", "ATP transport", "ATP <=> ATPex", 0),
    ("EXATP", "ATP exchange", "ATPex <=>", 0),
    ("TADP", "ADP transport", "ADP <=> ADPex", 1),
    ("EXADP", "ADP exchange", "ADPex <=>", 1),
    ("TransPi", "Phosphate transport", "Pi <=> Piex", 1),
    ("EXPi", "Phosphate exchange", "Piex <=>", 1),
    ("TH2O", "Water transport", "H2O <=> H2Oex", 1),
    ("EXH2O", "Water exchange", "H2Oex <=>", 1),
    ("TH", "Proton transport", "H <=> Hex", 1),
    ("EXH", "Proton exchange", "Hex <=>", 1),
]


def build_lactate():
    rxns = (
        EMP_CORE
        + [("LDH", "Lactate dehydrogenase", "PYR + NADH + H <=> LAC + NAD", 0)]
        + [("TGLC", "Glucose transport", "GLC <=> GLCex", 1),
           ("EXGLC", "Glucose exchange", "GLCex <=>", 1),
           ("TLAC", "Lactate transport", "LAC <=> LACex", 0),
           ("EXLAC", "Lactate exchange", "LACex <=>", 0)]
        + CURRENCY_EXCHANGES
    )
    mets = ["GLCex", "GLC", "G6P", "F6P", "FBP", "G3P", "DHAP", "BPG", "PG3",
            "PG2", "PEP", "PYR", "LAC", "LACex", "ATP", "ADP", "Pi", "NAD",
            "NADH", "H2O", "H", "ATPex", "ADPex", "Piex", "H2Oex", "Hex"]
    config = {
        "MODEL_NAME": "Example1_EMP_lactate",
        "SUBSTRATES": ["EXGLC"],
        "PRODUCTS": ["EXATP", "EXLAC"],
        "CARBON_PRODUCTS": ["EXLAC"],
        "ENERGY_PRODUCT": "EXATP",
        "REV_ALLOWED": ["EXADP", "EXPi", "EXH2O", "EXH"],
        "PSEUDO_METS": PSEUDO_METS,
    }
    return met_df(mets), rxn_df(rxns), config


def build_etoh():
    """The notebook's Example1 (ethanol/acetate/formate), PTS import."""
    rxns = [
        ("PTS", "Phosphotransferase system", "GLCex + PEP <=> G6P + PYR", 0),
        ("PGI", "Glucose-6-phosphate isomerase", "G6P <=> F6P", 0),
        ("PFK", "Phosphofructokinase", "F6P + ATP <=> FBP + ADP", 0),
        ("ALDO", "Aldolase", "FBP <=> G3P + DHAP", 0),
        ("TPI", "Triose phosphate isomerase", "DHAP <=> G3P", 0),
        ("GAPDH", "Glyceraldehyde-3-phosphate dehydrogenase",
         "G3P + Pi + NAD <=> BPG + NADH + H", 0),
        ("PGK", "Phosphoglycerate kinase", "BPG + ADP <=> PG3 + ATP", 0),
        ("PGM", "Phosphoglycerate mutase", "PG3 <=> PG2", 0),
        ("ENO", "Enolase", "PG2 <=> PEP + H2O", 0),
        ("PYK", "Pyruvate kinase", "PEP + ADP <=> PYR + ATP", 0),
        ("PFL", "Pyruvate formate lyase", "PYR + COA <=> ACCOA + FOR", 0),
        ("PTA", "Phosphotransacetylase", "ACCOA + Pi <=> ACTP + COA", 0),
        ("ACK", "Acetate kinase", "ACTP + ADP <=> AC + ATP", 0),
        ("ACALD", "Acetaldehyde dehydrogenase",
         "ACCOA + NADH + H <=> ACALD + COA + NAD", 0),
        ("ADH", "Alcohol dehydrogenase", "NADH + H + ACALD <=> ETOH + NAD", 0),
        ("EXGLC", "Glucose exchange", "GLCex <=>", 1),
        ("TETOH", "Ethanol transport", "ETOH <=> ETOHex", 0),
        ("EXETOH", "Ethanol exchange", "ETOHex <=>", 0),
        ("TAC", "Acetate transport", "AC <=> ACex", 0),
        ("EXAC", "Acetate exchange", "ACex <=>", 0),
        ("TFOR", "Formate transport", "FOR <=> FORex", 0),
        ("EXFOR", "Formate exchange", "FORex <=>", 0),
    ] + CURRENCY_EXCHANGES
    mets = ["GLCex", "GLC", "G6P", "F6P", "FBP", "G3P", "DHAP", "BPG", "PG3",
            "PG2", "PEP", "PYR", "FOR", "ACCOA", "ACTP", "AC", "ACALD", "ETOH",
            "ATP", "ADP", "Pi", "NAD", "NADH", "H2O", "H", "COA", "ETOHex",
            "ACex", "FORex", "ATPex", "ADPex", "Piex", "H2Oex", "Hex"]
    config = {
        "MODEL_NAME": "Example1_EMPglycolysis",
        "SUBSTRATES": ["EXGLC"],
        "PRODUCTS": ["EXATP", "EXETOH", "EXAC", "EXFOR"],
        "CARBON_PRODUCTS": ["EXETOH", "EXAC", "EXFOR"],
        "ENERGY_PRODUCT": "EXATP",
        "REV_ALLOWED": ["EXADP", "EXPi", "EXH2O", "EXH"],
        "PSEUDO_METS": PSEUDO_METS,
    }
    return met_df(mets), rxn_df(rxns), config


def add_pan_config():
    """Load the repo's pan_glycolysis.xlsx, attach a Config sheet, save here."""
    root = os.path.dirname(os.path.dirname(HERE))
    src = os.path.join(root, "pan_glycolysis.xlsx")
    df_met, df_rxn = io.load_excel(src)
    config = {
        "MODEL_NAME": "pan_glycolysis",
        "SUBSTRATES": ["EXGLC"],
        "PRODUCTS": ["EXATP", "EXETOH", "EXAC", "EXLAC", "EXFOR", "EXCO2"],
        "CARBON_PRODUCTS": ["EXETOH", "EXAC", "EXLAC", "EXFOR", "EXCO2"],
        "ENERGY_PRODUCT": "EXATP",
        "REV_ALLOWED": ["EXADP", "EXPi", "EXH2O", "EXH"],
        "PSEUDO_METS": PSEUDO_METS,
    }
    io.save_excel(df_met, df_rxn,
                  os.path.join(HERE, "pan_glycolysis.xlsx"), config)


def main():
    m, r, c = build_lactate()
    io.save_excel(m, r, os.path.join(HERE, "Example1_EMP_lactate.xlsx"), c)
    print("wrote Example1_EMP_lactate.xlsx")
    m, r, c = build_etoh()
    io.save_excel(m, r, os.path.join(HERE, "Example1_EMPglycolysis.xlsx"), c)
    print("wrote Example1_EMPglycolysis.xlsx")
    add_pan_config()
    print("wrote pan_glycolysis.xlsx")


if __name__ == "__main__":
    main()
