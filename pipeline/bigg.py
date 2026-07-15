"""KEGG-compound -> BiGG-id translation for human-readable working ids.

This is a *convenience* layer, never ground truth (CLAUDE.md keeps the KEGG id
load-bearing for balance + eQuilibrator).  It:

1. downloads BiGG's official ``bigg_models_metabolites.txt``/
   ``bigg_models_reactions.txt`` once and builds ``{KEGG id: [BiGG ids]}``
   maps, cached on disk like the KEGG cache;
2. picks a deterministic BiGG id per KEGG id (shortest, then alphabetical —
   the base form ``h2o`` beats the variant ``oh1``);
3. rewrites the working ids of a metabolite/reaction pair to those BiGG ids,
   **keeping the ``KEGG ID`` / ``KEGG Reaction ID`` columns intact** so
   thermodynamics and balance still resolve compounds by KEGG.

BiGG's own namespace files only carry a KEGG cross-reference for ids that
happen to appear in one of the ~100 genome-scale models BiGG hosts — most
KEGG reactions have no such annotation at all (only ~10% of BiGG's ~28k
reactions have a KEGG link). Reliability on what BiGG *does* cover is good
(~87% clean 1:1 on central-carbon metabolism), but coverage itself is the
bottleneck. As an optional, opt-in fallback (:func:`get_kegg_to_bigg_via_mnx`
/ :func:`get_kegg_to_bigg_reactions_via_mnx`), MetaNetX's ``chem_xref.tsv`` /
``reac_xref.tsv`` group many databases' local ids under a shared MNX id, so
joining the KEGG- and BiGG-prefixed rows that share an MNX id recovers many
KEGG<->BiGG links BiGG's own file never listed. This is a genuinely large
download (~700MB for compounds, ~80MB for reactions) — opt-in, not automatic,
and clearly labelled in the UI. Drop the files themselves at
``.cache/chem_xref.tsv`` / ``.cache/reac_xref.tsv`` (next to the pipeline
package; gitignored) to use a local copy and skip the download entirely.

Ids with no clean BiGG id (from either source) fall back to a slugified KEGG
name and are flagged in the report so the user can hand-edit.
"""

from __future__ import annotations

import os
import re
import urllib.request

import pandas as pd

from . import cache
from .idmap import apply_id_map, merge_duplicate_ids
from .kegg import get_compound_name

BIGG_URL = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"
_CACHE_KEY = "bigg:kegg2bigg"

BIGG_REACTIONS_URL = "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt"
_CACHE_KEY_RXN = "bigg:kegg2bigg_reactions"

# Optional fallback source: MetaNetX cross-references many databases' local
# ids under a shared MNX id, giving much broader (but indirect / two-hop)
# KEGG<->BiGG coverage than BiGG's own namespace files. Large downloads —
# opt-in via translate_to_bigg(..., use_mnx=True), not fetched by default.
MNX_CHEM_XREF_URL = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv"
MNX_REAC_XREF_URL = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv"
_CACHE_KEY_MNX_CHEM = "bigg:kegg2bigg_via_mnx"
_CACHE_KEY_MNX_REAC = "bigg:kegg2bigg_reactions_via_mnx"
_MNX_OBSOLETE = "secondary/obsolete/fantasy identifier"

# A manually-placed copy of the MetaNetX files (same directory the KEGG/BiGG
# on-disk cache lives in — see pipeline.cache — and equally gitignored) is
# used instead of downloading, when present. Override with env vars if you'd
# rather keep them elsewhere.
_CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), ".cache")
LOCAL_MNX_CHEM_PATH = os.environ.get(
    "MNX_CHEM_XREF_PATH", os.path.join(_CACHE_DIR, "chem_xref.tsv"))
LOCAL_MNX_REAC_PATH = os.environ.get(
    "MNX_REAC_XREF_PATH", os.path.join(_CACHE_DIR, "reac_xref.tsv"))

TIMEOUT = 90
MNX_TIMEOUT = 900  # the MetaNetX xref files are hundreds of MB

_MAP = None          # in-process cache of {kegg_compound_id: [bigg_id, ...]}
_MAP_RXN = None      # in-process cache of {kegg_reaction_id: [bigg_id, ...]}
_MAP_MNX = None      # in-process cache of the MetaNetX-derived compound map
_MAP_MNX_RXN = None  # in-process cache of the MetaNetX-derived reaction map


class BiggError(RuntimeError):
    """Raised when the BiGG namespace file can't be fetched."""


# --------------------------------------------------------------------------- #
# The KEGG -> BiGG map
# --------------------------------------------------------------------------- #
def get_kegg_to_bigg(force_refresh: bool = False):
    """Return ``{kegg_id: [bigg_id, ...]}``, building + caching it on first use.

    Raises :class:`BiggError` if the file must be downloaded but can't be
    reached (offline) and nothing is cached yet.
    """
    global _MAP
    if _MAP is not None and not force_refresh:
        return _MAP
    if not force_refresh:
        cached = cache.get(_CACHE_KEY)
        if cached:
            _MAP = cached
            return _MAP
    _MAP = _download_and_parse()
    cache.put(_CACHE_KEY, _MAP)
    return _MAP


def _download_and_parse():
    try:
        raw = urllib.request.urlopen(BIGG_URL, timeout=TIMEOUT).read().decode(
            "utf-8", "replace")
    except Exception as e:  # network, HTTP, timeout, ...
        raise BiggError(f"Could not download the BiGG namespace file: {e}") from e

    lines = raw.splitlines()
    if not lines:
        raise BiggError("BiGG namespace file was empty")
    header = lines[0].split("\t")
    try:
        li = header.index("database_links")
        bi = header.index("universal_bigg_id")
    except ValueError as e:
        raise BiggError(f"Unexpected BiGG file format: {e}") from e

    mapping: dict[str, list[str]] = {}
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) <= li:
            continue
        bigg = parts[bi].strip()
        if not bigg:
            continue
        for cid in re.findall(r"kegg\.compound/(C\d{5})", parts[li]):
            bag = mapping.setdefault(cid, [])
            if bigg not in bag:
                bag.append(bigg)
    return mapping


def get_kegg_to_bigg_reactions(force_refresh: bool = False):
    """Return ``{kegg_reaction_id: [bigg_reaction_id, ...]}``, cached like
    :func:`get_kegg_to_bigg`."""
    global _MAP_RXN
    if _MAP_RXN is not None and not force_refresh:
        return _MAP_RXN
    if not force_refresh:
        cached = cache.get(_CACHE_KEY_RXN)
        if cached:
            _MAP_RXN = cached
            return _MAP_RXN
    _MAP_RXN = _download_and_parse_reactions()
    cache.put(_CACHE_KEY_RXN, _MAP_RXN)
    return _MAP_RXN


def _download_and_parse_reactions():
    try:
        raw = urllib.request.urlopen(BIGG_REACTIONS_URL, timeout=TIMEOUT).read().decode(
            "utf-8", "replace")
    except Exception as e:  # network, HTTP, timeout, ...
        raise BiggError(f"Could not download the BiGG reaction namespace file: {e}") from e

    lines = raw.splitlines()
    if not lines:
        raise BiggError("BiGG reaction namespace file was empty")
    header = lines[0].split("\t")
    try:
        li = header.index("database_links")
        bi = header.index("bigg_id")
    except ValueError as e:
        raise BiggError(f"Unexpected BiGG reaction file format: {e}") from e

    mapping: dict[str, list[str]] = {}
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) <= li:
            continue
        bigg = parts[bi].strip()
        if not bigg:
            continue
        for rid in re.findall(r"kegg\.reaction/(R\d{5})", parts[li]):
            bag = mapping.setdefault(rid, [])
            if bigg not in bag:
                bag.append(bigg)
    return mapping


# --------------------------------------------------------------------------- #
# Optional fallback: MetaNetX cross-references (broader, large download)
# --------------------------------------------------------------------------- #
def get_kegg_to_bigg_via_mnx(force_refresh: bool = False):
    """Return ``{kegg_compound_id: [bigg_id, ...]}`` derived from MetaNetX.

    Joins the KEGG- and BiGG-prefixed rows of ``chem_xref.tsv`` that share an
    MNX id — much broader coverage than :func:`get_kegg_to_bigg`, at the cost
    of a large (~700MB) download (skipped if ``LOCAL_MNX_CHEM_PATH`` already
    exists on disk). Cached on disk after first use, like the other maps in
    this module.
    """
    global _MAP_MNX
    if _MAP_MNX is not None and not force_refresh:
        return _MAP_MNX
    if not force_refresh:
        cached = cache.get(_CACHE_KEY_MNX_CHEM)
        if cached:
            _MAP_MNX = cached
            return _MAP_MNX
    _MAP_MNX = _download_and_parse_mnx_xref(
        MNX_CHEM_XREF_URL, {"kegg.compound", "keggC"},
        {"bigg.metabolite", "biggM"}, re.compile(r"C\d{5}"),
        local_path=LOCAL_MNX_CHEM_PATH)
    cache.put(_CACHE_KEY_MNX_CHEM, _MAP_MNX)
    return _MAP_MNX


def get_kegg_to_bigg_reactions_via_mnx(force_refresh: bool = False):
    """Return ``{kegg_reaction_id: [bigg_id, ...]}`` derived from MetaNetX.

    Mirrors :func:`get_kegg_to_bigg_via_mnx` for reactions (``reac_xref.tsv``,
    ~80MB, skipped if ``LOCAL_MNX_REAC_PATH`` already exists on disk).
    """
    global _MAP_MNX_RXN
    if _MAP_MNX_RXN is not None and not force_refresh:
        return _MAP_MNX_RXN
    if not force_refresh:
        cached = cache.get(_CACHE_KEY_MNX_REAC)
        if cached:
            _MAP_MNX_RXN = cached
            return _MAP_MNX_RXN
    _MAP_MNX_RXN = _download_and_parse_mnx_xref(
        MNX_REAC_XREF_URL, {"kegg.reaction", "keggR"},
        {"bigg.reaction", "biggR"}, re.compile(r"R\d{5}"),
        local_path=LOCAL_MNX_REAC_PATH)
    cache.put(_CACHE_KEY_MNX_REAC, _MAP_MNX_RXN)
    return _MAP_MNX_RXN


def _download_and_parse_mnx_xref(url, kegg_prefixes, bigg_prefixes, id_pattern,
                                 local_path=None):
    """Parse a MetaNetX ``*_xref.tsv`` file into ``{kegg_id: [bigg_id, ...]}``.

    Each row is ``namespace:local_id <TAB> mnx_id <TAB> description``; every
    namespace's local id for the same underlying compound/reaction shares an
    ``mnx_id``, so this groups by that id and joins whichever rows carry a
    KEGG prefix against whichever carry a BiGG prefix. Deprecated aliases
    (``description == "secondary/obsolete/fantasy identifier"``) are skipped.

    If ``local_path`` exists on disk, it's read directly instead of hitting
    the network — drop a manually-downloaded copy there to skip the download.
    """
    if local_path and os.path.isfile(local_path):
        with open(local_path, "r", encoding="utf-8", errors="replace") as fh:
            raw = fh.read()
    else:
        try:
            raw = urllib.request.urlopen(url, timeout=MNX_TIMEOUT).read().decode(
                "utf-8", "replace")
        except Exception as e:  # network, HTTP, timeout, ...
            raise BiggError(f"Could not download the MetaNetX xref file: {e}") from e

    mnx_to_kegg: dict[str, set] = {}
    mnx_to_bigg: dict[str, set] = {}
    for ln in raw.splitlines():
        if not ln or ln[0] == "#":
            continue
        parts = ln.split("\t")
        if len(parts) < 3 or ":" not in parts[0]:
            continue
        source, mnx_id, description = parts[0], parts[1], parts[2]
        if description == _MNX_OBSOLETE:
            continue
        prefix, local_id = source.split(":", 1)
        if prefix in kegg_prefixes:
            if id_pattern.fullmatch(local_id):
                mnx_to_kegg.setdefault(mnx_id, set()).add(local_id)
        elif prefix in bigg_prefixes:
            mnx_to_bigg.setdefault(mnx_id, set()).add(local_id)

    mapping: dict[str, list[str]] = {}
    for mnx_id, kegg_ids in mnx_to_kegg.items():
        bigg_ids = mnx_to_bigg.get(mnx_id)
        if not bigg_ids:
            continue
        for kid in kegg_ids:
            mapping.setdefault(kid, sorted(bigg_ids))
    return mapping


def _merge_bigg_maps(primary, fallback):
    """Combine two ``{kegg_id: [bigg_id, ...]}`` maps; ``primary`` wins ties."""
    merged = dict(primary)
    for kid, biggs in fallback.items():
        if kid not in merged:
            merged[kid] = biggs
    return merged


def get_kegg_to_bigg_combined(force_refresh: bool = False):
    """:func:`get_kegg_to_bigg`, filled in with MetaNetX-derived entries for
    any KEGG compound BiGG's own file doesn't cover. Triggers the large
    MetaNetX download on first use — see :func:`get_kegg_to_bigg_via_mnx`."""
    return _merge_bigg_maps(get_kegg_to_bigg(force_refresh),
                            get_kegg_to_bigg_via_mnx(force_refresh))


def get_kegg_to_bigg_reactions_combined(force_refresh: bool = False):
    """Reaction counterpart of :func:`get_kegg_to_bigg_combined`."""
    return _merge_bigg_maps(get_kegg_to_bigg_reactions(force_refresh),
                            get_kegg_to_bigg_reactions_via_mnx(force_refresh))


def choose_bigg(kegg_id, mapping):
    """Deterministic pick: shortest BiGG id, ties broken alphabetically."""
    ids = mapping.get(kegg_id)
    if not ids:
        return None
    return sorted(ids, key=lambda s: (len(s), s))[0]


# --------------------------------------------------------------------------- #
# Slugify (fallback when there's no clean BiGG id)
# --------------------------------------------------------------------------- #
def slugify(name, fallback):
    """Turn a compound name into a compact, readable, id-safe token."""
    if not name or not str(name).strip():
        return str(fallback)
    s = str(name).strip()
    s = s.split(";")[0]                      # first synonym only
    s = re.sub(r"[^\w\s-]", "", s)           # drop punctuation (keep - and _)
    s = re.sub(r"\s+", "_", s.strip())       # spaces -> underscore
    s = s.strip("_-")
    return s or str(fallback)


# --------------------------------------------------------------------------- #
# Shared collision/merge resolution
# --------------------------------------------------------------------------- #
def _resolve_collision(old, new, kegg, kegg_valid, used, kegg_of_target):
    """If ``new`` is already claimed, either merge (same KEGG id) or suffix.

    Returns ``(new, status)`` where ``status`` is ``"merged"``, ``"collision"``,
    or ``None`` if there was no collision to resolve.
    """
    if new not in used or new == old:
        return new, None
    owner_kegg = kegg_of_target.get(new)
    if kegg_valid and owner_kegg == kegg:
        # Same KEGG id already claimed this target — e.g. a reaction/compound
        # added after an earlier BiGG conversion, still referencing the raw
        # KEGG id. Merge into the existing row instead of minting a "_2"
        # near-duplicate; apply_id_map / merge_duplicate_ids collapse the
        # now-identical rows.
        return new, "merged"
    # Genuine clash between two different KEGG ids: suffix it.
    base, n = new, 2
    while new in used:
        new = f"{base}_{n}"
        n += 1
    return new, "collision"


# --------------------------------------------------------------------------- #
# Translate a metabolite / reaction pair
# --------------------------------------------------------------------------- #
def translate_to_bigg(df_metabolites, df_reactions, mapping=None, use_mnx=False):
    """Rewrite working ids to BiGG ids (KEGG id kept in the ``KEGG ID`` column).

    Returns ``(df_met, df_rxn, report)``.  ``report`` is a list of dicts, one per
    metabolite: ``{old, new, kegg, status}`` with ``status`` in
    ``{"bigg", "bigg-ambiguous", "name-fallback", "kegg-fallback", "collision"}``.
    Reaction stoichiometry strings and ``EX{id}`` reaction ids are rewritten to
    match; the reaction that carries no mapped id is left untouched.

    ``use_mnx``: also consult the MetaNetX-derived fallback map for KEGG ids
    BiGG's own file doesn't cover (see :func:`get_kegg_to_bigg_combined`) —
    triggers a large one-time download on first use. Ignored if ``mapping``
    is given explicitly.
    """
    if mapping is None:
        mapping = get_kegg_to_bigg_combined() if use_mnx else get_kegg_to_bigg()
    df_met = df_metabolites.copy()
    df_rxn = df_reactions.copy()

    id_map, report = {}, []
    used = set()
    kegg_of_target: dict[str, str | None] = {}  # target id -> the KEGG id that claimed it
    for _, row in df_met.iterrows():
        old = str(row.get("ID", "")).strip()
        kegg = str(row.get("KEGG ID", "")).strip()
        if not old:
            continue
        kegg_valid = bool(re.fullmatch(r"C\d{5}", kegg))
        biggs = mapping.get(kegg, []) if kegg_valid else []
        if biggs:
            new = choose_bigg(kegg, mapping)
            status = "bigg" if len(biggs) == 1 else "bigg-ambiguous"
        elif kegg_valid:
            new = slugify(get_compound_name(kegg), kegg)
            status = "name-fallback"
        else:
            new = old
            status = "kegg-fallback"

        new, coll_status = _resolve_collision(old, new, kegg, kegg_valid, used,
                                              kegg_of_target)
        if coll_status:
            status = coll_status
        used.add(new)
        kegg_of_target.setdefault(new, kegg if kegg_valid else None)
        if new != old:
            id_map[old] = new
        report.append({"old": old, "new": new, "kegg": kegg, "status": status})

    df_met, df_rxn = apply_id_map(df_met, df_rxn, id_map)
    return df_met, df_rxn, report


def translate_reactions_to_bigg(df_reactions, mapping=None, use_mnx=False):
    """Rewrite reaction working ids to BiGG reaction ids.

    Mirrors :func:`translate_to_bigg` but for reactions: a reaction id has no
    cross-references elsewhere in the tables (only metabolite ids appear
    inside stoichiometries), so this only ever rewrites the ``ID`` column —
    no equation rewriting needed. The originating KEGG reaction id is kept
    (backfilled if blank) in the ``KEGG Reaction ID`` column, the same
    pattern :func:`translate_to_bigg` uses for metabolites' ``KEGG ID``.

    ``use_mnx``: see :func:`translate_to_bigg` — same large-download tradeoff,
    the reaction file (~80MB) is smaller than the compound one (~700MB).

    Returns ``(df_rxn, report)`` with the same ``status`` vocabulary as
    :func:`translate_to_bigg`.
    """
    if mapping is None:
        mapping = (get_kegg_to_bigg_reactions_combined() if use_mnx
                  else get_kegg_to_bigg_reactions())
    df_rxn = df_reactions.copy()
    if "KEGG Reaction ID" not in df_rxn.columns:
        df_rxn["KEGG Reaction ID"] = pd.Series(dtype=object)

    id_map, report = {}, []
    used = set()
    kegg_of_target: dict[str, str | None] = {}
    row_kegg = {}  # row index -> resolved KEGG reaction id, for backfill

    for idx, row in df_rxn.iterrows():
        old = str(row.get("ID", "")).strip()
        if not old:
            continue
        kegg = str(row.get("KEGG Reaction ID", "") or "").strip()
        if not kegg and re.fullmatch(r"R\d{5}", old):
            kegg = old  # not yet translated: the working id is still the KEGG id
        row_kegg[idx] = kegg

        kegg_valid = bool(re.fullmatch(r"R\d{5}", kegg))
        biggs = mapping.get(kegg, []) if kegg_valid else []
        if biggs:
            new = choose_bigg(kegg, mapping)
            status = "bigg" if len(biggs) == 1 else "bigg-ambiguous"
        elif kegg_valid:
            new = slugify(row.get("Name") or kegg, kegg)
            status = "name-fallback"
        else:
            new = old
            status = "kegg-fallback"

        new, coll_status = _resolve_collision(old, new, kegg, kegg_valid, used,
                                              kegg_of_target)
        if coll_status:
            status = coll_status
        used.add(new)
        kegg_of_target.setdefault(new, kegg if kegg_valid else None)
        if new != old:
            id_map[old] = new
        report.append({"old": old, "new": new, "kegg": kegg, "status": status})

    for idx, kegg in row_kegg.items():
        if kegg:
            df_rxn.loc[idx, "KEGG Reaction ID"] = kegg
    df_rxn["ID"] = df_rxn["ID"].map(lambda x: id_map.get(str(x).strip(), x))
    df_rxn = merge_duplicate_ids(df_rxn, id_col="ID")
    return df_rxn, report
