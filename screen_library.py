#!/usr/bin/env python3
"""
screen_library.py
=================

Step 2b of the Alpha-Ligand workflow: the REAL-LIBRARY alternative to
`ligand_builder.py`.

WHY THIS EXISTS
---------------
`ligand_builder.py` assembles molecules from a small hand-written kit of
ring cores and side groups (~640 molecules). That is enough to see the
pipeline work, but the molecules are small fragments, so their docking
scores land well behind a real drug.

This script does the thing a real discovery project does instead: it takes
a big collection of REAL, already-synthesised compounds and filters it down
to the handful that best match the pocket's shopping list, then docks those.

WHAT IT CANNOT DO
-----------------
It cannot screen "every molecule that works". Chemical space is far too
large to enumerate (estimates run to 10^60 drug-like molecules). What it
can do is screen a real, finite catalogue of compounds that chemists have
actually made. That is the honest version of "everything".

THE LIBRARY
-----------
Default source: ChEMBL (https://www.ebi.ac.uk/chembl/), the EBI's open
database of bioactive molecules extracted from the medicinal-chemistry
literature. The file used is `chembl_<N>_chemreps.txt.gz`, a plain
tab-separated list of ChEMBL ID + SMILES for every compound in the
release (~2.9 million molecules in ChEMBL 37, ~293 MB compressed).

These are molecules that have been synthesised and assayed in published
work. They are NOT guaranteed to be purchasable, and they are NOT a
"things that bind your protein" list - they bind all sorts of targets, or
nothing much. Being in ChEMBL only means someone made it and measured it.

ChEMBL data is released under CC BY-SA 3.0.

HOW THE FILTERING WORKS (and why it is crude)
---------------------------------------------
Stages, in order, each one cheaper than the next:

  1. Text prefilter   - reject on SMILES string length and on elements we
                        cannot dock sensibly (metals, etc). No chemistry.
  2. RDKit parse      - turn the survivors into real molecules, keep the
                        largest fragment (strips salt counter-ions).
  3. Property windows - molecular weight / cLogP / H-bond donors /
                        acceptors / rotatable bonds, taken from the pocket
                        profile itself.
  4. Structural alerts- optional PAINS + Brenk filters (RDKit FilterCatalog),
                        which flag notorious assay-interfering or unstable
                        groups.
  5. Checklist score  - the SAME 0-100 pharmacophore match used by
                        ligand_builder.py (shared code in
                        pharmacophore_score.py).
  6. Scaffold spread  - cap how many hits may share one Murcko scaffold, so
                        the top of the list is not ten near-twins.
  7. 3D conformers    - for the final shortlist only.

Stages 3-5 are a CRUDE PREFILTER. They ask "is this molecule roughly the
right size and does it carry roughly the right kinds of group", with no
protein, no shape matching and no energy. Their only job is to cut
millions of molecules down to a number that can actually be docked.
The docking step is what produces a score worth reading, and even that is
a weak predictor of real binding.

OUTPUT
------
Writes `<CODE>_ligands.sdf` in exactly the same format ligand_builder.py
produces, so `dock_ligands.py` and `make_report.py` work unchanged, plus a
`<CODE>_ligand_report.txt` with the stage-by-stage counts.
"""

import argparse
import gzip
import io
import os
import random
import sys
import time
import urllib.request

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from pharmacophore_score import (
    load_pharmacophore,
    summarise_demands,
    measure,
    score_candidate,
    make_3d,
)

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# SECTION 1: Library sources
# ---------------------------------------------------------------------------

CHEMBL_URL = ("https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
              "chembl_37_chemreps.txt.gz")
CHEMBL_NAME = "ChEMBL 37 (chemreps)"

# Bundled offline fallback: a small diverse slice of the same ChEMBL release,
# shipped in the repo so the script still does something useful with no
# internet (or if EBI is down).
FALLBACK_FILE = os.path.join(HERE, "data", "chembl_fallback_subset.smi")

# Elements we are willing to dock. Anything else (metals, boron clusters,
# silicon, selenium...) is dropped: Vina has no sensible parameters for them.
ALLOWED_ELEMENT_CHARS = set("BCNOPSFIHbcnopsfi")   # within-SMILES letters
BANNED_TOKENS = ("[Pt", "[Fe", "[Zn", "[Cu", "[Mn", "[Co", "[Ni", "[Ru",
                 "[Ag", "[Au", "[Hg", "[As", "[Se", "[Si", "[Sn", "[Li",
                 "[Na", "[K+", "[Ca", "[Mg", "[Al", "[Ti", "[Gd", "[Tc",
                 "[B-", "[B@", "[BH", "B(", "[Te", "[Sb", "[Bi", "[Zr",
                 "[Cr", "[Mo", "[W]", "[V]", "[Re", "[Rh", "[Pd", "[Ir",
                 "[La", "[Y]", "[Sr", "[Ba", "[Cs", "[Rb", "[Be")


def cheap_text_ok(smiles, max_len):
    """A no-chemistry string test, run on millions of lines, so it has to be
    fast. Rejects obviously-too-big molecules and unsupported elements."""
    if len(smiles) > max_len or len(smiles) < 6:
        return False
    for token in BANNED_TOKENS:
        if token in smiles:
            return False
    return True


def stream_library(url, cache_path, max_rows, max_len, verbose=True):
    """Yield (id, smiles) from a gzipped ChEMBL chemreps file.

    Reads straight from the network and decompresses on the fly, so when
    `max_rows` is reached we stop and never download the rest of the file.
    If `cache_path` already exists it is read from disk instead.
    """
    if cache_path and os.path.isfile(cache_path):
        if verbose:
            print("   using cached library file: %s (%.0f MB)"
                  % (cache_path, os.path.getsize(cache_path) / 1e6))
        handle = gzip.open(cache_path, "rt")
        source = "cache"
    else:
        if verbose:
            print("   streaming %s" % url)
            print("   (download stops as soon as enough rows are read)")
        request = urllib.request.Request(url, headers={"User-Agent": "alpha-ligand"})
        raw = urllib.request.urlopen(request, timeout=60)
        handle = io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8",
                                  errors="replace")
        source = "network"
    with handle:
        first = True
        read = 0
        for line in handle:
            if first:
                first = False
                if line.startswith("chembl_id"):
                    continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            read += 1
            if max_rows and read > max_rows:
                break
            yield parts[0], parts[1].strip()
    if verbose:
        print("   read %s rows from %s" % (f"{read:,}", source))


def read_plain_smiles(path):
    """Read a two-column '<smiles> <id>' or '<id>\\t<smiles>' file."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            bits = line.replace("\t", " ").split()
            if len(bits) == 1:
                yield "mol%d" % id(line), bits[0]
            elif bits[0].startswith("CHEMBL") or bits[0].startswith("ZINC"):
                yield bits[0], bits[1]
            else:
                yield bits[1], bits[0]


# ---------------------------------------------------------------------------
# SECTION 2: Filters
# ---------------------------------------------------------------------------

def largest_fragment(mol):
    """Salts arrive as 'drug.HCl'. Keep the biggest piece."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) == 1:
        return mol
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def windows_from_pharmacophore(data, slack):
    """Property windows suggested by pocket_profiler.py, optionally widened.

    slack = 1.0 means 'use them exactly as given'. Larger values widen each
    window around its midpoint, which lets more of the library through.
    """
    w = dict(data.get("suggested_ligand_properties", {}))
    if slack == 1.0:
        return w
    out = dict(w)
    for key in ("molecular_weight_Da", "clogp", "hbond_donors", "hbond_acceptors"):
        if key in w and isinstance(w[key], (list, tuple)) and len(w[key]) == 2:
            lo, hi = float(w[key][0]), float(w[key][1])
            mid = 0.5 * (lo + hi)
            half = 0.5 * (hi - lo) * slack
            out[key] = [mid - half, mid + half]
    if "rotatable_bonds_max" in w:
        out["rotatable_bonds_max"] = w["rotatable_bonds_max"] * slack
    return out


def inside_windows(props, windows):
    mw_lo, mw_hi = windows.get("molecular_weight_Da", [0, 1e9])
    if not (mw_lo <= props["mw"] <= mw_hi):
        return False
    lp_lo, lp_hi = windows.get("clogp", [-99, 99])
    if not (lp_lo <= props["clogp"] <= lp_hi):
        return False
    hbd_lo, hbd_hi = windows.get("hbond_donors", [0, 99])
    if not (hbd_lo <= props["hbd"] <= hbd_hi):
        return False
    hba_lo, hba_hi = windows.get("hbond_acceptors", [0, 99])
    if not (hba_lo <= props["hba"] <= hba_hi):
        return False
    if props["rotb"] > windows.get("rotatable_bonds_max", 99):
        return False
    return True


def make_alert_catalog(use_pains, use_brenk):
    if not (use_pains or use_brenk):
        return None
    from rdkit.Chem import FilterCatalog
    params = FilterCatalog.FilterCatalogParams()
    if use_pains:
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    if use_brenk:
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)


# ---------------------------------------------------------------------------
# SECTION 3: Main
# ---------------------------------------------------------------------------

def human(n):
    return f"{n:,}"


def main():
    ap = argparse.ArgumentParser(
        description="Screen a real compound library (ChEMBL by default) "
                    "against a pocket_profiler.py pharmacophore, and write "
                    "the best-matching molecules as 3D structures for docking.")
    ap.add_argument("pharmacophore_json")
    ap.add_argument("--outdir", default="output")
    ap.add_argument("--library", default="chembl",
                    choices=["chembl", "fallback", "file"],
                    help="chembl = stream ChEMBL from EBI (default); "
                         "fallback = the small bundled subset; "
                         "file = your own SMILES file via --library-file")
    ap.add_argument("--library-file", dest="library_file", default=None,
                    help="SMILES file to screen when --library file")
    ap.add_argument("--cache", default=None,
                    help="Path to a downloaded chemreps .gz to reuse "
                         "(default: <outdir>/downloads/chembl_chemreps.txt.gz "
                         "if it exists)")
    ap.add_argument("--max-scan", dest="max_scan", type=int, default=250000,
                    help="How many library rows to look at. 0 = all of them "
                         "(ChEMBL 37 has ~2.9 million; that takes much "
                         "longer). Default: 250000")
    ap.add_argument("--max-smiles-len", dest="max_smiles_len", type=int, default=120,
                    help="Text prefilter: skip SMILES longer than this "
                         "(default: 120 characters)")
    ap.add_argument("--shortlist", type=int, default=500,
                    help="How many molecules survive the chemistry prefilter "
                         "and get kept in the report (default: 500)")
    ap.add_argument("--molecules", type=int, default=15,
                    help="How many top molecules get 3D structures and go on "
                         "to docking (default: 15)")
    ap.add_argument("--confs", type=int, default=5,
                    help="Conformers per molecule (default: 5)")
    ap.add_argument("--top", type=int, default=12,
                    help="How many pharmacophore features form the shopping "
                         "list (default: 12)")
    ap.add_argument("--window-slack", dest="window_slack", type=float, default=1.0,
                    help="Widen the pocket's property windows by this factor "
                         "(1.0 = as given; 1.5 = 50%% wider)")
    ap.add_argument("--auto-relax", dest="auto_relax", action="store_true",
                    default=True,
                    help="If too few molecules pass the property windows, "
                         "widen them once and rescan (default: on)")
    ap.add_argument("--no-auto-relax", dest="auto_relax", action="store_false")
    ap.add_argument("--max-per-scaffold", dest="max_per_scaffold", type=int, default=2,
                    help="Cap on molecules sharing one Murcko scaffold in the "
                         "final shortlist (default: 2)")
    ap.add_argument("--pains", dest="pains", action="store_true", default=True,
                    help="Drop molecules matching PAINS/Brenk structural "
                         "alerts (default: on)")
    ap.add_argument("--no-pains", dest="pains", action="store_false")
    ap.add_argument("--seed", type=int, default=0xF00D)
    ap.add_argument("--progress-every", dest="progress_every", type=int,
                    default=25000, help="Print a progress line every N rows")
    args = ap.parse_args()

    t0 = time.time()
    print("=" * 72)
    print("LIBRARY SCREEN (real compound collection)")
    print("=" * 72)

    data = load_pharmacophore(args.pharmacophore_json)
    code = data.get("structure", "ligands")
    features = data.get("pharmacophore_features", [])
    demands = summarise_demands(features, args.top)
    if not demands:
        sys.exit("ERROR: no usable chemistry classes in the top %d features."
                 % args.top)
    print("Pocket  : %s   (%d pharmacophore features, %d shell residues)"
          % (code, len(features), len(data.get("pocket_shell_residues", []))))
    print("Chemistry the pocket asks for:")
    for cls, demand in sorted(demands.items(), key=lambda kv: -kv[1]):
        print("   %-16s weight %6.1f" % (cls, demand))

    windows = windows_from_pharmacophore(data, args.window_slack)
    print("Property windows used (from the pocket profile%s):"
          % ("" if args.window_slack == 1.0 else
             ", widened x%.2f" % args.window_slack))
    print("   MW %s Da | cLogP %s | HBD %s | HBA %s | RotB <= %s"
          % (windows.get("molecular_weight_Da"), windows.get("clogp"),
             windows.get("hbond_donors"), windows.get("hbond_acceptors"),
             windows.get("rotatable_bonds_max")))

    # --- choose the library --------------------------------------------
    os.makedirs(args.outdir, exist_ok=True)
    downloads = os.path.join(args.outdir, "downloads")
    os.makedirs(downloads, exist_ok=True)
    cache = args.cache
    if cache is None:
        for guess in ("chembl_37_chemreps.txt.gz", "chembl_chemreps.txt.gz"):
            path = os.path.join(downloads, guess)
            if os.path.isfile(path):
                cache = path
                break

    print("")
    print("-" * 72)
    print("LIBRARY")
    print("-" * 72)

    def chembl_rows():
        return stream_library(CHEMBL_URL, cache, args.max_scan,
                              args.max_smiles_len)

    if args.library == "file":
        if not args.library_file or not os.path.isfile(args.library_file):
            sys.exit("ERROR: --library file needs an existing --library-file")
        library_name = "custom file: %s" % args.library_file
        rows_factory = lambda: read_plain_smiles(args.library_file)
    elif args.library == "fallback":
        if not os.path.isfile(FALLBACK_FILE):
            sys.exit("ERROR: bundled fallback subset missing: %s" % FALLBACK_FILE)
        library_name = "bundled ChEMBL fallback subset"
        rows_factory = lambda: read_plain_smiles(FALLBACK_FILE)
    else:
        library_name = CHEMBL_NAME
        rows_factory = chembl_rows

    print("Source: %s" % library_name)
    if args.library == "chembl":
        print("Rows to scan: %s%s"
              % ("ALL (~2.9 million)" if not args.max_scan
                 else human(args.max_scan),
                 "" if args.max_scan else "  - this will take a while"))

    # --- the scan -------------------------------------------------------
    def scan(active_windows):
        """One pass over the library. Returns (counts, scored_list)."""
        counts = {"rows": 0, "text_ok": 0, "parsed": 0, "in_windows": 0,
                  "alerts_dropped": 0, "scored": 0}
        catalog = make_alert_catalog(args.pains, args.pains)
        kept = []
        last = time.time()
        rows = rows_factory()
        try:
            for ident, smiles in rows:
                counts["rows"] += 1
                if counts["rows"] % args.progress_every == 0:
                    now = time.time()
                    print("   ...%s rows | %s passed properties | %s scored "
                          "| %.0f s elapsed"
                          % (human(counts["rows"]), human(counts["in_windows"]),
                             human(counts["scored"]), now - t0))
                    last = now
                if not cheap_text_ok(smiles, args.max_smiles_len):
                    continue
                counts["text_ok"] += 1
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                if "." in smiles:
                    mol = largest_fragment(mol)
                    try:
                        Chem.SanitizeMol(mol)
                    except Exception:
                        continue
                counts["parsed"] += 1
                try:
                    props = measure(mol)
                except Exception:
                    continue
                if not inside_windows(props, active_windows):
                    continue
                counts["in_windows"] += 1
                if catalog is not None and catalog.HasMatch(mol):
                    counts["alerts_dropped"] += 1
                    continue
                try:
                    info = score_candidate(mol, demands, active_windows)
                except Exception:
                    continue
                counts["scored"] += 1
                kept.append({"name": ident, "smiles": Chem.MolToSmiles(mol),
                             "mol": mol, "props": props, "score_info": info,
                             "origin": library_name})
        except KeyboardInterrupt:
            print("   interrupted by user - working with what we have")
        except Exception as err:
            print("   library read stopped early (%s)" % err)
        return counts, kept

    counts, kept = scan(windows)

    used_windows = windows
    relaxed = False
    if args.auto_relax and len(kept) < max(50, args.molecules * 5) \
            and args.library != "fallback":
        print("")
        print("Only %d molecules passed the property windows - too few to "
              "choose from." % len(kept))
        print("Widening the windows x1.6 and scanning again "
              "(this is a deliberate loosening, not a better model).")
        used_windows = windows_from_pharmacophore(data, args.window_slack * 1.6)
        relaxed = True
        counts, kept = scan(used_windows)

    if not kept:
        sys.exit("ERROR: nothing survived the filters. Try --window-slack 2 "
                 "or a larger --max-scan.")

    kept.sort(key=lambda c: -c["score_info"]["score"])

    # --- scaffold spread + shortlist ------------------------------------
    shortlist = []
    scaffold_use = {}
    for cand in kept:
        if len(shortlist) >= args.shortlist:
            break
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=cand["mol"])
        except Exception:
            scaf = cand["smiles"]
        if scaffold_use.get(scaf, 0) >= args.max_per_scaffold:
            continue
        scaffold_use[scaf] = scaffold_use.get(scaf, 0) + 1
        cand["scaffold"] = scaf
        shortlist.append(cand)

    print("")
    print("-" * 72)
    print("STAGED FILTER COUNTS  (each stage is cheaper than the next)")
    print("-" * 72)
    print("   library rows read             : %s" % human(counts["rows"]))
    print("   passed text prefilter         : %s" % human(counts["text_ok"]))
    print("   parsed by RDKit               : %s" % human(counts["parsed"]))
    print("   inside pocket property windows: %s%s"
          % (human(counts["in_windows"]), "  (windows widened)" if relaxed else ""))
    print("   dropped by PAINS/Brenk alerts : %s"
          % (human(counts["alerts_dropped"]) if args.pains else "off"))
    print("   chemistry-scored              : %s" % human(counts["scored"]))
    print("   shortlist after scaffold cap  : %s" % human(len(shortlist)))
    print("   going to 3D + docking         : %s" % human(min(args.molecules,
                                                              len(shortlist))))
    print("")
    print("NOTE: everything above the docking step is a CRUDE PREFILTER.")
    print("It checks size and the presence of chemical groups. There is no")
    print("protein in it, no shape fitting and no energy. Its only job is to")
    print("cut a huge catalogue down to something dockable.")

    # --- 3D for the winners ---------------------------------------------
    sdf_path = os.path.join(args.outdir, "%s_ligands.sdf" % code)
    txt_path = os.path.join(args.outdir, "%s_ligand_report.txt" % code)
    print("")
    print("Generating 3D conformers for the top %d molecules..." % args.molecules)
    writer = Chem.SDWriter(sdf_path)
    built = []
    total_confs = 0
    for cand in shortlist:
        if len(built) >= args.molecules:
            break
        mol_h, conf_ids, energies, field = make_3d(cand["mol"], args.confs,
                                                   seed=args.seed)
        if mol_h is None:
            print("   %-28s SKIPPED (%s)" % (cand["name"], field))
            continue
        order = sorted(range(len(conf_ids)), key=lambda i: energies[i])
        sc = cand["score_info"]
        for position, idx in enumerate(order, start=1):
            conf_id = conf_ids[idx]
            mol_h.SetProp("_Name", "%s_conf%d" % (cand["name"], position))
            mol_h.SetProp("molecule_name", cand["name"])
            mol_h.SetProp("smiles", cand["smiles"])
            mol_h.SetProp("source", "%s (real compound library)" % library_name)
            mol_h.SetProp("checklist_score_0_100", "%.1f" % sc["score"])
            mol_h.SetProp("chemistry_match_fraction", "%.3f" % sc["chemistry_fraction"])
            mol_h.SetProp("property_window_fraction", "%.3f" % sc["property_fraction"])
            mol_h.SetProp("forcefield", field)
            mol_h.SetProp("forcefield_energy_kcal_per_mol", "%.3f" % energies[idx])
            mol_h.SetProp("conformer_rank_by_strain", str(position))
            mol_h.SetProp("WARNING", "checklist_score is chemical feature "
                                     "coverage only - NOT a binding affinity "
                                     "and NOT a docking score")
            writer.write(mol_h, confId=conf_id)
            total_confs += 1
        cand["field"] = field
        built.append(cand)
        print("   %-28s score %5.1f  MW %6.1f  %2d confs  %s"
              % (cand["name"][:28], sc["score"], cand["props"]["mw"],
                 len(conf_ids), field))
    writer.close()

    if not built:
        sys.exit("ERROR: no molecule could be embedded in 3D.")

    # --- report ----------------------------------------------------------
    elapsed = time.time() - t0
    with open(txt_path, "w") as fh:
        w = fh.write
        w("=" * 72 + "\n")
        w("LIBRARY SCREEN REPORT - %s\n" % code)
        w("=" * 72 + "\n\n")
        w("Library                 : %s\n" % library_name)
        if args.library == "chembl":
            w("Library source URL      : %s\n" % CHEMBL_URL)
            w("What it contains        : real, published, synthesised compounds\n")
            w("                          with measured bioactivity data. Being in\n")
            w("                          ChEMBL does NOT mean a compound binds your\n")
            w("                          protein, and does NOT mean it is for sale.\n")
        w("Rows read               : %s\n" % human(counts["rows"]))
        w("Pharmacophore file      : %s\n" % args.pharmacophore_json)
        w("Runtime (screen only)   : %.1f s\n\n" % elapsed)
        w("-" * 72 + "\nSTAGED FILTER COUNTS\n" + "-" * 72 + "\n")
        w("   library rows read             : %s\n" % human(counts["rows"]))
        w("   passed text prefilter         : %s\n" % human(counts["text_ok"]))
        w("   parsed by RDKit               : %s\n" % human(counts["parsed"]))
        w("   inside pocket property windows: %s%s\n"
          % (human(counts["in_windows"]), " (widened)" if relaxed else ""))
        w("   dropped by PAINS/Brenk alerts : %s\n"
          % (human(counts["alerts_dropped"]) if args.pains else "off"))
        w("   chemistry-scored              : %s\n" % human(counts["scored"]))
        w("   shortlist (scaffold-capped)   : %s\n" % human(len(shortlist)))
        w("   built in 3D                   : %s\n\n" % human(len(built)))
        w("Property windows used:\n")
        for key, value in used_windows.items():
            w("   %-24s %s\n" % (key, value))
        if relaxed:
            w("\n(The windows were widened automatically because too few\n"
              "library molecules fitted the original ones.)\n")
        w("\n" + "-" * 72 + "\nWHAT THE POCKET ASKED FOR\n" + "-" * 72 + "\n")
        for cls, demand in sorted(demands.items(), key=lambda kv: -kv[1]):
            w("   %-16s demand weight %6.1f\n" % (cls, demand))
        w("\n" + "-" * 72 + "\nMOLECULES SENT TO DOCKING\n" + "-" * 72 + "\n")
        for rank, cand in enumerate(built, start=1):
            sc = cand["score_info"]
            p = cand["props"]
            w("\n%2d. %s   checklist %.1f/100\n" % (rank, cand["name"], sc["score"]))
            w("    SMILES   : %s\n" % cand["smiles"])
            w("    MW %.1f  cLogP %.2f  HBD %d  HBA %d  RotB %d  rings %d\n"
              % (p["mw"], p["clogp"], p["hbd"], p["hba"], p["rotb"], p["rings"]))
            w("    features : %s\n"
              % ", ".join("%s x%d" % (k, v)
                          for k, v in sorted(sc["feature_counts"].items())
                          if v))
            if cand.get("scaffold"):
                w("    scaffold : %s\n" % cand["scaffold"])
        w("\n" + "-" * 72 + "\nTHE HONEST BIT\n" + "-" * 72 + "\n")
        w("The checklist score is feature coverage plus a size check. It is\n"
          "NOT an affinity. The library screen cannot look at shape, and a\n"
          "molecule with perfect chemistry can still be unable to enter the\n"
          "pocket. Only the docking step (dock_ligands.py) puts a molecule\n"
          "inside the protein, and even a good docking score is a weak\n"
          "predictor of real binding. Treat the output as a shortlist worth\n"
          "testing, not a prediction.\n")
        if args.max_scan:
            w("\nThis run looked at the FIRST %s rows of the library, not a\n"
              "random sample of it and not all of it. ChEMBL rows are roughly\n"
              "in order of when the compound entered the database, so a\n"
              "capped scan is a historical slice, not a representative one.\n"
              "Use --max-scan 0 to screen the whole collection.\n"
              % human(args.max_scan))

    print("")
    print("-" * 72)
    print("DONE in %.1f s." % elapsed)
    print("  molecules written  : %d" % len(built))
    print("  conformers written : %d" % total_confs)
    print("  SDF (3D molecules) : %s" % sdf_path)
    print("  report (plain text): %s" % txt_path)
    print("-" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
