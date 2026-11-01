#!/usr/bin/env python3
"""
ligand_builder.py
=================

Step 2 of the Alpha-Ligand workflow.

WHAT IT DOES (in plain English)
-------------------------------
`pocket_profiler.py` looks at a protein pocket and writes a "pharmacophore"
JSON file: a shopping list of the chemistry a good binder should have
(a positive group here, a ring there, a greasy bit over there), plus
sensible size limits for the molecule (molecular weight, logP, etc.).

This script takes that shopping list and produces REAL 3D MOLECULES:

1. It reads the JSON file.
2. It builds a set of candidate molecules out of a small, hand-written
   library of chemical scaffolds (ring cores) and substituent groups
   (the bits you hang off them), plus some whole known fragments.
3. For every candidate it uses RDKit to ask factual questions:
   does it have a positive group? a carboxylate? an aromatic ring?
   how greasy is it? how much does it weigh?
4. It scores each candidate by HOW WELL THOSE FACTS MATCH THE SHOPPING LIST.
5. It keeps the best ones, generates 3D conformers (actual x/y/z atom
   coordinates) with RDKit's ETKDG embedder, energy-minimises them with
   MMFF (or UFF if MMFF has no parameters), and writes them to an SDF file
   you can open in PyMOL, Chimera, Avogadro or feed to a docking program.
6. It writes a plain-English report.

WHAT THE SCORE IS *NOT*
------------------------
The score is NOT a binding affinity. It is NOT a docking score. It is NOT
an energy. Nothing in this script puts the molecule inside the protein or
calculates any interaction with the protein. The score is only a
"does this molecule's chemistry tick the boxes on the list" percentage.
Two molecules with the same score can bind a million-fold differently.
To find out whether anything actually binds you need docking, then
simulation, then a real experiment.

Only RDKit is used here (no Biopython needed).
"""

import argparse
import json
import math
import os
import sys

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

# RDKit prints a lot of chatter; quieten it so our own output is readable.
RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# SECTION 1: Shared chemistry code
#
# The SMARTS patterns, the pharmacophore reader, the property measurements,
# the 0-100 checklist score and the 3D conformer generator all live in
# pharmacophore_score.py, so that this script and screen_library.py (which
# screens a real compound library) score molecules in EXACTLY the same way.
# ---------------------------------------------------------------------------

from pharmacophore_score import (        # noqa: E402
    SMARTS_PATTERNS,
    COMPILED_PATTERNS,
    FEATURE_ALIASES,
    load_pharmacophore,
    summarise_demands,
    count_feature,
    measure,
    property_fit,
    score_candidate,
    make_3d,
)



# ---------------------------------------------------------------------------
# SECTION 2: The building-block library
#
# This is a hand-written, deliberately SMALL library. It is not a real
# commercial compound collection - it is a starter kit, so you can see the
# whole pipeline work end to end. Everything is written as SMILES, the
# standard text format for a molecule.
#
# A "[*:1]" in a SMILES string is an attachment point: "glue a group on here".
# Scaffolds have one or two attachment points. Groups have exactly one.
# ---------------------------------------------------------------------------

# Ring cores. These supply the aromatic / flat / stacking character and
# give the molecule a rigid shape.
SCAFFOLDS = [
    {"name": "benzene",              "smiles": "[*:1]c1ccccc1"},
    {"name": "phenol",               "smiles": "[*:1]c1ccc(O)cc1"},
    {"name": "para-disubst-benzene", "smiles": "[*:1]c1ccc([*:2])cc1"},
    {"name": "meta-disubst-benzene", "smiles": "[*:1]c1cccc([*:2])c1"},
    {"name": "pyridine",             "smiles": "[*:1]c1ccncc1"},
    {"name": "indole",               "smiles": "[*:1]c1cc2ccccc2[nH]1"},
    {"name": "naphthalene",          "smiles": "[*:1]c1ccc2ccccc2c1"},
    {"name": "biphenyl",             "smiles": "[*:1]c1ccc(-c2ccccc2)cc1"},
    {"name": "benzofuran",           "smiles": "[*:1]c1cc2ccccc2o1"},
    {"name": "quinoline",            "smiles": "[*:1]c1ccc2ccccc2n1"},
    {"name": "imidazole",            "smiles": "[*:1]c1cnc[nH]1"},
    {"name": "cyclohexane",          "smiles": "[*:1]C1CCCCC1"},
    {"name": "piperidine-ring",      "smiles": "[*:1]N1CCCCC1"},
    {"name": "indole-phenol",        "smiles": "[*:1]c1cc2cc(O)ccc2[nH]1"},
]

# Substituent groups. Each is tagged with the pharmacophore feature
# classes it provides, so the report can explain itself.
GROUPS = [
    {"name": "ethylamine",        "smiles": "[*:1]CCN",                 "provides": ["cationic", "hbond_donor"]},
    {"name": "propylamine",       "smiles": "[*:1]CCCN",                "provides": ["cationic", "hbond_donor"]},
    {"name": "dimethylaminoethyl","smiles": "[*:1]CCN(C)C",             "provides": ["cationic"]},
    {"name": "piperidinylethyl",  "smiles": "[*:1]CCN1CCCCC1",          "provides": ["cationic", "hydrophobic"]},
    {"name": "guanidinomethyl",   "smiles": "[*:1]CNC(=N)N",            "provides": ["cationic", "hbond_donor"]},
    {"name": "acetic-acid",       "smiles": "[*:1]CC(=O)O",             "provides": ["anionic", "hbond_acceptor"]},
    {"name": "propanoic-acid",    "smiles": "[*:1]CCC(=O)O",            "provides": ["anionic", "hbond_acceptor"]},
    {"name": "tetrazolylmethyl",  "smiles": "[*:1]Cc1nnn[nH]1",         "provides": ["anionic", "aromatic"]},
    {"name": "hydroxyl",          "smiles": "[*:1]O",                   "provides": ["hbond_donor", "hbond_acceptor"]},
    {"name": "hydroxyethyl",      "smiles": "[*:1]CCO",                 "provides": ["hbond_donor", "hbond_acceptor"]},
    {"name": "methoxy",           "smiles": "[*:1]OC",                  "provides": ["hbond_acceptor"]},
    {"name": "amide",             "smiles": "[*:1]C(=O)N",             "provides": ["hbond_donor", "hbond_acceptor"]},
    {"name": "sulfonamide",       "smiles": "[*:1]S(=O)(=O)N",          "provides": ["hbond_donor", "hbond_acceptor"]},
    {"name": "isobutyl",          "smiles": "[*:1]CC(C)C",              "provides": ["hydrophobic"]},
    {"name": "tert-butyl",        "smiles": "[*:1]C(C)(C)C",            "provides": ["hydrophobic"]},
    {"name": "cyclopropyl",       "smiles": "[*:1]C1CC1",               "provides": ["hydrophobic"]},
    {"name": "phenyl",            "smiles": "[*:1]c1ccccc1",            "provides": ["aromatic", "hydrophobic"]},
    {"name": "benzyl",            "smiles": "[*:1]Cc1ccccc1",           "provides": ["aromatic", "hydrophobic"]},
    {"name": "fluorophenyl",      "smiles": "[*:1]c1ccc(F)cc1",         "provides": ["aromatic", "hydrophobic"]},
    {"name": "methylsulfanyl",    "smiles": "[*:1]SC",                  "provides": ["hydrophobic"]},
    {"name": "hydrogen",          "smiles": "[*:1][H]",                 "provides": []},
]

# A few complete, well-known small molecules, included so the output is not
# purely combinatorial. These are public textbook structures.
KNOWN_FRAGMENTS = [
    {"name": "4-hydroxytamoxifen-core", "smiles": "CCC(=C(c1ccccc1)c1ccc(O)cc1)c1ccc(OCCN(C)C)cc1"},
    {"name": "raloxifene-like-fragment", "smiles": "Oc1ccc2cc(-c3ccc(OCCN4CCCCC4)cc3)sc2c1"},
    {"name": "genistein",                "smiles": "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12"},
    {"name": "diethylstilbestrol",       "smiles": "CCC(=C(CC)c1ccc(O)cc1)c1ccc(O)cc1"},
    {"name": "benzamidine",              "smiles": "NC(=N)c1ccccc1"},
    {"name": "tryptamine",               "smiles": "NCCc1c[nH]c2ccccc12"},
    {"name": "tyramine",                 "smiles": "NCCc1ccc(O)cc1"},
    {"name": "indole-3-acetic-acid",     "smiles": "OC(=O)Cc1c[nH]c2ccccc12"},
    {"name": "phenylalanine-like",       "smiles": "NC(Cc1ccccc1)C(=O)O"},
]


# (FEATURE_ALIASES, load_pharmacophore and summarise_demands live in
# pharmacophore_score.py)


# ---------------------------------------------------------------------------
# SECTION 4: Building the candidate molecules
# ---------------------------------------------------------------------------

def relabel_attachment(smiles, new_label):
    """Change a group's attachment point from [*:1] to [*:<new_label>]."""
    return smiles.replace("[*:1]", "[*:%d]" % new_label)


def join_fragments(piece_smiles_list):
    """Glue SMILES pieces together at their matching [*:n] attachment points.

    RDKit's Chem.molzip does the actual work: it finds dummy atoms that share
    a map number and fuses them into a real bond.
    """
    mols = []
    for smi in piece_smiles_list:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            return None
        mols.append(mol)
    combined = mols[0]
    for extra in mols[1:]:
        combined = Chem.CombineMols(combined, extra)
    try:
        zipped = Chem.molzip(combined)
    except Exception:
        return None
    if zipped is None:
        return None
    # molzip hands back a molecule that has not been "sanitized" (RDKit has
    # not yet worked out rings/aromaticity). Round-tripping through SMILES is
    # the simplest reliable way to get a clean, fully-initialised molecule.
    try:
        clean = Chem.MolFromSmiles(Chem.MolToSmiles(zipped))
    except Exception:
        return None
    return clean


def build_candidate_library():
    """Make every scaffold x group(s) combination, plus the known fragments.

    Returns a list of dicts: {"name", "smiles", "mol", "origin"}.
    Duplicates (the same molecule reached two different ways) are removed by
    comparing canonical SMILES.
    """
    candidates = []
    seen_smiles = set()

    def add(name, mol, origin):
        if mol is None:
            return
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return
        canonical = Chem.MolToSmiles(mol)
        if canonical in seen_smiles:
            return
        seen_smiles.add(canonical)
        candidates.append({"name": name, "smiles": canonical,
                           "mol": mol, "origin": origin})

    # (a) The hand-picked whole molecules.
    for frag in KNOWN_FRAGMENTS:
        add(frag["name"], Chem.MolFromSmiles(frag["smiles"]), "known fragment")

    # (b) Combinatorial scaffold + group assembly.
    for scaffold in SCAFFOLDS:
        two_points = "[*:2]" in scaffold["smiles"]
        for group_a in GROUPS:
            if not two_points:
                mol = join_fragments([scaffold["smiles"], group_a["smiles"]])
                add("%s + %s" % (scaffold["name"], group_a["name"]),
                    mol, "assembled")
            else:
                for group_b in GROUPS:
                    mol = join_fragments([
                        scaffold["smiles"],
                        group_a["smiles"],
                        relabel_attachment(group_b["smiles"], 2),
                    ])
                    add("%s + %s / %s" % (scaffold["name"], group_a["name"],
                                          group_b["name"]),
                        mol, "assembled")
    return candidates


# (feature counting, property measurement and scoring live in
# pharmacophore_score.py)


# (make_3d lives in pharmacophore_score.py)


# ---------------------------------------------------------------------------
# SECTION 7: Writing the report
# ---------------------------------------------------------------------------

def write_report(path, data, demands, results, args):
    lines = []
    add = lines.append
    code = data.get("structure", "unknown")

    add("=" * 72)
    add("LIGAND BUILDER REPORT - %s" % code)
    add("=" * 72)
    add("")
    add("Input pharmacophore file : %s" % args.pharmacophore_json)
    add("Active-site residues     : %s" % ", ".join(data.get("active_site_residues", [])))
    add("Pocket shell residues    : %d" % len(data.get("pocket_shell_residues", [])))
    add("Pocket max atom span     : %s A" % data.get("geometry", {}).get("max_atom_span_A"))
    add("Conformers requested     : %d per molecule" % args.confs)
    add("")

    add("-" * 72)
    add("1. WHAT THE POCKET ASKED FOR (top %d features, by priority)" % args.top)
    add("-" * 72)
    add("Each number is the total priority pocket_profiler.py gave that kind")
    add("of chemistry. Higher = the pocket wall asks for it more insistently.")
    add("")
    for cls, demand in sorted(demands.items(), key=lambda kv: -kv[1]):
        add("   %-16s demand weight %6.1f" % (cls, demand))
    add("")

    windows = data.get("suggested_ligand_properties", {})
    add("Suggested property windows taken from the profiler:")
    add("   molecular weight : %s Da" % (windows.get("molecular_weight_Da"),))
    add("   cLogP            : %s" % (windows.get("clogp"),))
    add("   H-bond donors    : %s" % (windows.get("hbond_donors"),))
    add("   H-bond acceptors : %s" % (windows.get("hbond_acceptors"),))
    add("   rotatable bonds  : max %s" % (windows.get("rotatable_bonds_max"),))
    add("")

    add("-" * 72)
    add("2. THE MOLECULES, BEST CHECKLIST MATCH FIRST")
    add("-" * 72)
    add("")
    for rank, item in enumerate(results, start=1):
        props = item["props"]
        sc = item["score_info"]
        add("#%-3d %s" % (rank, item["name"]))
        add("     SMILES          : %s" % item["smiles"])
        add("     source          : %s" % item["origin"])
        add("     checklist score : %.1f / 100   (chemistry %.0f%%, size-window fit %.0f%%)"
            % (sc["score"], 100 * sc["chemistry_fraction"], 100 * sc["property_fraction"]))
        add("     properties      : MW %.1f | cLogP %.2f | HBD %d | HBA %d | RotB %d | TPSA %.1f | charge %+d"
            % (props["mw"], props["clogp"], props["hbd"], props["hba"],
               props["rotb"], props["tpsa"], props["formal_charge"]))
        matched = ", ".join("%s x%d" % (k, v)
                            for k, v in sorted(sc["feature_counts"].items())
                            if v > 0)
        add("     features found  : %s" % (matched if matched else "none"))
        missing = ", ".join(k for k, v in sorted(sc["feature_counts"].items()) if v == 0)
        if missing:
            add("     features MISSING: %s" % missing)
        if sc["window_misses"]:
            add("     outside window  : %s" % "; ".join(sc["window_misses"]))
        add("     3D conformers   : %d written, minimised with %s"
            % (len(item["conf_ids"]), item["field"]))
        if item["energies"]:
            add("     force-field energy range: %.1f to %.1f kcal/mol "
                "(INTERNAL strain of the molecule alone - NOT binding energy)"
                % (min(item["energies"]), max(item["energies"])))
        add("")

    add("-" * 72)
    add("3. WHAT THE SCORE MEANS, AND WHAT IT DOES NOT")
    add("-" * 72)
    add("")
    add("The checklist score is built from exactly two things:")
    add("  (a) 70% - chemistry match. For each feature class the pocket asked")
    add("      for (cationic, aromatic, hydrophobic, ...), does this molecule")
    add("      contain a group of that kind? Found via RDKit SMARTS pattern")
    add("      matching. Having one counts for most of the credit; extra")
    add("      copies add a little more, with diminishing returns.")
    add("  (b) 30% - property-window fit. The fraction of the five suggested")
    add("      windows (MW, cLogP, HBD, HBA, rotatable bonds) the molecule")
    add("      falls inside.")
    add("")
    add("The score is therefore a measure of CHEMICAL SHOPPING-LIST COVERAGE.")
    add("")
    add("It is NOT:")
    add("  - a binding affinity, Kd, Ki or IC50;")
    add("  - a docking score;")
    add("  - an interaction energy with the protein;")
    add("  - evidence that the molecule fits in the pocket at all.")
    add("")
    add("Nothing in this script ever puts a molecule inside the protein. The")
    add("3D coordinates are the molecule's own relaxed shape in a vacuum, not")
    add("a pose in the binding site. The geometry of the pharmacophore (the")
    add("3D anchor points in the JSON) is NOT used to place atoms - only the")
    add("list of feature TYPES is used. Two molecules with identical scores")
    add("can differ a thousand-fold in real affinity, and a high-scoring")
    add("molecule may be sterically unable to enter the pocket.")
    add("")
    add("The force-field energies are internal strain of the isolated")
    add("molecule, useful only for comparing conformers of the SAME molecule.")
    add("")
    add("The candidate library is a small hand-written starter kit (%d scaffolds,"
        % len(SCAFFOLDS))
    add("%d substituent groups, %d known fragments), not a real screening" % (len(GROUPS), len(KNOWN_FRAGMENTS)))
    add("collection. Many of the assembled molecules may be hard or impossible")
    add("to synthesise; no synthetic accessibility check is done.")
    add("")
    add("HONEST NEXT STEP: feed the SDF file to a docking program")
    add("(AutoDock Vina or smina) together with the protein structure. That is")
    add("the first step that actually tests fit. Even then, docking scores are")
    add("weak predictors and only an experiment settles it.")
    add("")

    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# SECTION 8: Command line
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Turn a pocket_profiler.py pharmacophore JSON into real "
                    "3D candidate molecules (SDF) plus a plain-English report.")
    ap.add_argument("pharmacophore_json",
                    help="Path to <CODE>_pharmacophore.json from pocket_profiler.py")
    ap.add_argument("--outdir", default="output",
                    help="Output folder (default: output)")
    ap.add_argument("--top", type=int, default=12,
                    help="How many pharmacophore features to take as the "
                         "chemistry shopping list (default: 12)")
    ap.add_argument("--molecules", type=int, default=10,
                    help="How many best-scoring molecules to build in 3D "
                         "(default: 10)")
    ap.add_argument("--confs", type=int, default=10,
                    help="3D conformers to generate per molecule (default: 10)")
    ap.add_argument("--seed", type=int, default=0xF00D,
                    help="Random seed for conformer embedding (default: 61453)")
    args = ap.parse_args()

    print("=" * 72)
    print("LIGAND BUILDER")
    print("=" * 72)

    # --- step 1: read the blueprint -------------------------------------
    data = load_pharmacophore(args.pharmacophore_json)
    code = data.get("structure", "ligands")
    features = data.get("pharmacophore_features", [])
    windows = data.get("suggested_ligand_properties", {})
    print("Read pharmacophore for %s: %d features, %d pocket residues."
          % (code, len(features), len(data.get("pocket_shell_residues", []))))

    demands = summarise_demands(features, args.top)
    if not demands:
        sys.exit("ERROR: the top %d features contain no usable chemistry "
                 "classes. Try a larger --top." % args.top)
    print("Chemistry the pocket asks for (top %d features):" % args.top)
    for cls, demand in sorted(demands.items(), key=lambda kv: -kv[1]):
        print("   %-16s weight %6.1f" % (cls, demand))

    # --- step 2: build and score candidates -----------------------------
    print("")
    print("Building candidate library...")
    candidates = build_candidate_library()
    print("   %d unique candidate molecules assembled." % len(candidates))

    print("Scoring candidates against the pharmacophore...")
    scored = []
    for cand in candidates:
        try:
            info = score_candidate(cand["mol"], demands, windows)
        except Exception as err:
            print("   skipped %s (%s)" % (cand["name"], err))
            continue
        cand["score_info"] = info
        cand["props"] = measure(cand["mol"])
        scored.append(cand)
    scored.sort(key=lambda c: -c["score_info"]["score"])
    print("   %d candidates scored. Best checklist score: %.1f / 100."
          % (len(scored), scored[0]["score_info"]["score"]))

    # --- step 3: 3D conformers for the winners --------------------------
    os.makedirs(args.outdir, exist_ok=True)
    sdf_path = os.path.join(args.outdir, "%s_ligands.sdf" % code)
    txt_path = os.path.join(args.outdir, "%s_ligand_report.txt" % code)

    print("")
    print("Generating 3D conformers for the top %d molecules..." % args.molecules)
    writer = Chem.SDWriter(sdf_path)
    kept = []
    total_confs = 0
    # Keep the output varied: at most MAX_PER_SCAFFOLD molecules built from
    # any single ring core, otherwise the top of the list is ten near-twins.
    MAX_PER_SCAFFOLD = 2
    used_scaffold = {}
    for cand in scored:
        if len(kept) >= args.molecules:
            break
        core = cand["name"].split(" + ")[0]
        if used_scaffold.get(core, 0) >= MAX_PER_SCAFFOLD:
            continue
        mol_h, conf_ids, energies, field = make_3d(cand["mol"], args.confs,
                                                   seed=args.seed)
        if mol_h is None:
            print("   %-45s SKIPPED (%s)" % (cand["name"][:45], field))
            continue
        # Sort this molecule's conformers from lowest to highest strain.
        order = sorted(range(len(conf_ids)), key=lambda i: energies[i])
        sc = cand["score_info"]
        for position, idx in enumerate(order, start=1):
            conf_id = conf_ids[idx]
            mol_h.SetProp("_Name", "%s_conf%d" % (cand["name"].replace(" ", "_"), position))
            mol_h.SetProp("molecule_name", cand["name"])
            mol_h.SetProp("smiles", cand["smiles"])
            mol_h.SetProp("source", cand["origin"])
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
        cand["conf_ids"] = conf_ids
        cand["energies"] = energies
        cand["field"] = field
        used_scaffold[core] = used_scaffold.get(core, 0) + 1
        kept.append(cand)
        print("   %-45s score %5.1f  %2d confs  %s"
              % (cand["name"][:45], sc["score"], len(conf_ids), field))
    writer.close()

    if not kept:
        sys.exit("ERROR: no molecule could be embedded in 3D.")

    # --- step 4: report --------------------------------------------------
    write_report(txt_path, data, demands, kept, args)

    print("")
    print("-" * 72)
    print("DONE.")
    print("  molecules written : %d" % len(kept))
    print("  conformers written: %d" % total_confs)
    print("  SDF (3D molecules): %s" % sdf_path)
    print("  report (plain text): %s" % txt_path)
    print("")
    print("REMINDER: the checklist score measures how well a molecule's")
    print("chemistry matches the pocket's shopping list. It is NOT a binding")
    print("affinity and nothing here has been docked into the protein.")
    print("-" * 72)


if __name__ == "__main__":
    main()
