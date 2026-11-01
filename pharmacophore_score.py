#!/usr/bin/env python3
"""
pharmacophore_score.py
======================

Shared chemistry code for the Alpha-Ligand workflow.

This module holds the pieces that BOTH `ligand_builder.py` (the small
hand-written starter library) and `screen_library.py` (screening a real
compound collection) need, so the scoring is defined in exactly one place
and the two paths cannot drift apart:

  * the SMARTS patterns used to recognise chemistry,
  * reading a pocket_profiler.py pharmacophore JSON,
  * turning its ranked features into a per-class "demand",
  * measuring a molecule's standard properties,
  * the 0-100 checklist score,
  * 3D conformer generation.

WHAT THE SCORE IS *NOT*: it is not a binding affinity, not a docking score
and not an energy. Nothing in this file puts a molecule inside a protein.
It only asks "does this molecule's chemistry tick the boxes the pocket
asked for, and is it roughly the right size".
"""

import json
import math
import os
import sys

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# SECTION 1: How we RECOGNISE chemistry, using SMARTS patterns
#
# SMARTS is a little pattern language for molecules, like a regular
# expression for chemical structure. Each pattern below asks a yes/no
# question such as "is there a basic amine in this molecule?".
# These are deliberately simple, readable patterns - not perfect ones.
# ---------------------------------------------------------------------------

SMARTS_PATTERNS = {
    # Groups that carry a POSITIVE charge at body pH: plain amines,
    # amidines, guanidines, and already-charged nitrogens.
    "cationic": [
        "[NX3;H2,H1;!$(N[!#6;!#1]);!$(N-C=[O,N,S]);!$(N-a)]",  # 1deg/2deg aliphatic amine
        "[NX3;H0;!$(N[!#6;!#1]);!$(N-C=[O,N,S]);!$(N-a);!R]",  # 3deg aliphatic amine
        "[$([NX3][CX3]=[NX2]),$([NX3][CX3](=[NX2])[NX3])]",    # amidine / guanidine
        "[N+;!$([N+][O-])]",                                    # permanent cation
    ],
    # Groups that carry a NEGATIVE charge at body pH.
    "anionic": [
        "[CX3](=O)[OX2H1]",            # carboxylic acid
        "[CX3](=O)[OX1-]",             # carboxylate
        "c1nnn[nH]1",                  # tetrazole (a carboxylate mimic)
        "[SX4](=O)(=O)[OX2H1,OX1-]",   # sulfonic acid
        "[PX4](=O)([OX2H1,OX1-])",     # phosphonate
    ],
    # Hydrogen-bond donors and acceptors (fine detail is handled by
    # RDKit's own counters further down; these just flag presence).
    "hbond_donor": ["[N;!H0]", "[O;H1]", "[S;H1]"],
    "hbond_acceptor": ["[O;!$(O=*);!$([O-])]", "[o]", "[nX2;!H1]", "[O]=[C,S,P]"],
}

# Pre-compile the patterns once (faster, and it fails loudly if a pattern
# has a typo in it).
COMPILED_PATTERNS = {}
for feature_name, smarts_list in SMARTS_PATTERNS.items():
    COMPILED_PATTERNS[feature_name] = []
    for smarts in smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            raise ValueError("Bad SMARTS pattern in the table: %s" % smarts)
        COMPILED_PATTERNS[feature_name].append(patt)

# ---------------------------------------------------------------------------
# SECTION 3: Reading the pharmacophore JSON
# ---------------------------------------------------------------------------

# pocket_profiler.py uses some combined feature names. Map them onto the
# simple classes this script knows how to detect.
FEATURE_ALIASES = {
    "cationic": ["cationic"],
    "anionic": ["anionic"],
    "aromatic": ["aromatic"],
    "hydrophobic": ["hydrophobic"],
    "hbond_donor": ["hbond_donor"],
    "hbond_acceptor": ["hbond_acceptor"],
    # Combined types the profiler uses: either/both chemistries would work.
    "anionic_or_aromatic": ["anionic", "aromatic"],
    "aromatic_and_hbond": ["aromatic", "hbond_acceptor"],
    "hbond_mixed": ["hbond_acceptor", "hbond_donor"],
    "hbond_acceptor_or_covalent": ["hbond_acceptor"],
    # Types that describe SPACE rather than chemistry. We cannot check these
    # with a substructure pattern, so they are deliberately ignored here.
    "shape": [],
    "open_space": [],
}


def load_pharmacophore(json_path):
    """Read the JSON file written by pocket_profiler.py."""
    if not os.path.isfile(json_path):
        sys.exit("ERROR: no such file: %s\n"
                 "Run pocket_profiler.py first to create it." % json_path)
    with open(json_path) as handle:
        data = json.load(handle)
    for required_key in ("pharmacophore_features", "suggested_ligand_properties"):
        if required_key not in data:
            sys.exit("ERROR: %s does not look like a pocket_profiler.py "
                     "pharmacophore file (missing '%s')." % (json_path, required_key))
    return data


def summarise_demands(features, top_n):
    """Turn the top-ranked pharmacophore features into a 'demand' per class.

    Returns a dict like {"cationic": 15.0, "aromatic": 24.0, ...} where the
    number is the total priority score pocket_profiler.py assigned to that
    kind of feature. Bigger number = the pocket asks for it more loudly.
    """
    demands = {}
    for feature in features[:top_n]:
        raw_type = feature.get("ligand_feature_type", "")
        score = float(feature.get("priority_score", 0.0))
        if raw_type not in FEATURE_ALIASES:
            # An unrecognised feature type: say so rather than silently
            # scoring every molecule against a pattern we do not have.
            print("   note: ignoring unrecognised feature type '%s'" % raw_type)
            continue
        classes = FEATURE_ALIASES[raw_type]
        if not classes:
            continue
        # A combined type like "anionic_or_aromatic" splits its score, since
        # either chemistry would do.
        share = score / float(len(classes))
        for cls in classes:
            demands[cls] = demands.get(cls, 0.0) + share
    return demands

# ---------------------------------------------------------------------------
# SECTION 5: Measuring each candidate
# ---------------------------------------------------------------------------

def count_feature(mol, feature_name):
    """How many times does this molecule show a given pharmacophore feature?"""
    if feature_name == "aromatic":
        return rdMolDescriptors.CalcNumAromaticRings(mol)
    if feature_name == "hydrophobic":
        # Carbon atoms with no polar neighbour: a crude greasiness count.
        patt = Chem.MolFromSmarts("[#6;!$([#6]~[#7,#8,#9,#15,#16])]")
        return len(mol.GetSubstructMatches(patt))
    total = 0
    for patt in COMPILED_PATTERNS.get(feature_name, []):
        total += len(mol.GetSubstructMatches(patt))
    return total


def measure(mol):
    """Collect the standard drug-like property numbers RDKit can calculate."""
    return {
        "mw": Descriptors.MolWt(mol),
        "clogp": Crippen.MolLogP(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
    }


def property_fit(props, windows):
    """Fraction (0-1) of the pocket's suggested property windows this molecule
    sits inside, plus a list of the ones it misses, in plain English."""
    checks = []

    def window_check(label, value, low, high, fmt="%.1f"):
        inside = (low <= value <= high)
        checks.append((label, inside,
                       ("%s " + fmt + " (want " + fmt + "-" + fmt + ")")
                       % (label, value, low, high)))
        return inside

    mw_low, mw_high = windows.get("molecular_weight_Da", [0, 10000])
    window_check("MW", props["mw"], mw_low, mw_high)

    lp_low, lp_high = windows.get("clogp", [-10, 10])
    window_check("cLogP", props["clogp"], lp_low, lp_high, "%.2f")

    hbd_low, hbd_high = windows.get("hbond_donors", [0, 99])
    window_check("HBD", props["hbd"], hbd_low, hbd_high, "%.0f")

    hba_low, hba_high = windows.get("hbond_acceptors", [0, 99])
    window_check("HBA", props["hba"], hba_low, hba_high, "%.0f")

    rot_max = windows.get("rotatable_bonds_max", 99)
    window_check("RotB", props["rotb"], 0, rot_max, "%.0f")

    passed = sum(1 for _, ok, _ in checks if ok)
    misses = [text for _, ok, text in checks if not ok]
    return passed / float(len(checks)), misses, checks


def score_candidate(mol, demands, windows):
    """Score how well one molecule matches the pharmacophore. 0-100.

    Chemistry match (70% of the score): for every feature class the pocket
    asks for, we check whether the molecule actually has that group. Having
    one counts for most of the credit; having several adds a little more
    (diminishing returns, via a square root), because a pocket only has so
    many places to use them.

    Property fit (30% of the score): the fraction of the suggested size /
    greasiness / H-bond windows the molecule falls inside.

    THIS IS A CHECKLIST SCORE, NOT AN AFFINITY.
    """
    counts = {}
    earned = 0.0
    possible = 0.0
    for cls, demand in demands.items():
        n = count_feature(mol, cls)
        counts[cls] = n
        possible += demand
        if n > 0:
            # 1 copy -> 100% of credit for the first unit, extra copies add
            # less and less, capped at 1.0 of the demand.
            credit = min(1.0, 0.75 + 0.25 * math.sqrt(min(n, 4) - 1) / math.sqrt(3))
            earned += demand * credit
    chemistry = (earned / possible) if possible > 0 else 0.0

    fit, misses, checks = property_fit(measure(mol), windows)
    total = 100.0 * (0.70 * chemistry + 0.30 * fit)
    return {
        "score": total,
        "chemistry_fraction": chemistry,
        "property_fraction": fit,
        "feature_counts": counts,
        "window_misses": misses,
        "window_checks": checks,
    }

# ---------------------------------------------------------------------------
# SECTION 6: 3D conformer generation
# ---------------------------------------------------------------------------

def make_3d(mol, n_confs, seed=0xF00D):
    """Turn a flat molecule into real 3D coordinates and tidy them up.

    Steps:
      1. Add explicit hydrogens (3D geometry needs them).
      2. EmbedMultipleConfs: RDKit's ETKDGv3 distance-geometry method guesses
         several different 3D shapes ("conformers").
      3. Minimise each shape with the MMFF94 force field, or UFF when MMFF
         has no parameters for an atom type.

    Returns (molecule_with_Hs, list_of_conformer_ids, energies, force_field_used)
    or (None, [], [], reason) if the embedding failed.
    """
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.5          # throw away near-identical shapes
    params.useSmallRingTorsions = True
    conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs,
                                               params=params))
    if not conf_ids:
        # Fall back to a less fussy embedder before giving up.
        params.useRandomCoords = True
        conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs,
                                                   params=params))
    if not conf_ids:
        return None, [], [], "embedding failed"

    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        results = AllChem.MMFFOptimizeMoleculeConfs(mol_h, maxIters=2000)
        field = "MMFF94"
    else:
        results = AllChem.UFFOptimizeMoleculeConfs(mol_h, maxIters=2000)
        field = "UFF"
    # results is a list of (converged_flag, energy) pairs, one per conformer.
    energies = [energy for _converged, energy in results]
    return mol_h, conf_ids, energies, field
