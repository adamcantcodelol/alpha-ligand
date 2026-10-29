#!/usr/bin/env python3
"""
pocket_profiler.py
==================

A beginner-friendly, structure-based "work backwards" ligand design helper.

WHAT IT DOES (in plain English)
-------------------------------
1. You give it a protein structure (a 4-letter PDB code that it downloads from
   RCSB, or a local .pdb / .cif file) plus the active-site residues you care
   about (e.g. "A:45,A:88,A:102").
2. It finds every residue that has an atom within a cutoff distance (8 A by
   default) of the side-chain atoms of those residues. That set of residues is
   the "pocket shell" - the wall of the binding site.
3. It measures the chemistry of that wall: how hydrophobic it is, how charged,
   how many hydrogen-bond donors and acceptors it offers, how big the pocket
   is, and (if your Biopython has Shrake-Rupley) how exposed it is.
4. It then flips every pocket feature into the COMPLEMENTARY ligand feature.
   A positively charged wall wants a negatively charged ligand group, an
   aromatic wall wants a ring to stack against, and so on. That list of
   complementary features, each with an approximate 3D position, is a
   "pharmacophore": a blueprint of what a high-affinity repressor ligand
   should look like.

IMPORTANT: this is a HYPOTHESIS GENERATOR. It does no docking, no energy
calculation and no validation. It tells you what to go look for; it does not
tell you that any molecule will bind.

Only Biopython is used.
"""

import argparse
import json
import os
import sys
from itertools import combinations

from Bio.PDB import PDBList, PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa

# ---------------------------------------------------------------------------
# SECTION 1: Chemistry lookup tables
# These are plain dictionaries describing each amino acid. Nothing clever here;
# they are just the facts we need later.
# ---------------------------------------------------------------------------

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Kyte-Doolittle hydropathy: positive = water-hating, negative = water-loving.
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

HYDROPHOBIC = set("AVLIMFWCP")
AROMATIC = set("FWYH")
POSITIVE = set("KR")          # formally + at pH 7
NEGATIVE = set("DE")          # formally - at pH 7
POLAR_NEUTRAL = set("STNQCYH")

# How many H-bond donors / acceptors each SIDE CHAIN offers (backbone ignored).
SIDECHAIN_DONORS = {"R": 3, "K": 1, "W": 1, "N": 1, "Q": 1,
                    "S": 1, "T": 1, "Y": 1, "H": 1, "C": 1}
SIDECHAIN_ACCEPTORS = {"D": 2, "E": 2, "N": 1, "Q": 1,
                       "S": 1, "T": 1, "Y": 1, "H": 1}

# Atoms that make up each side chain (used to anchor pharmacophore points).
KEY_SIDECHAIN_ATOMS = {
    "R": ["CZ", "NH1", "NH2"], "K": ["NZ"], "D": ["CG", "OD1", "OD2"],
    "E": ["CD", "OE1", "OE2"], "F": ["CZ"], "Y": ["OH", "CZ"],
    "W": ["CZ2", "NE1"], "H": ["NE2", "ND1"], "S": ["OG"], "T": ["OG1"],
    "N": ["OD1", "ND2"], "Q": ["OE1", "NE2"], "C": ["SG"], "M": ["SD"],
    "L": ["CD1"], "I": ["CD1"], "V": ["CG1"], "A": ["CB"], "P": ["CG"],
    "G": ["CA"],
}


# ---------------------------------------------------------------------------
# SECTION 2: Small helpers
# ---------------------------------------------------------------------------

def one_letter(residue):
    """Turn a Biopython residue into its one-letter code, or None if it is
    not one of the 20 standard amino acids (water, ions, ligands -> None)."""
    return THREE_TO_ONE.get(residue.get_resname().strip().upper())


def sidechain_atoms(residue):
    """Return the side-chain atoms of a residue (everything that is not
    backbone N/CA/C/O). Glycine has no side chain, so we fall back to CA."""
    backbone = {"N", "CA", "C", "O", "OXT"}
    atoms = [a for a in residue if a.get_id() not in backbone]
    if not atoms:
        atoms = [a for a in residue if a.get_id() == "CA"]
    return atoms


def residue_key(chain_id, residue):
    """A readable label like 'A:45' (insertion codes kept if present)."""
    hetflag, resseq, icode = residue.get_id()
    return "%s:%d%s" % (chain_id, resseq, icode.strip())


def parse_residue_spec(spec, default_chain):
    """Turn the command-line residue list into [(chain, number), ...].

    Accepts e.g.  "A:45,A:88,102"  or  "45,88,102" with --chain A.
    """
    wanted = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" in token:
            chain, number = token.split(":", 1)
        else:
            chain, number = default_chain, token
        if chain is None:
            raise ValueError(
                "Residue '%s' has no chain. Either write it as A:%s or pass "
                "--chain A." % (token, token))
        try:
            wanted.append((chain.strip(), int(number.strip())))
        except ValueError:
            raise ValueError("Could not read residue number from '%s'." % token)
    if not wanted:
        raise ValueError("No residues given.")
    return wanted


# ---------------------------------------------------------------------------
# SECTION 3: Loading the structure
# ---------------------------------------------------------------------------

def load_structure(source, workdir):
    """Load a structure from a 4-character PDB code or a local file.

    Returns (structure_object, short_name).
    """
    # Case 1: it is an existing file on disk.
    if os.path.isfile(source):
        name = os.path.splitext(os.path.basename(source))[0]
        if source.lower().endswith((".cif", ".mmcif")):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        try:
            return parser.get_structure(name, source), name
        except Exception as exc:
            raise SystemExit("ERROR: could not parse '%s': %s" % (source, exc))

    # Case 2: it looks like a PDB code -> download it from RCSB.
    code = source.strip().lower()
    if len(code) != 4 or not code[0].isdigit():
        raise SystemExit(
            "ERROR: '%s' is neither an existing file nor a valid 4-character "
            "PDB code (codes look like 1LMB or 3ERT)." % source)

    os.makedirs(workdir, exist_ok=True)
    pdbl = PDBList(verbose=False)
    try:
        path = pdbl.retrieve_pdb_file(code, pdir=workdir, file_format="mmCif")
    except Exception as exc:
        raise SystemExit("ERROR: download of '%s' failed: %s" % (source, exc))
    if not path or not os.path.isfile(path):
        raise SystemExit(
            "ERROR: PDB code '%s' could not be downloaded. Check the code and "
            "your internet connection." % source)
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(code, path), code.upper()


def first_model(structure):
    """Structures can hold several NMR models; we simply use the first."""
    for model in structure:
        return model
    raise SystemExit("ERROR: the structure contains no models.")


# ---------------------------------------------------------------------------
# SECTION 4: Selecting the active site and the pocket shell
# ---------------------------------------------------------------------------

def collect_active_site(model, wanted):
    """Find the residues the user asked for. Errors clearly if missing."""
    found = []
    for chain_id, number in wanted:
        if chain_id not in [c.get_id() for c in model]:
            raise SystemExit(
                "ERROR: chain '%s' is not in this structure. Chains present: %s"
                % (chain_id, ", ".join(c.get_id() for c in model)))
        chain = model[chain_id]
        hit = None
        for residue in chain:
            if residue.get_id()[1] == number and is_aa(residue, standard=False):
                hit = residue
                break
        if hit is None:
            raise SystemExit(
                "ERROR: residue %d was not found in chain %s." % (number, chain_id))
        if one_letter(hit) is None:
            raise SystemExit(
                "ERROR: residue %s:%d is '%s', not a standard amino acid."
                % (chain_id, number, hit.get_resname()))
        found.append((chain_id, hit))
    return found


def collect_pocket_shell(model, active_site, cutoff):
    """Every standard amino acid with any atom within `cutoff` angstroms of an
    active-site side-chain atom. This is the wall of the pocket."""
    probe_atoms = []
    for _, residue in active_site:
        probe_atoms.extend(sidechain_atoms(residue))

    shell = {}
    for chain in model:
        for residue in chain:
            if one_letter(residue) is None:
                continue  # skip water, ions, bound ligands
            for atom in residue:
                if any((atom - probe) <= cutoff for probe in probe_atoms):
                    shell[residue_key(chain.get_id(), residue)] = (chain.get_id(), residue)
                    break
    return shell


# ---------------------------------------------------------------------------
# SECTION 5: Profiling the pocket chemistry and geometry
# ---------------------------------------------------------------------------

def profile_pocket(shell):
    """Count the chemical character of the pocket residues."""
    letters = [one_letter(res) for _, res in shell.values()]
    n = len(letters)

    def frac(group):
        hits = [c for c in letters if c in group]
        return len(hits), (len(hits) / n if n else 0.0)

    hydrophobic_n, hydrophobic_f = frac(HYDROPHOBIC)
    aromatic_n, aromatic_f = frac(AROMATIC)
    positive_n, positive_f = frac(POSITIVE)
    negative_n, negative_f = frac(NEGATIVE)
    polar_n, polar_f = frac(POLAR_NEUTRAL)

    composition = {}
    for c in letters:
        composition[c] = composition.get(c, 0) + 1

    return {
        "n_residues": n,
        "composition": dict(sorted(composition.items(), key=lambda kv: -kv[1])),
        "hydrophobic": {"count": hydrophobic_n, "fraction": round(hydrophobic_f, 3)},
        "aromatic": {"count": aromatic_n, "fraction": round(aromatic_f, 3)},
        "positive": {"count": positive_n, "fraction": round(positive_f, 3)},
        "negative": {"count": negative_n, "fraction": round(negative_f, 3)},
        "polar_neutral": {"count": polar_n, "fraction": round(polar_f, 3)},
        "mean_kd_hydropathy": round(
            sum(KD_HYDROPATHY[c] for c in letters) / n, 3) if n else 0.0,
        "net_formal_charge_pH7": positive_n - negative_n,
        "sidechain_hbond_donors": sum(SIDECHAIN_DONORS.get(c, 0) for c in letters),
        "sidechain_hbond_acceptors": sum(SIDECHAIN_ACCEPTORS.get(c, 0) for c in letters),
    }


def pocket_geometry(shell):
    """Bounding box, centroid and the longest distance between any two pocket
    atoms. A rough but honest measure of how big the site is."""
    coords = []
    for _, residue in shell.values():
        for atom in residue:
            coords.append(atom.get_coord())
    if not coords:
        return {}

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    centroid = [round(float(sum(xs) / len(xs)), 2),
                round(float(sum(ys) / len(ys)), 2),
                round(float(sum(zs) / len(zs)), 2)]
    box = [round(float(max(xs) - min(xs)), 2),
           round(float(max(ys) - min(ys)), 2),
           round(float(max(zs) - min(zs)), 2)]

    # Longest span. For big pockets we subsample so this stays fast.
    sample = coords if len(coords) <= 400 else coords[::max(1, len(coords) // 400)]
    span = 0.0
    for a, b in combinations(sample, 2):
        d = sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5
        span = max(span, float(d))

    return {
        "centroid": centroid,
        "bounding_box_xyz": box,
        "max_atom_span_A": round(float(span), 2),
        "n_pocket_atoms": len(coords),
    }


def solvent_accessibility(model, shell):
    """Use Shrake-Rupley if the INSTALLED Biopython provides it. We check for
    it rather than assuming, because older versions do not have it."""
    try:
        from Bio.PDB import ShrakeRupley
    except ImportError:
        return {"available": False,
                "note": "ShrakeRupley not present in this Biopython version."}
    try:
        sr = ShrakeRupley()
        sr.compute(model, level="R")
    except Exception as exc:
        return {"available": False, "note": "Shrake-Rupley failed: %s" % exc}

    per_residue = {}
    for key, (_, residue) in shell.items():
        sasa = getattr(residue, "sasa", None)
        if sasa is not None:
            per_residue[key] = round(float(sasa), 1)
    total = round(sum(per_residue.values()), 1)
    return {"available": True,
            "total_pocket_sasa_A2": total,
            "per_residue_sasa_A2": per_residue}


# ---------------------------------------------------------------------------
# SECTION 6: Working backwards -> the complementary pharmacophore
# ---------------------------------------------------------------------------

def anchor_point(residue):
    """A representative 3D point for a side chain: the average position of its
    chemically important atoms (falls back to the whole side chain)."""
    letter = one_letter(residue)
    names = KEY_SIDECHAIN_ATOMS.get(letter, [])
    atoms = [a for a in residue if a.get_id() in names] or sidechain_atoms(residue)
    n = len(atoms)
    coord = [sum(a.get_coord()[i] for a in atoms) / n for i in range(3)]
    # round(float(...)) keeps these as plain Python floats so they can be
    # written to JSON (numpy float32 cannot be).
    return [round(float(v), 2) for v in coord]


# For each pocket residue type: what the LIGAND should bring to face it.
COMPLEMENT_RULES = {
    "R": ("anionic", "carboxylate / phosphate / sulfonate",
          "Arg guanidinium: salt bridge + bidentate H-bond", 10),
    "K": ("anionic", "carboxylate / phosphate",
          "Lys ammonium: salt bridge", 10),
    "H": ("anionic_or_aromatic", "carboxylate, or aryl ring for stacking",
          "His can be + at pH 7 and also stacks", 7),
    "D": ("cationic", "amine / amidine / guanidine",
          "Asp carboxylate: salt bridge", 10),
    "E": ("cationic", "amine / amidine / guanidine",
          "Glu carboxylate: salt bridge", 10),
    "F": ("aromatic", "aryl / heteroaryl ring",
          "Phe ring: pi-pi stacking", 8),
    "W": ("aromatic", "aryl ring or cation for cation-pi",
          "Trp indole: stacking / cation-pi", 9),
    "Y": ("aromatic_and_hbond", "aryl ring plus an H-bond acceptor",
          "Tyr: stacking ring + donor OH", 9),
    "S": ("hbond_acceptor", "carbonyl / ester / ether / N of heteroaryl",
          "Ser OH donates an H-bond", 5),
    "T": ("hbond_acceptor", "carbonyl / ether oxygen",
          "Thr OH donates an H-bond", 5),
    "N": ("hbond_mixed", "H-bond acceptor facing the NH2, donor facing the C=O",
          "Asn amide donates and accepts", 6),
    "Q": ("hbond_mixed", "H-bond acceptor facing the NH2, donor facing the C=O",
          "Gln amide donates and accepts", 6),
    "C": ("hbond_acceptor_or_covalent", "H-bond acceptor; or a mild electrophile",
          "Cys thiol - only consider covalent chemistry deliberately", 4),
    "M": ("hydrophobic", "alkyl / thioether-friendly greasy group",
          "Met side chain is soft and greasy", 4),
    "L": ("hydrophobic", "alkyl group", "Leu greasy wall", 4),
    "I": ("hydrophobic", "alkyl group", "Ile greasy wall", 4),
    "V": ("hydrophobic", "small alkyl group", "Val greasy wall", 3),
    "A": ("hydrophobic", "methyl-sized contact", "Ala small greasy contact", 2),
    "P": ("hydrophobic", "rigid hydrophobic contact", "Pro rigid greasy contact", 2),
    "G": ("shape", "room to grow - Gly leaves space",
          "Gly has no side chain; the pocket is open here", 2),
}


def build_pharmacophore(active_site, shell):
    """Turn each pocket residue into the complementary ligand feature."""
    active_keys = {residue_key(cid, res) for cid, res in active_site}
    features = []
    for key, (chain_id, residue) in shell.items():
        letter = one_letter(residue)
        rule = COMPLEMENT_RULES.get(letter)
        if rule is None:
            continue
        feature_type, chemistry, why, base_score = rule
        # Residues the user named as the active site matter most, so we boost
        # their score. Everything else is pocket wall.
        score = base_score * (1.5 if key in active_keys else 1.0)
        features.append({
            "pocket_residue": "%s%s" % (letter, key),
            "pocket_residue_name": residue.get_resname(),
            "is_named_active_site": key in active_keys,
            "ligand_feature_type": feature_type,
            "suggested_ligand_chemistry": chemistry,
            "rationale": why,
            "approx_xyz": anchor_point(residue),
            "priority_score": round(score, 2),
        })
    features.sort(key=lambda f: -f["priority_score"])
    return features


def suggest_properties(profile, geometry):
    """Rough, honest ligand property windows derived from pocket size and
    chemistry. These are rules-of-thumb, not measurements."""
    span = geometry.get("max_atom_span_A", 15)
    # A ligand should span roughly half to two thirds of the pocket.
    target_len = max(6.0, span * 0.55)
    # ~ 15-20 Da of ligand per angstrom of length is a common rough scaling.
    mw_low = round(max(150, target_len * 16), -1)
    mw_high = round(min(650, target_len * 30), -1)
    if mw_high <= mw_low:
        mw_high = mw_low + 150

    kd = profile["mean_kd_hydropathy"]
    # Greasy pocket -> higher logP; polar pocket -> lower logP.
    logp_mid = 1.5 + kd  # kd roughly ranges -2 .. +2 for real pockets
    logp_low = round(max(-1.0, logp_mid - 1.5), 1)
    logp_high = round(min(5.0, logp_mid + 1.5), 1)

    # Ligand charge should oppose the pocket's net charge, but stay modest.
    net = profile["net_formal_charge_pH7"]
    ligand_charge = max(-3, min(3, -net // 2)) if net else 0

    return {
        "molecular_weight_Da": [mw_low, mw_high],
        "clogp": [logp_low, logp_high],
        "rotatable_bonds_max": max(3, min(10, int(target_len / 2.5))),
        "hbond_donors": [0, max(1, min(5, profile["sidechain_hbond_acceptors"] // 3))],
        "hbond_acceptors": [1, max(2, min(10, profile["sidechain_hbond_donors"] // 2))],
        "suggested_net_charge": ligand_charge,
        "basis": ("Derived from pocket span %.1f A, mean Kyte-Doolittle %.2f and "
                  "net pocket charge %+d. Rules of thumb only."
                  % (span, kd, net)),
    }


# ---------------------------------------------------------------------------
# SECTION 7: Reporting
# ---------------------------------------------------------------------------

DISCLAIMER = (
    "NOTE: This is a HYPOTHESIS-GENERATING pharmacophore derived purely from\n"
    "the geometry and chemistry of the pocket residues. Nothing here has been\n"
    "docked, scored, simulated or experimentally validated. Treat it as a\n"
    "starting blueprint for designing or searching for a repressor ligand."
)


def render_report(name, args, active_site, shell, profile, geometry, sasa,
                  features, properties):
    L = []
    add = L.append
    add("=" * 72)
    add("POCKET PROFILE AND COMPLEMENTARY LIGAND BLUEPRINT")
    add("Structure : %s" % name)
    add("Cutoff    : %.1f A around active-site side chains" % args.cutoff)
    add("=" * 72)
    add("")
    add("1. ACTIVE SITE RESIDUES YOU NAMED")
    for cid, res in active_site:
        add("   %s%s (%s)" % (one_letter(res), residue_key(cid, res), res.get_resname()))
    add("")
    add("2. POCKET SHELL (%d residues within %.1f A)" % (profile["n_residues"], args.cutoff))
    add("   " + ", ".join(sorted(shell.keys())))
    add("")
    add("3. CHEMICAL CHARACTER OF THE POCKET")
    add("   Composition          : " + ", ".join(
        "%s=%d" % (k, v) for k, v in profile["composition"].items()))
    for label, key in [("Hydrophobic", "hydrophobic"), ("Aromatic", "aromatic"),
                       ("Positive (K/R)", "positive"), ("Negative (D/E)", "negative"),
                       ("Polar neutral", "polar_neutral")]:
        d = profile[key]
        add("   %-21s: %d  (%.0f%%)" % (label, d["count"], 100 * d["fraction"]))
    add("   Mean KD hydropathy   : %+.2f  (%s)" % (
        profile["mean_kd_hydropathy"],
        "greasy pocket" if profile["mean_kd_hydropathy"] > 0 else "polar pocket"))
    add("   Net formal charge    : %+d at pH 7" % profile["net_formal_charge_pH7"])
    add("   Side-chain H-bond    : %d donors, %d acceptors" % (
        profile["sidechain_hbond_donors"], profile["sidechain_hbond_acceptors"]))
    add("")
    add("4. POCKET SIZE AND SHAPE")
    add("   Centroid (x,y,z)     : %s" % geometry.get("centroid"))
    add("   Bounding box (A)     : %s" % geometry.get("bounding_box_xyz"))
    add("   Max atom span (A)    : %s" % geometry.get("max_atom_span_A"))
    add("   Pocket atoms counted : %s" % geometry.get("n_pocket_atoms"))
    if sasa.get("available"):
        add("   Total pocket SASA    : %.1f A^2 (Shrake-Rupley)" % sasa["total_pocket_sasa_A2"])
    else:
        add("   SASA                 : not computed (%s)" % sasa.get("note"))
    add("")
    add("5. WORKING BACKWARDS: COMPLEMENTARY LIGAND FEATURES (ranked)")
    add("   %-12s %-22s %-34s %s" % ("POCKET RES", "LIGAND FEATURE", "SUGGESTED CHEMISTRY", "APPROX XYZ"))
    for f in features[:args.top]:
        star = "*" if f["is_named_active_site"] else " "
        add("  %s%-12s %-22s %-34s %s" % (
            star, f["pocket_residue"], f["ligand_feature_type"],
            f["suggested_ligand_chemistry"][:34],
            f["approx_xyz"]))
    if len(features) > args.top:
        add("   ... %d more features in the JSON file." % (len(features) - args.top))
    add("   (* = residue you named as active site; these are weighted highest.)")
    add("")
    add("6. SUGGESTED LIGAND PROPERTY WINDOWS")
    add("   Molecular weight (Da): %s - %s" % tuple(properties["molecular_weight_Da"]))
    add("   cLogP                : %s - %s" % tuple(properties["clogp"]))
    add("   Rotatable bonds      : <= %d" % properties["rotatable_bonds_max"])
    add("   H-bond donors        : %s - %s" % tuple(properties["hbond_donors"]))
    add("   H-bond acceptors     : %s - %s" % tuple(properties["hbond_acceptors"]))
    add("   Net charge           : %+d" % properties["suggested_net_charge"])
    add("   Basis                : %s" % properties["basis"])
    add("")
    add("7. HOW TO READ THIS")
    add("   Each line in section 5 says: 'the pocket offers X here, so your")
    add("   ligand should present the opposite of X at roughly these")
    add("   coordinates'. Build or search for a molecule that places those")
    add("   groups near those points, then dock and test it.")
    add("")
    add(DISCLAIMER)
    return "\n".join(L)


# ---------------------------------------------------------------------------
# SECTION 8: Command line
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Profile a protein pocket and propose a complementary "
                    "(repressor) ligand pharmacophore.")
    ap.add_argument("structure",
                    help="4-character PDB code (downloaded from RCSB) or a local .pdb/.cif file")
    ap.add_argument("residues",
                    help="Active-site residues, e.g. 'A:45,A:88,A:102' or '45,88,102' with --chain A")
    ap.add_argument("--chain", default=None,
                    help="Default chain for residue numbers given without one")
    ap.add_argument("--cutoff", type=float, default=8.0,
                    help="Distance in angstroms defining the pocket shell (default 8)")
    ap.add_argument("--top", type=int, default=25,
                    help="How many ranked features to print (default 25)")
    ap.add_argument("--outdir", default="output", help="Output folder (default: output)")
    args = ap.parse_args()

    try:
        wanted = parse_residue_spec(args.residues, args.chain)
    except ValueError as exc:
        raise SystemExit("ERROR: %s" % exc)

    structure, name = load_structure(args.structure, os.path.join(args.outdir, "downloads"))
    model = first_model(structure)

    active_site = collect_active_site(model, wanted)
    shell = collect_pocket_shell(model, active_site, args.cutoff)
    if not shell:
        raise SystemExit("ERROR: no residues found within the cutoff. Try a larger --cutoff.")

    profile = profile_pocket(shell)
    geometry = pocket_geometry(shell)
    sasa = solvent_accessibility(model, shell)
    features = build_pharmacophore(active_site, shell)
    properties = suggest_properties(profile, geometry)

    report = render_report(name, args, active_site, shell, profile, geometry,
                           sasa, features, properties)
    print(report)

    os.makedirs(args.outdir, exist_ok=True)
    txt_path = os.path.join(args.outdir, "%s_pocket_report.txt" % name)
    json_path = os.path.join(args.outdir, "%s_pharmacophore.json" % name)
    with open(txt_path, "w") as fh:
        fh.write(report + "\n")
    with open(json_path, "w") as fh:
        json.dump({
            "structure": name,
            "cutoff_A": args.cutoff,
            "active_site_residues": ["%s%s" % (one_letter(r), residue_key(c, r))
                                     for c, r in active_site],
            "pocket_shell_residues": sorted(shell.keys()),
            "profile": profile,
            "geometry": geometry,
            "solvent_accessibility": sasa,
            "pharmacophore_features": features,
            "suggested_ligand_properties": properties,
            "disclaimer": DISCLAIMER.replace("\n", " "),
        }, fh, indent=2)

    print("\nWrote %s" % txt_path)
    print("Wrote %s" % json_path)


if __name__ == "__main__":
    main()
