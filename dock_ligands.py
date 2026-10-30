#!/usr/bin/env python3
"""
dock_ligands.py - STEP 3 of the alpha-ligand pipeline.

The three steps, in order:

    1. pocket_profiler.py   protein  -> output/<CODE>_pharmacophore.json
                            (what chemistry does this pocket ask for?)
    2. ligand_builder.py    that JSON -> output/<CODE>_ligands.sdf
                            (real 3D molecules whose chemistry matches)
    3. dock_ligands.py      <- YOU ARE HERE
                            molecules + protein -> output/<CODE>_docked.sdf
                            (does each molecule physically FIT the pocket,
                             and how tightly, according to AutoDock Vina?)

Step 2 only ever asked "does this molecule have the right chemical groups?".
It never checked whether the molecule can physically get into the pocket.
This script answers that: it puts each molecule inside the pocket in many
orientations and shapes, and reports the best score Vina can find.

PLEASE READ THE HONESTY NOTE AT THE BOTTOM OF THE REPORT. A good docking
score is NOT proof that a molecule binds.

What you need installed (this is the Python step - it cannot run in the
browser web app):
    * RDKit          (reads/writes molecules)
    * Biopython      (reads the protein structure)
    * meeko          (converts molecules and proteins to Vina's PDBQT format)
    * the `vina` program itself, on your PATH

Example use:

    python dock_ligands.py output/3ERT_pharmacophore.json \
        --ligands output/3ERT_ligands.sdf \
        --receptor 3ERT

Author's note: every section below is commented for beginners.
"""

# ---------------------------------------------------------------------------
# SECTION 1: Imports
# ---------------------------------------------------------------------------

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

# RDKit: the chemistry toolkit. We use it to read the SDF of candidate
# molecules and to write the docked poses back out as an SDF.
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")          # keep the console tidy

# Biopython bits for writing a cleaned-up copy of the protein.
from Bio.PDB import PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

# meeko: turns an RDKit molecule into the PDBQT text Vina wants.
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate

# We reuse the structure-loading code from step 1 so that a 4-character PDB
# code like "3ERT" works here exactly as it does there.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pocket_profiler import load_structure, first_model


# ---------------------------------------------------------------------------
# SECTION 2: Finding the external programs we depend on
# ---------------------------------------------------------------------------

def find_vina():
    """Locate the `vina` executable, or stop with a helpful message."""
    path = shutil.which("vina")
    if path:
        return path
    sys.exit(
        "ERROR: could not find the 'vina' program on your PATH.\n"
        "       Download it from:\n"
        "         https://github.com/ccsb-scripps/AutoDock-Vina/releases\n"
        "       On Windows, put vina.exe in this folder or add it to PATH.")


def find_receptor_prep():
    """Locate meeko's mk_prepare_receptor.py script that ships with meeko."""
    # meeko installs it next to the python interpreter (Scripts/ on Windows).
    candidates = [
        os.path.join(os.path.dirname(sys.executable), "mk_prepare_receptor.py"),
        os.path.join(os.path.dirname(sys.executable), "mk_prepare_receptor"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    found = shutil.which("mk_prepare_receptor.py") or shutil.which("mk_prepare_receptor")
    if found:
        return found
    sys.exit("ERROR: meeko's mk_prepare_receptor.py was not found. "
             "Install meeko with:  pip install meeko")


# ---------------------------------------------------------------------------
# SECTION 3: Reading the blueprint from step 1
# ---------------------------------------------------------------------------

def load_pharmacophore(json_path):
    """Read the JSON written by pocket_profiler.py."""
    if not os.path.isfile(json_path):
        sys.exit("ERROR: '%s' does not exist. Run pocket_profiler.py first."
                 % json_path)
    with open(json_path) as fh:
        try:
            return json.load(fh)
        except Exception as exc:
            sys.exit("ERROR: '%s' is not valid JSON: %s" % (json_path, exc))


def docking_box(data, padding, max_side):
    """Work out the search box Vina should look inside.

    pocket_profiler.py already measured the pocket's centre point (centroid)
    and how wide it is in x, y and z (the bounding box). We start from that,
    add a little padding so a ligand can stick out slightly, and then cap each
    side: a huge box makes docking slow and sloppy, because Vina then has to
    search a lot of empty space.
    """
    geom = data.get("geometry", {})
    centre = geom.get("centroid")
    box = geom.get("bounding_box_xyz")
    if not centre or not box:
        sys.exit("ERROR: the pharmacophore JSON has no geometry/centroid. "
                 "Re-run pocket_profiler.py to regenerate it.")
    size = [min(float(b) + 2.0 * padding, float(max_side)) for b in box]
    # Vina also refuses boxes that are too small to hold anything.
    size = [max(s, 16.0) for s in size]
    return [float(c) for c in centre], size


# ---------------------------------------------------------------------------
# SECTION 4: Preparing the protein (the "receptor")
# ---------------------------------------------------------------------------

class ProteinOnly(Select):
    """A filter telling Biopython which atoms to keep when saving the protein.

    Docking wants the protein and nothing else:
      * waters are dropped - they are usually not part of the pocket we model
      * any other HETATM (ions, sugars, and crucially the ORIGINAL LIGAND)
        is dropped, otherwise the pocket would already be occupied
      * hydrogens are dropped; meeko adds its own, consistently
      * only one chain is kept, the one your active-site residues live on
      * for atoms modelled in two positions (altloc) we keep the first
    """

    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return 1 if chain.get_id() == self.chain_id else 0

    def accept_residue(self, residue):
        # is_aa(standard=True) keeps the 20 normal amino acids only, which
        # automatically throws away water (HOH) and all other heteroatoms.
        return 1 if is_aa(residue, standard=True) else 0

    def accept_atom(self, atom):
        if atom.element == "H":
            return 0
        if atom.get_altloc() not in (" ", "A"):
            return 0
        return 1


def guess_chain(data):
    """Which chain did the user profile in step 1? Read it off the JSON."""
    for spec in data.get("active_site_residues", []):
        # entries look like "DA:351" (one-letter aa + chain + ':' + number)
        if ":" in spec:
            left = spec.split(":")[0]
            if left:
                return left[-1]
    for spec in data.get("pocket_shell_residues", []):
        if ":" in spec:
            return spec.split(":")[0][-1]
    return "A"


def prepare_receptor(receptor_source, chain_id, centre, size, workdir, verbose):
    """Turn a protein structure into the PDBQT file Vina reads.

    Two stages:
      1. Biopython writes a cleaned PDB (protein only, see ProteinOnly above).
      2. meeko's mk_prepare_receptor.py adds hydrogens and partial charges and
         writes PDBQT.
    """
    print("Loading receptor structure: %s" % receptor_source)
    structure, name = load_structure(receptor_source, os.path.join(workdir, "downloads"))
    model = first_model(structure)

    chains = [c.get_id() for c in model]
    if chain_id not in chains:
        sys.exit("ERROR: chain '%s' is not in this structure. Chains present: %s"
                 % (chain_id, ", ".join(chains)))

    clean_pdb = os.path.join(workdir, "receptor_clean.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(clean_pdb, select=ProteinOnly(chain_id))

    n_kept = sum(1 for r in model[chain_id] if is_aa(r, standard=True))
    print("  Kept chain %s: %d standard amino acids "
          "(waters, ions and any bound ligand removed)." % (chain_id, n_kept))

    basename = os.path.join(workdir, "receptor")
    cmd = [
        sys.executable, find_receptor_prep(),
        "--read_pdb", clean_pdb,
        "-o", basename,
        "-p",                       # write the PDBQT
        "--box_center", "%.3f" % centre[0], "%.3f" % centre[1], "%.3f" % centre[2],
        "--box_size", "%.1f" % size[0], "%.1f" % size[1], "%.1f" % size[2],
        "--default_altloc", "A",
        "-x",                       # allow_bad_res: skip broken residues
    ]
    print("  Converting to PDBQT with meeko ...")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    pdbqt = basename + ".pdbqt"
    if not os.path.isfile(pdbqt):
        print((proc.stdout or "")[-3000:])
        print((proc.stderr or "")[-3000:])
        sys.exit("ERROR: receptor preparation failed (no PDBQT was written).")

    # meeko tells us which residues it had to drop; pass that on honestly.
    combined_log = (proc.stdout or "") + (proc.stderr or "")
    dropped = re.findall(r"Template matching failed for: \[(.*?)\]", combined_log)
    if dropped:
        names = [d.strip().strip("'") for d in dropped[-1].split(",")]
        print("  NOTE: meeko skipped %d residue(s) it could not model "
              "(usually incomplete side chains in the crystal structure):"
              % len(names))
        print("        %s" % ", ".join(names))
    if verbose:
        print(combined_log)
    print("  Receptor ready: %s" % pdbqt)
    return pdbqt, name, dropped[-1] if dropped else ""


# ---------------------------------------------------------------------------
# SECTION 5: Preparing each ligand
# ---------------------------------------------------------------------------

def read_candidates(sdf_path):
    """Read the SDF from step 2 and keep ONE 3D structure per molecule.

    The builder writes several conformers (shapes) of each molecule. Docking
    does its own conformational search - it rotates every rotatable bond
    itself - so handing it ten shapes of the same molecule would just repeat
    the same work. We therefore keep the first (lowest-strain) conformer of
    each distinct molecule and let Vina reshape it.
    """
    if not os.path.isfile(sdf_path):
        sys.exit("ERROR: '%s' does not exist. Run ligand_builder.py first."
                 % sdf_path)
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    seen = {}
    order = []
    for mol in supplier:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        name = props.get("molecule_name") or mol.GetProp("_Name")
        if name in seen:
            continue
        seen[name] = mol
        order.append(name)
    if not order:
        sys.exit("ERROR: no readable molecules in '%s'." % sdf_path)
    return [(n, seen[n]) for n in order]


def molecule_from_smiles(smiles, name, seed):
    """Build a 3D molecule from a SMILES string (used for the control ligand)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.exit("ERROR: '%s' is not a valid SMILES string." % smiles)
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        sys.exit("ERROR: RDKit could not generate 3D coordinates for %s." % name)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    mol.SetProp("_Name", name)
    mol.SetProp("molecule_name", name)
    mol.SetProp("smiles", smiles)
    return mol


def ligand_to_pdbqt(mol, path):
    """Write one molecule as a PDBQT file. Returns True on success."""
    work = Chem.Mol(mol)
    if work.GetNumConformers() == 0:
        return False, "molecule has no 3D coordinates"
    # Vina needs explicit hydrogens to assign atom types correctly.
    try:
        work = Chem.AddHs(work, addCoords=True)
    except Exception as exc:
        return False, "AddHs failed: %s" % exc
    try:
        setups = MoleculePreparation().prepare(work)
    except Exception as exc:
        return False, "meeko preparation failed: %s" % exc
    if not setups:
        return False, "meeko returned no setup"
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return False, "PDBQT writer refused it: %s" % err
    with open(path, "w") as fh:
        fh.write(text)
    return True, ""


# ---------------------------------------------------------------------------
# SECTION 6: Running Vina
# ---------------------------------------------------------------------------

SCORE_RE = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+\.\d+)", re.MULTILINE)


def run_vina(vina_exe, receptor, ligand, out_path, centre, size,
             exhaustiveness, n_poses, seed, cpu):
    """Call the vina program once, for one ligand. Returns (scores, log)."""
    cmd = [
        vina_exe,
        "--receptor", receptor,
        "--ligand", ligand,
        "--out", out_path,
        "--center_x", "%.3f" % centre[0],
        "--center_y", "%.3f" % centre[1],
        "--center_z", "%.3f" % centre[2],
        "--size_x", "%.1f" % size[0],
        "--size_y", "%.1f" % size[1],
        "--size_z", "%.1f" % size[2],
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--seed", str(seed),
        "--cpu", str(cpu),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log = (proc.stdout or "") + (proc.stderr or "")
    if not os.path.isfile(out_path):
        return [], log
    with open(out_path) as fh:
        text = fh.read()
    scores = [float(s) for s in SCORE_RE.findall(text)]
    return scores, log


def poses_to_rdkit(out_pdbqt):
    """Turn Vina's PDBQT output back into one RDKit molecule PER POSE.

    meeko hands us the molecule once, with each docked pose stored as a
    separate 'conformer' (a set of coordinates) inside it. SDF readers are
    much easier to use if each pose is its own entry, so we split them here.
    """
    poses = []
    try:
        pdbqt_mol = PDBQTMolecule.from_file(out_pdbqt, skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    except Exception:
        return []
    for mol in mols:
        if mol is None:
            continue
        for conf in mol.GetConformers():
            single = Chem.Mol(mol)
            single.RemoveAllConformers()
            single.AddConformer(Chem.Conformer(conf), assignId=True)
            poses.append(single)
    return poses


# ---------------------------------------------------------------------------
# SECTION 7: The honesty note and the report
# ---------------------------------------------------------------------------

HONESTY = """
HOW MUCH SHOULD YOU TRUST THESE NUMBERS?
----------------------------------------
Read this before you spend any money on the molecules above.

* A Vina score is an EMPIRICAL ESTIMATE, not a measurement. It is a formula
  fitted to a set of known protein-ligand complexes. It is reported in
  kcal/mol and more negative means "Vina likes it more".

* Docking scores correlate only LOOSELY with real binding affinity. In
  published benchmarks the typical error is on the order of 2-3 kcal/mol,
  which is a 30-100x error in affinity terms. Vina is reasonably good at
  telling "roughly fits" from "definitely does not fit", and much worse at
  ranking two molecules that both score well.

* A good score is NOT a validated binder. Molecules that dock beautifully
  frequently do not bind at all in the test tube. The only way to know is to
  obtain the compound and measure it.

* No Kd, Ki or IC50 is quoted anywhere in this report, because this
  calculation cannot produce one. Do not convert these scores into one.

* The protein was treated as RIGID. Real side chains move to accommodate a
  ligand, so a molecule scored as "too bulky" might still fit in reality,
  and vice versa.

* Water was removed. Some pockets rely on bridging water molecules, and
  removing them changes the answer.

* Only one protonation state was used, the default one meeko assigns at
  roughly pH 7. Charged groups in particular may behave differently.

* Nothing here says anything about whether a molecule is stable, soluble,
  selective, non-toxic, cell-permeable, synthesisable, or purchasable.

WHAT THIS LIST ACTUALLY IS: a shortlist, ranked by one cheap computational
filter, of molecules that are worth looking at next. Treat it as a starting
point for experiments, not as a result.
"""


def write_report(path, code, args, centre, size, results, control,
                 dropped_residues, n_input_molecules):
    """Write the plain-English ranking report."""
    lines = []
    add = lines.append
    bar = "=" * 72

    add(bar)
    add("DOCKING REPORT - %s" % code)
    add(bar)
    add("")
    add("Produced by dock_ligands.py (step 3 of the alpha-ligand pipeline).")
    add("Docking engine: AutoDock Vina. Scores are in kcal/mol; more negative")
    add("means Vina thinks the fit is better.")
    add("")

    add("-" * 72)
    add("1. WHAT WAS DOCKED, AND WHERE")
    add("-" * 72)
    add("Receptor               : %s" % args.receptor)
    add("Ligand file            : %s" % args.ligands)
    add("Pharmacophore blueprint: %s" % args.pharmacophore_json)
    add("Distinct molecules     : %d (one 3D start shape each, from %d SDF entries)"
        % (len(results), n_input_molecules))
    add("Search box centre      : x %.2f  y %.2f  z %.2f" % tuple(centre))
    add("Search box size (Angst): x %.1f  y %.1f  z %.1f" % tuple(size))
    add("   The centre is the pocket centroid measured in step 1. The size is")
    add("   the pocket's bounding box plus %.1f A padding on each side, capped" % args.box_padding)
    add("   at %.1f A so the search does not waste time on empty space." % args.max_box)
    add("Exhaustiveness         : %d  (higher = slower, more thorough search)" % args.exhaustiveness)
    add("Poses kept per molecule: %d" % args.poses)
    add("Random seed            : %d  (same seed = same numbers next time)" % args.seed)
    if dropped_residues:
        add("")
        add("Residues meeko could not model and therefore left out of the")
        add("receptor (usually incomplete side chains in the crystal file):")
        add("   %s" % dropped_residues)
    add("")

    if control:
        add("-" * 72)
        add("2. CONTROL / SANITY CHECK")
        add("-" * 72)
        add("Control molecule : %s" % control["name"])
        add("SMILES           : %s" % control["smiles"])
        if control.get("best") is None:
            add("Result           : DOCKING FAILED - see console output.")
        else:
            add("Best Vina score  : %+.2f kcal/mol" % control["best"])
            add("All poses        : %s"
                % ", ".join("%+.2f" % s for s in control["scores"]))
        add("")
        add("Why this matters: if a molecule already KNOWN to bind this pocket")
        add("scores poorly here, the setup (box, receptor, protonation) is")
        add("probably wrong and every other number below is suspect. If it")
        add("scores well, the setup is at least self-consistent - that still")
        add("does not make the generated molecules below real binders.")
        add("")

    add("-" * 72)
    add("3. RANKING BY BEST DOCKING SCORE")
    add("-" * 72)
    add("")
    add("%-4s %-42s %10s %10s %7s" % ("rank", "molecule", "best kcal", "checklist", "MW"))
    add("%-4s %-42s %10s %10s %7s" % ("", "", "/mol", "0-100", ""))
    add("-" * 78)
    for i, r in enumerate(results, 1):
        best = "%+.2f" % r["best"] if r["best"] is not None else "failed"
        chk = "%.1f" % r["checklist"] if r["checklist"] is not None else "-"
        mw = "%.0f" % r["mw"] if r["mw"] is not None else "-"
        name = r["name"]
        add("%-4d %-42s %10s %10s %7s" % (i, name[:42], best, chk, mw))
    add("")
    add("'checklist' is the step-2 score: how well the molecule's chemical")
    add("groups match what the pocket asks for. It knows nothing about shape.")
    add("'best kcal/mol' is this step: can the molecule physically fit, and")
    add("how favourably. The two often DISAGREE, and that disagreement is the")
    add("whole point of doing this step.")
    add("")

    # How much did the ranking move?
    ranked_ok = [r for r in results if r["best"] is not None and r["checklist"] is not None]
    if len(ranked_ok) > 2:
        by_chk = sorted(ranked_ok, key=lambda r: -r["checklist"])
        chk_rank = {r["name"]: i for i, r in enumerate(by_chk, 1)}
        add("-" * 72)
        add("4. HOW THE RANKING CHANGED VERSUS THE CHEMISTRY-ONLY SCORE")
        add("-" * 72)
        add("%-42s %10s %10s %8s" % ("molecule", "dock rank", "chem rank", "moved"))
        add("-" * 74)
        for i, r in enumerate(ranked_ok, 1):
            old = chk_rank[r["name"]]
            move = old - i
            arrow = "same" if move == 0 else ("up %d" % move if move > 0 else "down %d" % -move)
            add("%-42s %10d %10d %8s" % (r["name"][:42], i, old, arrow))
        add("")

    add(HONESTY)
    add(bar)
    add("END OF REPORT")
    add(bar)

    text = "\n".join(lines) + "\n"
    with open(path, "w") as fh:
        fh.write(text)
    return text


# ---------------------------------------------------------------------------
# SECTION 8: Main - glue it all together
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Dock the molecules built by ligand_builder.py into the "
                    "pocket profiled by pocket_profiler.py, using AutoDock "
                    "Vina. Writes scored poses (SDF) and a plain-English "
                    "ranking report.")
    ap.add_argument("pharmacophore_json",
                    help="Path to <CODE>_pharmacophore.json from pocket_profiler.py "
                         "(supplies the docking box centre and size)")
    ap.add_argument("--ligands", required=True,
                    help="Path to <CODE>_ligands.sdf from ligand_builder.py")
    ap.add_argument("--receptor", required=True,
                    help="The protein: a 4-character PDB code (e.g. 3ERT) or a "
                         "local .pdb/.cif file")
    ap.add_argument("--outdir", default="output",
                    help="Output folder (default: output)")
    ap.add_argument("--exhaustiveness", type=int, default=8,
                    help="How hard Vina searches, per molecule (default: 8)")
    ap.add_argument("--poses", type=int, default=5,
                    help="How many poses to keep per molecule (default: 5)")
    ap.add_argument("--box-padding", dest="box_padding", type=float, default=2.0,
                    help="Angstroms added to each side of the pocket bounding "
                         "box (default: 2.0)")
    ap.add_argument("--max-box", dest="max_box", type=float, default=26.0,
                    help="Largest allowed box side in Angstroms (default: 26)")
    ap.add_argument("--chain", default=None,
                    help="Protein chain to dock into (default: taken from the "
                         "pharmacophore JSON)")
    ap.add_argument("--control-smiles", dest="control_smiles", default=None,
                    help="Optional SMILES of a molecule already known to bind "
                         "this pocket, docked as a sanity check")
    ap.add_argument("--control-name", dest="control_name", default="control",
                    help="Name for the control molecule (default: control)")
    ap.add_argument("--seed", type=int, default=42,
                    help="Random seed, so runs are repeatable (default: 42)")
    ap.add_argument("--cpu", type=int, default=0,
                    help="CPU cores for Vina; 0 = let Vina decide (default: 0)")
    ap.add_argument("--keep-temp", dest="keep_temp", action="store_true",
                    help="Keep the intermediate PDBQT files for inspection")
    ap.add_argument("--verbose", action="store_true",
                    help="Print the full receptor-preparation log")
    args = ap.parse_args()

    print("=" * 72)
    print("DOCK LIGANDS (AutoDock Vina)")
    print("=" * 72)

    vina_exe = find_vina()
    ver = subprocess.run([vina_exe, "--version"], capture_output=True, text=True)
    print("Using docking engine: %s" % (ver.stdout.strip() or vina_exe))

    # --- step 1: blueprint and box -------------------------------------
    data = load_pharmacophore(args.pharmacophore_json)
    code = data.get("structure", "docked")
    centre, size = docking_box(data, args.box_padding, args.max_box)
    print("Pocket from %s: centre (%.2f, %.2f, %.2f), box %.1f x %.1f x %.1f A"
          % (code, centre[0], centre[1], centre[2], size[0], size[1], size[2]))

    chain_id = args.chain or guess_chain(data)

    os.makedirs(args.outdir, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix="dock_%s_" % code)

    try:
        # --- step 2: receptor ------------------------------------------
        receptor_pdbqt, _, dropped = prepare_receptor(
            args.receptor, chain_id, centre, size, workdir, args.verbose)

        # --- step 3: ligands -------------------------------------------
        n_entries = sum(1 for m in Chem.SDMolSupplier(args.ligands, removeHs=False) if m)
        candidates = read_candidates(args.ligands)
        print("Read %d distinct molecules from %s (%d SDF entries total)."
              % (len(candidates), args.ligands, n_entries))

        jobs = list(candidates)
        control_info = None
        if args.control_smiles:
            control_mol = molecule_from_smiles(
                args.control_smiles, args.control_name, args.seed)
            jobs.append((args.control_name, control_mol))
            print("Plus 1 control molecule: %s" % args.control_name)

        # --- step 4: dock each one -------------------------------------
        print("")
        print("Docking (exhaustiveness %d, %d poses each) ..." % (
            args.exhaustiveness, args.poses))
        results = []
        all_poses = []
        for i, (name, mol) in enumerate(jobs, 1):
            tag = "lig%02d" % i
            lig_pdbqt = os.path.join(workdir, tag + ".pdbqt")
            out_pdbqt = os.path.join(workdir, tag + "_out.pdbqt")

            ok, why = ligand_to_pdbqt(mol, lig_pdbqt)
            if not ok:
                print("  [%2d/%2d] %-44s SKIPPED (%s)" % (i, len(jobs), name[:44], why))
                scores = []
            else:
                scores, log = run_vina(vina_exe, receptor_pdbqt, lig_pdbqt,
                                       out_pdbqt, centre, size,
                                       args.exhaustiveness, args.poses,
                                       args.seed, args.cpu)
                if not scores:
                    print("  [%2d/%2d] %-44s FAILED" % (i, len(jobs), name[:44]))
                    if args.verbose:
                        print(log[-1500:])
                else:
                    print("  [%2d/%2d] %-44s best %+7.2f kcal/mol"
                          % (i, len(jobs), name[:44], scores[0]))

            props = mol.GetPropsAsDict()
            checklist = props.get("checklist_score_0_100")
            try:
                mw = Descriptors.MolWt(Chem.RemoveHs(Chem.Mol(mol)))
            except Exception:
                mw = None

            record = {
                "name": name,
                "smiles": props.get("smiles", Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)))),
                "scores": scores,
                "best": scores[0] if scores else None,
                "checklist": checklist,
                "mw": mw,
                "is_control": bool(args.control_smiles and name == args.control_name),
            }

            # Collect the actual 3D poses so we can write them to SDF.
            if scores:
                for rank, pose in enumerate(poses_to_rdkit(out_pdbqt), 1):
                    if pose is None:
                        continue
                    pose.SetProp("_Name", "%s_pose%d" % (name, rank))
                    pose.SetProp("molecule_name", name)
                    pose.SetProp("pose_rank", str(rank))
                    if rank <= len(scores):
                        pose.SetProp("vina_score_kcal_per_mol", "%.3f" % scores[rank - 1])
                    pose.SetProp("vina_best_score_kcal_per_mol", "%.3f" % scores[0])
                    if checklist is not None:
                        pose.SetProp("checklist_score_0_100", str(checklist))
                    pose.SetProp("receptor", str(args.receptor))
                    pose.SetProp("is_control", "yes" if record["is_control"] else "no")
                    pose.SetProp("WARNING",
                                 "Vina score is an empirical estimate, not a "
                                 "measured affinity. A good score is not a "
                                 "validated binder.")
                    all_poses.append((name, rank, pose))

            if record["is_control"]:
                control_info = record
            else:
                results.append(record)

        # --- step 5: rank and write ------------------------------------
        results.sort(key=lambda r: (r["best"] is None, r["best"] if r["best"] is not None else 0))

        sdf_path = os.path.join(args.outdir, "%s_docked.sdf" % code)
        writer = Chem.SDWriter(sdf_path)
        order = {r["name"]: i for i, r in enumerate(results)}
        if control_info:
            order[control_info["name"]] = -1
        for name, rank, pose in sorted(all_poses, key=lambda t: (order.get(t[0], 999), t[1])):
            writer.write(pose)
        writer.close()

        report_path = os.path.join(args.outdir, "%s_docking_report.txt" % code)
        control_payload = None
        if control_info:
            control_payload = {
                "name": control_info["name"],
                "smiles": control_info["smiles"],
                "best": control_info["best"],
                "scores": control_info["scores"],
            }
        write_report(report_path, code, args, centre, size, results,
                     control_payload, dropped, n_entries)

        print("")
        print("Wrote %d poses to %s" % (len(all_poses), sdf_path))
        print("Wrote report to      %s" % report_path)
        print("")
        print("Best molecules by Vina score:")
        for r in results[:5]:
            if r["best"] is None:
                continue
            print("   %+7.2f kcal/mol   %s" % (r["best"], r["name"]))
        if control_info and control_info["best"] is not None:
            print("Control (%s): %+.2f kcal/mol"
                  % (control_info["name"], control_info["best"]))
        print("")
        print("REMINDER: Vina scores are empirical estimates that correlate only")
        print("loosely with real affinity. A good score is not a validated")
        print("binder. This is a shortlist to test experimentally, nothing more.")

    finally:
        if args.keep_temp:
            print("Intermediate files kept in: %s" % workdir)
        else:
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
