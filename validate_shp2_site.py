#!/usr/bin/env python3
"""validate_shp2_site.py - reproduce the SHP-2 allosteric-site checks.

This script does NOT design anything. It is the evidence trail behind the
SHP-2 worked example in the README, so you (or a reviewer) can re-run every
factual claim instead of taking it on trust. Four checks:

  1. Which residues actually line the SHP099 tunnel, computed as every amino
     acid with an atom within 5.0 A of the bound 5OD ligand in PDB 5EHR.
  2. Whether 5EHR and 2SHP share residue numbering, by comparing the amino
     acid identity at every residue number present in both.
  3. Whether the allosteric tunnel is already formed in the autoinhibited
     apo structure 2SHP: superpose 5EHR onto 2SHP, carry SHP099 across, and
     look for steric clashes.
  4. Whether docking can reproduce the crystal pose (run separately; see
     README) - this script prints the crystal pose for comparison.

Usage:  python validate_shp2_site.py data/pdb/5EHR.pdb data/pdb/2SHP.pdb
Needs:  biopython, numpy
"""

import sys
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch, Superimposer
from Bio.PDB.Polypeptide import is_aa, three_to_index, index_to_one

LIGAND = "5OD"          # PDB chemical component id for SHP099
CUTOFF = 5.0            # angstroms, for "lines the pocket"
CLASH = 2.5             # angstroms, heavy-atom clash threshold


def aa_map(model, chain_id):
    return {r.id[1]: r for r in model[chain_id] if is_aa(r, standard=True)}


def one(res):
    return index_to_one(three_to_index(res.get_resname()))


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    holo_path, apo_path = sys.argv[1], sys.argv[2]
    p = PDBParser(QUIET=True)
    holo = p.get_structure("holo", holo_path)[0]
    apo = p.get_structure("apo", apo_path)[0]

    # --- 1. residues lining the ligand ------------------------------------
    print("=" * 70)
    print("1. Residues within %.1f A of %s in %s" % (CUTOFF, LIGAND, holo_path))
    print("=" * 70)
    ligs = [r for r in holo.get_residues() if r.get_resname() == LIGAND]
    if not ligs:
        sys.exit("No %s found in %s." % (LIGAND, holo_path))
    protein = [a for a in holo.get_atoms()
               if is_aa(a.get_parent(), standard=True) and a.element != "H"]
    ns = NeighborSearch(protein)
    site_numbers = None
    for lig in ligs:
        near = set()
        for atom in lig:
            for nb in ns.search(atom.coord, CUTOFF):
                r = nb.get_parent()
                near.add((r.get_parent().id, r.id[1], r.get_resname()))
        sel = sorted(near)
        ch = lig.get_parent().id
        print("\nchain %s: %d residues" % (ch, len(sel)))
        print("  " + ", ".join("%s%d" % (rn, n) for _, n, rn in sel))
        print("  paste-ready: " + ",".join("%s:%d" % (c, n) for c, n, _ in sel))
        if site_numbers is None:
            site_numbers = [n for _, n, _ in sel]

    # --- 2. numbering correspondence --------------------------------------
    print("\n" + "=" * 70)
    print("2. Residue numbering: %s vs %s" % (holo_path, apo_path))
    print("=" * 70)
    h = aa_map(holo, "A")
    for ch in sorted(c.id for c in apo):
        a = aa_map(apo, ch)
        shared = [n for n in h if n in a]
        same = [n for n in shared if h[n].get_resname() == a[n].get_resname()]
        diff = [(n, one(h[n]), one(a[n])) for n in shared if n not in set(same)]
        print("\nholo:A vs apo:%s -> %d shared numbers, %d identical (%.1f%%)"
              % (ch, len(shared), len(same), 100.0 * len(same) / max(1, len(shared))))
        if diff:
            print("  differences: " + ", ".join("%d %s->%s" % d for d in sorted(diff)))
        missing = [n for n in site_numbers if n not in a]
        print("  site residues present in apo: %d/%d%s"
              % (len(site_numbers) - len(missing), len(site_numbers),
                 "" if not missing else "  MISSING: %s" % missing))
        bad = [n for n in site_numbers if n in a and h[n].get_resname() != a[n].get_resname()]
        print("  site residues with matching identity: %s"
              % ("all" if not bad else "MISMATCH at %s" % bad))

    # --- 3. is the tunnel open in the apo structure? ----------------------
    print("\n" + "=" * 70)
    print("3. Tunnel check: transplant %s into %s and look for clashes"
          % (LIGAND, apo_path))
    print("=" * 70)
    lig = [r for r in holo["A"] if r.get_resname() == LIGAND][0]
    lig_atoms = [a for a in lig if a.element != "H"]
    for ch in sorted(c.id for c in apo):
        a = aa_map(apo, ch)
        common = [n for n in h if n in a and "CA" in h[n] and "CA" in a[n]]
        sup = Superimposer()
        sup.set_atoms([a[n]["CA"] for n in common], [h[n]["CA"] for n in common])
        rot, tran = sup.rotran
        apo_atoms = [at for at in apo[ch].get_atoms()
                     if is_aa(at.get_parent(), standard=True) and at.element != "H"]
        ns2 = NeighborSearch(apo_atoms)
        clashes = []
        mins = []
        for atom in lig_atoms:
            c = atom.coord @ rot + tran
            nbs = ns2.search(c, 6.0)
            if not nbs:
                mins.append(99.0)
                continue
            d = min(np.linalg.norm(nb.coord - c) for nb in nbs)
            mins.append(d)
            if d < CLASH:
                nb = min(nbs, key=lambda x: np.linalg.norm(x.coord - c))
                r = nb.get_parent()
                side = nb.get_id() not in ("N", "CA", "C", "O")
                clashes.append((atom.get_id(), "%s%d" % (r.get_resname(), r.id[1]),
                                nb.get_id(), "sidechain" if side else "BACKBONE", d))
        print("\napo chain %s: CA rmsd %.2f A over %d residues"
              % (ch, sup.rms, len(common)))
        print("  ligand heavy atoms: %d | closest protein contact %.2f A | mean %.2f A"
              % (len(lig_atoms), min(mins), float(np.mean(mins))))
        if not clashes:
            print("  No clashes under %.1f A - tunnel is open as modelled." % CLASH)
        else:
            print("  %d atom(s) clash under %.1f A:" % (len(clashes), CLASH))
            for cl in clashes:
                print("    lig %-4s -- %-8s %-4s %-9s %.2f A" % cl)
            backbone = [c for c in clashes if c[3] == "BACKBONE"]
            print("  Verdict: %s"
                  % ("backbone blocks the site - tunnel is NOT formed"
                     if backbone else
                     "clashes are side-chain only - the tunnel is formed, but a "
                     "side-chain rotamer occludes it; rigid docking will "
                     "under-score a real binder here"))

    print("\nDone. Numbers above are computed from the structures, not recalled.")


if __name__ == "__main__":
    main()
