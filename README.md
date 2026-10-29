# binder-app

Two small, beginner-friendly Python tools:

1. **`pocket_profiler.py`** – reads a protein structure, looks at an active site, and works backwards to a blueprint ("pharmacophore") for a ligand that should bind there and repress the protein.
2. **`plasmid_builder.py`** – takes protein sequences in a FASTA file and designs codon-optimised, His-tagged E. coli genes ready for NdeI/XhoI cloning.

## Setup

```bash
cd binder-app
python3 -m venv venv
./venv/bin/pip install biopython dnachisel
```

Nothing else is needed. (Installed and tested with Biopython 1.88 and DNA Chisel 3.2.16.)

## Running them

```bash
./venv/bin/python pocket_profiler.py 3ERT A:351,A:383,A:524 --cutoff 6
./venv/bin/python plasmid_builder.py test_proteins.fasta
```

Everything is written into `output/`.

---

## pocket_profiler.py, section by section

**Section 1 – Chemistry lookup tables.** Plain dictionaries of facts about the 20 amino acids: three-letter to one-letter codes, Kyte-Doolittle hydropathy (how water-hating each one is), which are hydrophobic / aromatic / positive / negative / polar, how many hydrogen-bond donors and acceptors each side chain offers, and which atoms matter chemically in each side chain.

**Section 2 – Small helpers.** Short utilities: get a residue's one-letter code, grab only its side-chain atoms (everything that isn't backbone N/CA/C/O), make a readable label like `A:351`, and read your residue list from the command line (`A:45,A:88` or `45,88` with `--chain A`).

**Section 3 – Loading the structure.** If you pass a filename it parses it (`.cif` with MMCIFParser, otherwise PDBParser). If you pass a 4-character code like `3ERT` it downloads it from RCSB with Biopython's `PDBList`. Bad codes and unreadable files produce a plain-English error instead of a crash.

**Section 4 – Selecting the site.** `collect_active_site` finds exactly the residues you named and complains clearly if a chain or residue doesn't exist. `collect_pocket_shell` then sweeps the whole structure and keeps every standard amino acid that has any atom within the cutoff (8 Å by default) of those side chains. That set is the wall of the pocket. Water, ions and bound ligands are skipped.

**Section 5 – Profiling.** Counts what the pocket is made of (composition, % hydrophobic/aromatic/charged/polar), averages the hydropathy, adds up formal charges to get the net charge at pH 7, totals the H-bond donors and acceptors, and measures geometry: centroid, bounding box and the longest distance between any two pocket atoms. Solvent accessibility is computed with Shrake-Rupley **only if** the installed Biopython has it – the code imports it inside a `try` and reports honestly if it is missing.

**Section 6 – Working backwards.** This is the heart of it. `COMPLEMENT_RULES` maps each pocket residue to the opposite thing a ligand should present: Arg/Lys → an anionic group (carboxylate/phosphate), Asp/Glu → a cationic group (amine/amidine), Phe/Trp/Tyr → an aryl ring for π-stacking, Ser/Thr → an H-bond acceptor, Leu/Ile/Val/Met → a greasy alkyl contact, Gly → open space. Each feature gets an approximate 3D coordinate (the average position of the chemically important atoms of that side chain) and a priority score; residues you named as the active site are weighted 1.5× so they rank first. `suggest_properties` then turns pocket size and chemistry into rough ligand property windows – molecular weight scaled from the pocket span, logP shifted by the mean hydropathy, and a suggested net charge opposite the pocket's.

**Section 7 – Reporting.** Formats everything into a readable seven-part report, ending with an explicit disclaimer that this is a hypothesis, not a docked or validated binder.

**Section 8 – Command line.** `argparse` wiring, then writes `output/<name>_pocket_report.txt` and `output/<name>_pharmacophore.json`.

### Honest limits
No docking, no energy function, no conformational sampling, no estimate of affinity. Coordinates are side-chain anchor points, not computed ligand atom positions. The property windows are rules of thumb. It tells you what to look for; it cannot tell you that anything will bind.

---

## plasmid_builder.py, section by section

**Section 1 – Constants.** The 20 standard amino acids, the `LEHHHHHH` tag, the NdeI (`CATATG`) and XhoI (`CTCGAG`) sites, the 20–300 aa length limits, and the longest allowed run of a single base.

**Section 2 – Validation and tagging.** Each sequence is rejected (with a clear message, then skipped) if it contains a stop codon `*`, any non-standard letter, or is outside 20–300 aa. Valid ones get a starting Met if they don't already have one, and the `LEHHHHHH` tag appended.

**Section 3 – Codon optimisation.** `reverse_translate` makes a first rough DNA sequence. That goes into a DNA Chisel `DnaOptimizationProblem` with hard constraints: `EnforceTranslation` (the DNA must still code for your exact protein), `AvoidPattern` for NdeI and XhoI so no extra sites appear inside the gene, `EnforceGCContent(mini=0.40, maxi=0.60, window=50)`, and `AvoidPattern("7xA")` etc. to kill long homopolymer runs. The soft objective is `CodonOptimize(species="e_coli")`. `resolve_constraints()` makes the rules true, `optimize()` then improves codon usage. `add_flanks` prepends `CAT` so that `CAT + ATG…` spells `CATATG` with the ATG being the start codon, and appends `TAA` + `CTCGAG` at the end, keeping the reading frame a multiple of 3.

**Section 4 – Report and checks.** Slices the CDS back out, translates it with Biopython, and confirms it matches the tagged protein character for character. Prints length, GC%, frame check, and the counts of NdeI and XhoI sites (should be exactly one each).

**Section 5 – Writing files.** Writes `output/<name>.dna.fasta`, and a GenBank file `output/<name>.gb` annotated with the NdeI site, the CDS (with translation), the His-tag region, the TAA stop and the XhoI site.

**Section 6 – Main.** Reads the FASTA, loops over every record, and prints a summary of how many succeeded.

### Honest limits
The output is the **insert only** – no plasmid backbone is invented anywhere. You paste this into a real NdeI/XhoI vector you already have (e.g. a pET plasmid). The DNA Chisel codon table for E. coli comes bundled with the library.
