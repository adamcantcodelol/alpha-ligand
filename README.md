# binder-app

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adamcantcodelol/alpha-ligand/blob/main/Alpha_Ligand_Pipeline.ipynb)

**Want the whole pipeline with nothing to install?** Click the badge above - it opens
`Alpha_Ligand_Pipeline.ipynb` in Google Colab and runs all three steps (profile the
pocket, build molecules, dock them) on Google's computers, in your browser. See
[Run it in Google Colab](#run-it-in-google-colab-no-install) below.

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

## ligand_builder.py - turning the blueprint into real 3D molecules

`pocket_profiler.py` tells you what chemistry a binder should have.
`ligand_builder.py` takes that answer and builds actual molecules with real
3D coordinates you can open in a viewer or hand to a docking program.

### Running it on Windows (PowerShell)

First profile a pocket (this writes the JSON the builder needs):

```
cd C:\Users\Adam\Downloads\binder-app\alpha-ligand
.\venv\Scripts\python pocket_profiler.py 3ERT A:351,A:383,A:524 --cutoff 6
```

Then build the molecules:

```
.\venv\Scripts\python ligand_builder.py output\3ERT_pharmacophore.json
```

That's it. Two files appear in `output\`:

* `3ERT_ligands.sdf` - the 3D molecules, several conformers each. Open it in
  PyMOL, ChimeraX or the free Avogadro.
* `3ERT_ligand_report.txt` - a plain-English write-up of what was built and why.

### Options

| option | what it does | default |
| --- | --- | --- |
| `--outdir` | where to write the two files | `output` |
| `--top` | how many pharmacophore features to use as the chemistry shopping list | `12` |
| `--molecules` | how many best-scoring molecules to build in 3D | `10` |
| `--confs` | 3D conformers (shapes) to generate per molecule | `10` |
| `--seed` | random seed, so runs are repeatable | `61453` |

Example with more of everything:

```
.\venv\Scripts\python ligand_builder.py output\3ERT_pharmacophore.json --molecules 20 --confs 25
```

You need RDKit in the venv for this script:

```
.\venv\Scripts\pip install rdkit
```

### How it works, in plain English

1. Reads the pharmacophore JSON and adds up, per chemistry class (cationic,
   anionic, aromatic, hydrophobic, H-bond donor/acceptor), how loudly the
   pocket asks for it. That is the "shopping list".
2. Builds a few hundred candidate molecules by gluing a small hand-written
   library of ring cores (benzene, indole, quinoline, piperidine, ...) to a
   small library of substituent groups (amines, carboxylates, tetrazoles,
   alcohols, greasy alkyls, ...), plus a handful of complete known fragments.
   RDKit's `Chem.molzip` does the gluing.
3. Asks RDKit, with SMARTS substructure patterns, which of those chemistry
   classes each candidate actually contains, and measures its properties
   (molecular weight, cLogP, H-bond donors/acceptors, rotatable bonds, TPSA).
4. Gives each candidate a **checklist score out of 100**: 70% for chemistry
   match against the shopping list, 30% for falling inside the profiler's
   suggested property windows.
5. Takes the best ones (at most two per ring core, so the list is varied),
   generates 3D conformers with RDKit's ETKDGv3 embedder, energy-minimises
   each one with the MMFF94 force field (UFF as fallback), and writes them
   to the SDF sorted so the least-strained conformer of each molecule is first.

### Honest limits - please read this bit

**The checklist score is not a binding affinity.** It is not a Kd, Ki, IC50 or
docking score. It only measures how many boxes on the chemical shopping list a
molecule ticks and whether it is roughly the right size. Two molecules with the
same score can bind a thousand-fold differently.

**Nothing is ever placed inside the protein.** The 3D coordinates are the
molecule's own relaxed shape on its own, not a pose in the pocket. The 3D
anchor points in the pharmacophore JSON are *not* used to position atoms - only
the list of feature types is used. A high-scoring molecule may be physically
unable to fit in the pocket at all.

**The force-field energies are internal strain** of the isolated molecule. They
are only useful for comparing different shapes of the *same* molecule, never
for comparing one molecule with another, and never as a binding number.

**The library is a starter kit**, not a real screening collection, and no
synthetic-accessibility check is done - some assembled molecules would be hard
or impossible to actually make.

The honest next step is docking, which is now implemented in
`dock_ligands.py` below. That is the first step that tests whether anything
physically fits.

---

## dock_ligands.py - does the molecule actually fit the pocket?

Steps 1 and 2 only ever asked *"does this molecule carry the right chemical
groups?"*. Neither one checked whether the molecule can physically get into the
pocket. `dock_ligands.py` answers that with **AutoDock Vina**: it places each
molecule inside the pocket in many orientations and shapes and reports the best
score it can find, in kcal/mol (more negative = Vina likes it more).

> **This step needs the Python install. It cannot run in the browser web app.**
> The web app in `docs/` stays the no-setup pocket profiler. RDKit, meeko and
> Vina are desktop programs - there is no way to run them inside a web page, so
> molecule building and docking are a local step.

### The full three-step pipeline on Windows (PowerShell)

Open PowerShell and run these one block at a time, from the repo folder.

**One-time setup**

```powershell
cd C:\Users\Adam\Downloads\binder-app\alpha-ligand
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install biopython dnachisel rdkit meeko scipy gemmi
```

Then download Vina itself (it is a single program, not a Python package) from
<https://github.com/ccsb-scripps/AutoDock-Vina/releases> - take the
`vina_1.2.5_win.exe` file, rename it to `vina.exe`, and put it in this repo
folder. Check it works:

```powershell
.\vina.exe --version
```

**Step 1 - profile the pocket**

```powershell
python pocket_profiler.py 3ERT A:351,A:383,A:524 --cutoff 6.0
```

Writes `output\3ERT_pharmacophore.json` and `output\3ERT_pocket_report.txt`.

**Step 2 - build candidate molecules**

```powershell
python ligand_builder.py output\3ERT_pharmacophore.json
```

Writes `output\3ERT_ligands.sdf` and `output\3ERT_ligand_report.txt`.

**Step 3 - dock them**

```powershell
python dock_ligands.py output\3ERT_pharmacophore.json --ligands output\3ERT_ligands.sdf --receptor 3ERT
```

Writes `output\3ERT_docked.sdf` (the poses, with scores attached) and
`output\3ERT_docking_report.txt` (the ranking in plain English).

If `vina.exe` is in the repo folder rather than on your PATH, add the folder to
PATH for that window first:

```powershell
$env:PATH = "$PWD;$env:PATH"
```

### Options

| flag | what it does |
|---|---|
| `--receptor` | 4-character PDB code (e.g. `3ERT`) **or** a local `.pdb`/`.cif` file |
| `--ligands` | the SDF from step 2 (required) |
| `--outdir` | where to write results (default `output`) |
| `--exhaustiveness` | how hard Vina searches per molecule (default 8; 16-32 is slower and more careful) |
| `--poses` | poses kept per molecule (default 5) |
| `--box-padding` | Angstroms added to each side of the pocket's bounding box (default 2.0) |
| `--max-box` | largest allowed box side (default 26) - stops the search wasting time on empty space |
| `--chain` | which protein chain to dock into (default: read from the JSON) |
| `--control-smiles` | SMILES of a molecule already known to bind, docked as a sanity check |
| `--seed` | random seed, so the same command gives the same numbers |
| `--keep-temp` | keep the intermediate PDBQT files |

### The sanity check, and why you should always use it

3ERT is the oestrogen receptor with **4-hydroxytamoxifen** bound, so we know
what a real binder for this pocket looks like. Docking it back in:

```powershell
python dock_ligands.py output\3ERT_pharmacophore.json --ligands output\3ERT_ligands.sdf --receptor 3ERT --control-smiles "CC/C(=C(\c1ccc(O)cc1)c1ccc(OCCN(C)C)cc1)c1ccccc1" --control-name "4-hydroxytamoxifen_(native)"
```

gives **-9.58 kcal/mol** for the native ligand, versus **-7.23 to -5.85** for
the ten generated molecules. That is the result you want from a control: the
known binder clearly beats the guesses, which means the box and receptor are
set up sensibly. If the control had scored *worse* than the generated
molecules, the setup would be wrong and every number would be suspect.

Note what it also tells you: our generated molecules are about 2.4 kcal/mol
behind the real drug. They are small fragments (MW 165-221 vs 388 for
tamoxifen), so they simply do not fill the pocket as well.

### How it works, in plain English

**Section 2 – finding the programs.** Looks for `vina` and for meeko's
`mk_prepare_receptor.py`, and gives a clear message with a download link if
either is missing.

**Section 3 – the docking box.** Vina only searches inside a box you give it.
The box centre is the pocket centroid that step 1 measured; the size is the
pocket bounding box plus padding, capped at `--max-box`.

**Section 4 – preparing the protein.** Biopython writes a cleaned copy:
one chain, the 20 standard amino acids only, which drops waters, ions and
crucially **the original bound ligand** (otherwise the pocket would already be
full). Then meeko adds hydrogens and partial charges and writes PDBQT.
Residues with incomplete side chains in the crystal file are skipped, and the
script prints and reports exactly which ones.

**Section 5 – preparing the ligands.** Reads the SDF and keeps one 3D start
shape per distinct molecule. Vina does its own conformational search, so
feeding it ten shapes of the same molecule would just repeat work.

**Section 6 – running Vina.** One call per molecule, reading back the
`REMARK VINA RESULT` scores and converting the posed PDBQT into RDKit
molecules, one entry per pose.

**Section 7 – the report.** Ranks molecules by their best score, shows the
step-2 chemistry score next to it, and shows how much the ranking *moved*
between the two.

### Honest limits - please read this bit

* A Vina score is an **empirical estimate**, not a measurement. Typical error
  in published benchmarks is 2-3 kcal/mol, which is a 30-100x error in
  affinity terms.
* **A good score is not a validated binder.** Molecules that dock beautifully
  often do not bind at all in the test tube.
* No Kd, Ki or IC50 appears anywhere, because this calculation cannot produce
  one. Do not convert the scores into one.
* The protein is treated as **rigid**. Real side chains move.
* Water is removed. Some pockets genuinely need bridging waters.
* One protonation state only, roughly pH 7.
* Nothing here says anything about stability, solubility, selectivity,
  toxicity, cell permeability, or whether the molecule can be made or bought.

What this gives you is a **shortlist**, ranked by one cheap filter, of
molecules worth testing next. It is a starting point for experiments, not a
result.

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

---

## Run it in Google Colab (no install)

`Alpha_Ligand_Pipeline.ipynb` is a notebook that runs the full three-step pipeline
in a browser tab. Nothing is installed on your computer: Colab gives you a free
Linux machine, the notebook installs RDKit, Biopython, meeko and the Linux build of
AutoDock Vina on it, clones this repo, and runs the three scripts.

**Link:** <https://colab.research.google.com/github/adamcantcodelol/alpha-ligand/blob/main/Alpha_Ligand_Pipeline.ipynb>

You click the play button on each cell, top to bottom. The default example is `3ERT`
with active-site residues `A:351,A:383,A:524`, and it docks 4-hydroxytamoxifen (the
drug that really binds that pocket) as a control. The last cell downloads the result
files - reports, the generated molecules and the docked poses - to your computer.

Things to know:

* The install cell takes a few minutes, and must be re-run if Colab disconnects your
  session.
* A free Colab machine is reclaimed when idle. Nothing is lost except time; re-run
  from the top.
* The same honesty applies as everywhere else here: Vina scores are empirical
  estimates, not affinities, and a good score is not a validated binder.

The web app in `docs/` still does step 1 alone, with no notebook at all:
<https://adamcantcodelol.github.io/alpha-ligand/>
