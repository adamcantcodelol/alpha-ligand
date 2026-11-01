# binder-app

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adamcantcodelol/alpha-ligand/blob/main/Alpha_Ligand_Pipeline.ipynb)

**Want the whole pipeline with nothing to install?** Click the badge above - it opens
`Alpha_Ligand_Pipeline.ipynb` in Google Colab and runs all three steps (profile the
pocket, build molecules, dock them, and write the readable report) on Google's computers, in your browser. See
[Run it in Google Colab](#run-it-in-google-colab-no-install) below.

Two small, beginner-friendly Python tools:

1. **`pocket_profiler.py`** – reads a protein structure, looks at an active site, and works backwards to a blueprint ("pharmacophore") for a ligand that should bind there and repress the protein.
2. **`ligand_builder.py`** – turns that blueprint into real molecules with 3D coordinates.
3. **`dock_ligands.py`** – places each molecule in the pocket with AutoDock Vina and scores the fit.
4. **`make_report.py`** – turns all of that into **one readable HTML page** (plus a CSV) that you can double-click, print, or email to someone who has none of this installed.
5. **`plasmid_builder.py`** – takes protein sequences in a FASTA file and designs codon-optimised, His-tagged E. coli genes ready for NdeI/XhoI cloning.

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

## make_report.py - the readable results page

The three steps above leave you with `.json`, `.sdf` and `.txt` files. Useful to a
computer, not to a person. `make_report.py` reads them all back and writes two files
that a human can actually use:

| File | What it is |
| --- | --- |
| `output/<NAME>_results.html` | **The one you read.** A single web page: ranked molecules with pictures of their structures, docking scores explained in plain English, the properties that matter, the SMILES string to paste into a vendor search, why the pocket asked for that chemistry, and an honest section on what the numbers are not. |
| `output/<NAME>_top_hits.csv` | The same table for Excel or Google Sheets: rank, name, SMILES, docking score, checklist score, MW, cLogP, HBD, HBA, rotatable bonds, TPSA. |

### Running it

```bash
# after pocket_profiler.py, ligand_builder.py and (optionally) dock_ligands.py
python make_report.py output/3ERT_pharmacophore.json
```

Windows (PowerShell), from the repo folder:

```powershell
.\venv\Scripts\python.exe make_report.py output\3ERT_pharmacophore.json
```

It finds `output/3ERT_ligands.sdf` and `output/3ERT_docked.sdf` next to the JSON on
its own. Point it somewhere else with `--ligands`, `--docked` and `--outdir`.

**If you never ran the docking step**, that is fine - the report is still written, it
just becomes a *chemistry-only* report, ranked by the checklist score, and it says so
in a yellow box at the top so nobody mistakes it for a fit prediction.

### Why the HTML file is one file

Everything is embedded inside it: the stylesheet, and every molecule picture as a
base64 PNG drawn by RDKit's `rdMolDraw2D`. There are no `<link>`s, no `<script>`s, no
CDNs and no separate image folder. So it:

* opens by double-clicking, on any computer, with no software and no internet;
* survives being emailed as an attachment or dropped in Slack;
* prints sensibly - the dark theme switches to a light one in a `@media print` block;
* is a permanent record of which protein, which residues and which settings gave
  those hits.

### What it deliberately does not contain

No Kd, no Ki, no IC50, no "predicted affinity in nM". None of those have been
measured, and a docking score cannot be converted into one. The report says this
on the page, in a section headed *What these numbers are not*, along with the
rigid-receptor, stripped-waters and +/-2-3 kcal/mol caveats. Molecules whose scores
are within 1 kcal/mol of the best are flagged as **statistically tied**, because
inside that window the ordering is noise.

### Honest limits

* The pictures are 2D sketches of the molecule, not the docked pose. For the pose,
  open `_docked.sdf` next to the structure in ChimeraX or PyMOL.
* Properties (MW, cLogP, HBD, HBA, rotatable bonds, TPSA) are RDKit calculations
  from the 2D structure, not measurements.
* If there was no control ligand, the report says so loudly: without a known binder
  to calibrate against, the scores are only comparable with each other.

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

You click the play button on each cell, top to bottom. **No code editing is
needed**: Step 2 is a Colab form with real text boxes and sliders for

* `PDB_CODE` (or tick `UPLOAD_MY_OWN_FILE` to upload your own `.pdb`/`.cif`),
* `ACTIVE_SITE_RESIDUES` as `chain:number`, comma separated,
* `CUTOFF_ANGSTROMS`, `NUMBER_OF_MOLECULES`, `CONFORMERS_PER_MOLECULE`,
* `CONTROL_LIGAND` (`auto` / `none` / `custom` + a SMILES box).

The boxes start on the `3ERT` example (`A:351,A:383,A:524`); overwrite them for
your own protein. Step 2b inspects the structure and suggests active-site
residues from whatever is bound in the crystal, printed in the exact
`chain:number` format to copy and paste. Step 2c echoes back the protein,
residues, cutoff and control that are about to be used, and stops with a plain
message if the residues box is empty or malformed. The last cell downloads the
result files - reports, the generated molecules and the docked poses - to your
computer.

Things to know:

* The control is tied to the protein: `auto` uses 4-hydroxytamoxifen **only**
  when the code really is `3ERT`, and no control at all otherwise, since a
  control from a different protein is meaningless. Use `custom` to supply one.
* Step 2b's suggestions come from the crystal contents, not from biology. If the
  structure is apo (nothing bound) it says so and points you at the literature or
  the RCSB entry page.
* The install cell takes a few minutes, and must be re-run if Colab disconnects your
  session.
* A free Colab machine is reclaimed when idle. Nothing is lost except time; re-run
  from the top.
* The same honesty applies as everywhere else here: Vina scores are empirical
  estimates, not affinities, and a good score is not a validated binder.

The web app in `docs/` still does step 1 alone, with no notebook at all:
<https://adamcantcodelol.github.io/alpha-ligand/>

---

## The peptide branch - `peptide_designer.py` and `make_peptide_report.py`

Everything above designs a **small molecule**: something you buy from a chemical
vendor. This branch is the other road. It takes the **same** pharmacophore JSON
from `pocket_profiler.py` and designs short **peptide** binder candidates, then
hands them to `plasmid_builder.py` so you finish with an expression-ready,
annotated GenBank insert.

| | Small molecule branch | Peptide branch |
|---|---|---|
| Output | ranked molecules + docking scores | ranked sequences + a cloneable DNA insert |
| How you obtain it | buy it from a catalogue | express it in *E. coli*, or order it synthesised |
| Docked? | yes (AutoDock Vina) | **no** - see the honest limits below |
| Scripts | `ligand_builder.py`, `dock_ligands.py`, `make_report.py` | `peptide_designer.py`, `plasmid_builder.py`, `make_peptide_report.py` |

### Running it

Linux / macOS:

```bash
./venv/bin/python peptide_designer.py output/3ERT_pharmacophore.json --n 10 --plasmid
./venv/bin/python make_peptide_report.py output/3ERT_peptides.json --constructs output
```

Windows (PowerShell):

```powershell
.\venv\Scripts\python peptide_designer.py output\3ERT_pharmacophore.json --n 10 --plasmid
.\venv\Scripts\python make_peptide_report.py output\3ERT_peptides.json --constructs output
```

`--plasmid` runs `plasmid_builder.py` for you. Without it, the one-line handoff
is printed at the end of the run and is simply:

```bash
python plasmid_builder.py output/3ERT_peptides.fasta
```

Needs Biopython, plus DNA Chisel for the plasmid step (`pip install biopython dnachisel`).
No RDKit, no Vina - the peptide branch is pure Biopython.

### Options for `peptide_designer.py`

| Flag | Default | What it does |
|---|---|---|
| `--n` | 10 | how many candidates to keep |
| `--length` | auto | force a peptide length; default is chosen from pocket span and feature count, clamped to 8-20 aa |
| `--tries` | 40000 | how many sequences to sample |
| `--seed` | 0 | random seed, so runs are reproducible |
| `--plasmid` | off | also run `plasmid_builder.py` on the designed FASTA |
| `--outdir` | `output` | where files go |

### How it works, in plain English

1. `pocket_profiler.py` already decided, residue by residue, what a binder should
   present back at the pocket, with a priority score per wish. This script reads
   that list straight out of the JSON. **The complementarity table lives in one
   place only** (`COMPLEMENT_RULES` in `pocket_profiler.py`); `peptide_designer.py`
   only maps each of those feature types onto the amino acids that can present
   that chemistry - `anionic` -> Asp/Glu, `cationic` -> Arg/Lys, `aromatic` ->
   Phe/Trp/Tyr, `hydrophobic` -> Leu/Ile/Val/Met, the H-bond types ->
   Ser/Thr/Asn/Gln/His, and so on. There is no second, contradictory table.
2. Length comes from the pocket: roughly half the pocket span divided by ~2.4 A
   per residue, averaged with ~1.6 residues per distinct chemical demand, clamped
   to 8-20 aa. For 3ERT that gives 11.
3. It samples tens of thousands of sequences from that weighted demand and scores
   each on **coverage**: the share of the pocket's weighted amino-acid demand the
   composition satisfies. Over-supplying one residue earns nothing extra, so an
   all-hydrophobic pocket does not collapse the answer to `LLLLLLLL`.
4. Sanity filters throw out anything unusable: free Cys, runs of more than 3
   identical residues, any residue over 40% of the peptide, too many prolines,
   GRAVY outside -2.0..1.2, |net charge at pH 7| above 6, a leading Met (a start
   Met gets added later, and two in a row invites internal-start and
   Met-aminopeptidase ambiguity) and a leading Pro.
5. Properties are computed with Biopython's `ProtParam`: length, net charge at
   pH 7, GRAVY, predicted pI, molecular weight, instability index.

### The 20 aa floor, handled openly

`plasmid_builder.py` refuses proteins shorter than 20 aa, and these peptides are
deliberately 8-20 aa, so most would be silently skipped. Instead of letting that
happen, `peptide_designer.py` pads any short design up to 20 aa with an inert,
flexible `GSGGS` linker on the C-terminal side, which also keeps the designed face
away from the His tag. The designed residues are never altered, and the padding is
recorded in the FASTA header, the JSON, the text report and the HTML page.

### `make_peptide_report.py`

The peptide twin of `make_report.py`, kept as a separate file so neither branch
can break the other. One self-contained HTML page: sequences in copyable,
chemistry-coloured form, every property, a per-peptide check that the construct's
DNA really translates back to the design, what to do at the bench, and the
limitations. Plus a CSV for Excel.

### Honest limits - please read this bit

This is **composition-based complementarity, not structure-based binder design.**
It does not fold the peptide, does not dock it, does not place any residue at any
coordinate, and computes **no** affinity, Kd, IC50 or probability of binding.
No such number appears anywhere in the output, because none was calculated.

* **Residue order is essentially arbitrary.** The score only sees composition, so
  two sequences with the same residues score identically however they are
  arranged. A scrambled version of a design is the right negative control, and it
  is a control this method cannot pass by construction.
* **A peptide may not reach the pocket at all.** Sites that comfortably fit a
  small molecule are often too enclosed for a peptide backbone. 3ERT's
  tamoxifen pocket is a good example of a deeply buried site.
* Solubility, aggregation, protease sensitivity and cell permeability are not
  modelled. GRAVY and the instability index are crude proxies.
* Cys is excluded by default, so disulfide-stapled and cyclic designs are out of
  scope.
* The pocket profile came from one crystal structure and the residues you named.
  Wrong residues in, wrong peptides out.
* Real structure-based mini-protein binder design uses **RFdiffusion,
  ProteinMPNN, Rosetta and AlphaFold2 / AF-Multimer**, which build and score an
  actual 3D interface. Nothing here substitutes for them. The honest prior is
  that most composition-only designs do nothing; a negative panel should push you
  toward those tools rather than toward more sampling.

### Tested end to end on 3ERT

`peptide_designer.py output/3ERT_pharmacophore.json --n 8 --plasmid` designs 11-aa
peptides (top hit `RLMWDEIVAKF`, 78.4% coverage), pads each to 20 aa, and
`plasmid_builder.py` builds all 8 as 99 bp NdeI/XhoI inserts whose DNA translates
back to `M` + peptide + `LEHHHHHH` exactly - 8 of 8 with all checks passed.

---

## Worked example 2: SHP-2 (PTPN11) and the SHP099 allosteric tunnel

The 3ERT example above uses 4-hydroxytamoxifen as a known-binder benchmark.
This second example does the same thing for **SHP-2**, a phosphatase whose
best-validated drug site is not its catalytic pocket but an **allosteric
tunnel** at the junction of the N-SH2, C-SH2 and PTP domains. The clinical-stage
chemistry (SHP099 and its descendants) binds there and locks the enzyme in its
closed, inactive state, so it is the right site to aim a repressor at.

### Where the site came from (not from memory)

* **Structure of the complex:** PDB **5EHR**, "SHP2 in complex with allosteric
  inhibitor SHP099", 1.7 Å — <https://www.rcsb.org/structure/5EHR>. The bound
  inhibitor is chemical component **5OD**.
* **Paper:** Chen *et al.*, *Nature* **535**, 148–152 (2016), "Allosteric
  inhibition of SHP2 phosphatase inhibits cancers driven by receptor tyrosine
  kinases" — <https://www.nature.com/articles/nature18621>. It reports hydrogen
  bonds to Arg111 and Phe113 on the SH2–SH2 linker and to Glu250 of the PTP
  domain, with the dichlorophenyl group packing against Leu254, Gln257, Pro491
  and Gln495.
* **Apo / autoinhibited structure:** PDB **2SHP** (Hof *et al.*, 1998) —
  <https://www.rcsb.org/structure/2SHP>.
* **SHP099 chemistry:** PubChem CID **118238298** —
  <https://pubchem.ncbi.nlm.nih.gov/compound/118238298> — and the RCSB entry for
  component 5OD — <https://www.rcsb.org/ligand/5OD>.
  SMILES `CC1(N)CCN(CC1)c2cnc(c(N)n2)c3cccc(Cl)c3Cl`
  (6-(4-amino-4-methylpiperidin-1-yl)-3-(2,3-dichlorophenyl)pyrazin-2-amine,
  C16H19Cl2N5, MW 352.3). PubChem and RCSB agree, and RDKit reproduces
  InChIKey `YGUFCDOEKKVKJK-UHFFFAOYSA-N` from both, so the structure is confirmed
  rather than assumed.

Rather than trusting any of those residue lists, `validate_shp2_site.py`
recomputes them from the coordinates:

```bash
./venv/bin/python validate_shp2_site.py data/pdb/5EHR.pdb data/pdb/2SHP.pdb
```

Every amino acid with an atom within 5 Å of SHP099 in 5EHR chain A:

```
A:108,A:110,A:111,A:112,A:113,A:114,A:218,A:219,A:249,A:250,A:253,A:254,A:257,A:489,A:491,A:492,A:495,A:511
```

(Thr108, Glu110, Arg111, Trp112, Phe113, His114, Thr218, Thr219, Glu249,
Glu250, Thr253, Leu254, Gln257, Asp489, Pro491, Lys492, Gln495, Tyr511.)
That contains every residue the paper names, plus Trp112, His114, Asp489 and
Tyr511 which are in contact but not singled out in the text — so computation
and literature agree.

### Numbering check (5EHR vs 2SHP)

Both deposit human PTPN11 with the same numbering, but the script verifies it
instead of assuming: over the **470 residue numbers present in both**, the amino
acid identity matches at **468 (99.6 %)**. The only two differences are
**41 F→L** and **513 F→S**, neither anywhere near the tunnel. **All 18 site
residues are present in 2SHP with identical identities.** The same residue
numbers therefore mean the same residues in both files, and the site definition
transfers directly.

### Is the tunnel even there in the autoinhibited structure?

This matters, because 2SHP has nothing bound — the tunnel might be collapsed.
Test: superpose 5EHR onto 2SHP (CA rmsd **0.72 Å** over 470 residues — the two
structures are nearly identical), carry SHP099 across, and look for clashes.

Result: the ligand fits, except that **four of its atoms clash with the Gln257
side chain** (closest contact 0.91 Å). **No backbone clashes at all.** Both
chains of 2SHP behave the same way.

So the honest answer is: **the tunnel is formed in the autoinhibited structure —
the protein backbone already leaves the space — but one side chain, Gln257, sits
in a rotamer that partly fills it.** That is ordinary induced fit, and it has a
direct practical consequence: **rigid docking into 2SHP will systematically
under-score anything that binds there**, including SHP099 itself.

### Running the pipeline on both structures

```bash
mkdir -p data/pdb
curl -sSL -o data/pdb/5EHR.pdb https://files.rcsb.org/download/5EHR.pdb
curl -sSL -o data/pdb/2SHP.pdb https://files.rcsb.org/download/2SHP.pdb

# the *_allo.pdb files are just copies of the downloads, so the pipeline
# names its outputs 5EHR_allo_* / 2SHP_allo_* and does not overwrite earlier runs
cp data/pdb/5EHR.pdb data/pdb/5EHR_allo.pdb
cp data/pdb/2SHP.pdb data/pdb/2SHP_allo.pdb

SITE="A:108,A:110,A:111,A:112,A:113,A:114,A:218,A:219,A:249,A:250,A:253,A:254,A:257,A:489,A:491,A:492,A:495,A:511"

for S in 5EHR 2SHP; do
  ./venv/bin/python pocket_profiler.py data/pdb/${S}_allo.pdb "$SITE" --cutoff 5
  ./venv/bin/python ligand_builder.py output/${S}_allo_pharmacophore.json --molecules 10
  ./venv/bin/python dock_ligands.py output/${S}_allo_pharmacophore.json \
      --ligands output/${S}_allo_ligands.sdf --receptor data/pdb/${S}.pdb \
      --chain A --exhaustiveness 16 --max-box 24 \
      --control-smiles 'CC1(N)CCN(CC1)c2cnc(c(N)n2)c3cccc(Cl)c3Cl' \
      --control-name 'SHP099_(known_allosteric_binder)'
  ./venv/bin/python make_report.py output/${S}_allo_pharmacophore.json
done
```

Both were profiled on purpose. 5EHR is the ligand-bound structure, so its tunnel
is definitely open and it gives the cleanest benchmark; 2SHP is the structure
Adam actually cares about. Comparing them is the point.

The two pocket profiles come out closely similar (5EHR: 67 shell residues, 40 %
hydrophobic, net −1; 2SHP: 76 shell residues, 42 % hydrophobic, net 0), which
independently supports the conclusion that the site exists in both.

### Results

| Molecule | 5EHR (open tunnel) | 2SHP (autoinhibited) |
|---|---|---|
| **SHP099 (known binder, control)** | **−9.95** | **−8.68** |
| para-benzene + guanidinomethyl / tetrazolylmethyl | −9.13 | −7.95 |
| meta-benzene + guanidinomethyl / tetrazolylmethyl | −8.52 | −7.55 |
| naphthalene + guanidinomethyl | −8.08 | −7.66 |
| para-benzene + guanidinomethyl / acetic-acid | −7.82 | −7.04 |
| quinoline + guanidinomethyl | −7.81 | −7.45 |
| indole + guanidinomethyl | −7.49 | −6.97 |
| indole-phenol + guanidinomethyl | −7.40 | −7.25 |
| meta-benzene + guanidinomethyl / acetic-acid | −7.38 | −6.96 |
| benzofuran + guanidinomethyl | −7.64 | −7.08 |
| phenylalanine-like | −6.27 | −5.90 |

All scores in kcal/mol; more negative is a better predicted fit.

**The docking setup validates.** Re-docking SHP099 into its own structure
(5EHR) reproduces the crystallographic pose at **1.18 Å heavy-atom RMSD** —
under the conventional 2 Å success threshold, so the box, the receptor prep and
the scoring function are doing their job. Cross-docked into 2SHP the best pose
lands **2.53 Å** from the crystal pose, slightly displaced, exactly as the
Gln257 occlusion predicts.

**The control wins, as it should.** SHP099 outscores every generated molecule in
both structures. The gap to the best generated fragment is **0.82 kcal/mol** in
5EHR and **0.73 kcal/mol** in 2SHP. That is a far smaller gap than the 2.4
kcal/mol seen for 3ERT/tamoxifen — but read it carefully, because it mostly
reflects the tunnel being a smaller, more fragment-friendly pocket, not the
generated molecules being good.

**And the 2SHP-vs-5EHR shift is systematic.** Every single molecule scores
0.4–1.3 kcal/mol worse in 2SHP than in 5EHR, control included. That is the
Gln257 rotamer, not chemistry. The practical rule: **compare scores only within
one structure**, and prefer 5EHR when you want numbers to mean something.

### Limitations specific to this example

* A ~0.8 kcal/mol gap to SHP099 does **not** mean these fragments are nearly as
  good as a 71 nM clinical compound. Vina's error is roughly ±2 kcal/mol, the
  gap is inside it, and scoring functions reward raw size and greasiness, which
  flatters charged fragments in a polar tunnel.
* The top hits are guanidine-bearing fragments of 165–240 Da. SHP099 is 352 Da
  with a dichlorophenyl group doing most of the real work. These are starting
  points, not analogues.
* The receptor is rigid. The Gln257 finding above is precisely a case where that
  assumption costs you, and it argues for flexible-side-chain docking or using
  5EHR as the receptor.
* Docking into the tunnel of an already-closed structure cannot tell you whether
  a molecule would *stabilise* the closed state, which is SHP099's actual
  mechanism. Only an assay can.
* The top six generated molecules are within ~1 kcal/mol of each other and
  should be treated as tied.

In the Colab notebook, `5EHR` and `2SHP` are now in `NATIVE_CONTROLS`, so
choosing either with the control set to `auto` docks SHP099 as the benchmark
automatically. For 2SHP that control is a **known binder from another
structure**, not a ligand present in that crystal — still a valid yardstick, but
worth knowing.
