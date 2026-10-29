#!/usr/bin/env python3
"""
plasmid_builder.py
==================

Takes a FASTA file of PROTEIN sequences and, for each one, designs a
codon-optimised E. coli gene ready to be cloned into an NdeI / XhoI
expression vector, with a C-terminal LEHHHHHH (His6) purification tag.

Steps for every sequence:
  1. Validate it (20 standard amino acids only, no stop codon, 20-300 aa).
  2. Add a starting Met if missing, then append the LEHHHHHH tag.
  3. Back-translate and codon-optimise for E. coli with DNA Chisel, while
     avoiding internal NdeI/XhoI sites, keeping GC between 40-60% in every
     50 bp window, and avoiding long single-base runs.
  4. Add cloning flanks: an NdeI site (CATATG) whose ATG *is* the start
     codon, and at the 3' end a TAA stop followed by an XhoI site (CTCGAG).
  5. Print a report and check the DNA translates back to exactly the protein.
  6. Write output/<name>.dna.fasta and output/<name>.gb (GenBank with
     CDS and His-tag features).

No plasmid backbone is invented: the output is only the insert.
"""

import argparse
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import gc_fraction

from dnachisel import (
    DnaOptimizationProblem,
    EnforceTranslation,
    CodonOptimize,
    AvoidPattern,
    EnforceGCContent,
    reverse_translate,
)

# ---------------------------------------------------------------------------
# SECTION 1: Constants
# ---------------------------------------------------------------------------

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")  # the 20 standard amino acids
HIS_TAG = "LEHHHHHH"                       # C-terminal tag we always append
NDEI = "CATATG"                            # NdeI recognition site
XHOI = "CTCGAG"                            # XhoI recognition site
MIN_LEN, MAX_LEN = 20, 300                 # allowed protein length
MAX_HOMOPOLYMER = 6                        # no run longer than this many bases


# ---------------------------------------------------------------------------
# SECTION 2: Validation
# ---------------------------------------------------------------------------

def validate_protein(name, seq):
    """Return (True, cleaned_sequence) if the protein is usable, otherwise
    print a clear reason and return (False, None)."""
    s = str(seq).strip().upper().replace(" ", "")

    if "*" in s or s.endswith("."):
        print("  ERROR: '%s' contains a stop codon ('*'). Remove it and rerun. SKIPPED." % name)
        return False, None

    bad = sorted({c for c in s if c not in STANDARD_AA})
    if bad:
        print("  ERROR: '%s' contains non-standard amino acid letters: %s. SKIPPED."
              % (name, ", ".join(bad)))
        return False, None

    if not (MIN_LEN <= len(s) <= MAX_LEN):
        print("  ERROR: '%s' is %d aa long; must be between %d and %d aa. SKIPPED."
              % (name, len(s), MIN_LEN, MAX_LEN))
        return False, None

    return True, s


def add_met_and_tag(seq):
    """Make sure the protein starts with Met, then append LEHHHHHH."""
    added_met = not seq.startswith("M")
    if added_met:
        seq = "M" + seq
    return seq + HIS_TAG, added_met


# ---------------------------------------------------------------------------
# SECTION 3: Codon optimisation with DNA Chisel
# ---------------------------------------------------------------------------

def optimise_gene(protein):
    """Back-translate the protein then let DNA Chisel fix the sequence so all
    constraints hold. Returns the optimised coding DNA (no flanks yet)."""

    # A first naive DNA sequence: any codon for each amino acid.
    starting_dna = reverse_translate(protein)

    constraints = [
        # Hard requirement: the DNA must still translate to our protein.
        EnforceTranslation(),
        # No internal NdeI / XhoI, so our cloning sites stay unique.
        AvoidPattern(NDEI),
        AvoidPattern(XHOI),
        # GC content between 40% and 60% in every 50 bp window.
        EnforceGCContent(mini=0.40, maxi=0.60, window=50),
    ]
    # Avoid long runs of the same base (e.g. "7xA" means 7 A's in a row).
    for base in "ATGC":
        constraints.append(AvoidPattern("%dx%s" % (MAX_HOMOPOLYMER + 1, base)))

    problem = DnaOptimizationProblem(
        sequence=starting_dna,
        constraints=constraints,
        # Objective (a soft goal): use E. coli's favourite codons.
        objectives=[CodonOptimize(species="e_coli")],
        logger=None,  # keep the terminal quiet
    )
    problem.resolve_constraints()   # make all the hard rules true
    problem.optimize()              # then improve codon usage
    return str(problem.sequence), problem


def add_flanks(coding_dna):
    """Add the cloning sites.

    5' end: NdeI is CATATG and its last three bases ARE the start codon ATG,
            so we put 'CAT' in front of the existing ATG:
            CAT + ATG... = CATATG...  (NdeI cuts CA^TATG)
    3' end: a TAA stop codon, then the XhoI site CTCGAG. The CDS length stays
            a multiple of 3, so the reading frame is preserved.
    """
    assert coding_dna.startswith("ATG"), "coding sequence must start with ATG"
    return "CAT" + coding_dna + "TAA" + XHOI


# ---------------------------------------------------------------------------
# SECTION 4: Reporting and checks
# ---------------------------------------------------------------------------

def report(name, protein_in, tagged_protein, final_dna, added_met):
    """Print all the sanity checks for one designed gene."""
    # The CDS starts at index 3 (after the 'CAT' of the NdeI site) and runs to
    # just before the XhoI site; it includes the TAA stop.
    cds_start = 3
    cds_end = len(final_dna) - len(XHOI)
    cds = final_dna[cds_start:cds_end]
    translated = str(Seq(cds).translate())        # ends with '*' for TAA
    protein_from_dna = translated.rstrip("*")

    gc = 100 * gc_fraction(final_dna)
    n_nde = final_dna.count(NDEI)
    n_xho = final_dna.count(XHOI)
    matches = (protein_from_dna == tagged_protein)

    print("  Input protein length     : %d aa" % len(protein_in))
    print("  Start Met added          : %s" % ("yes" if added_met else "no (already present)"))
    print("  Tagged protein length    : %d aa (+ %s)" % (len(tagged_protein), HIS_TAG))
    print("  Final DNA length         : %d bp (CDS %d bp incl. stop)" % (len(final_dna), len(cds)))
    print("  GC content               : %.1f%%" % gc)
    print("  Reading frame            : CDS is %s multiple of 3" %
          ("a" if len(cds) % 3 == 0 else "NOT a"))
    print("  Translation check        : %s" %
          ("MATCHES the tagged input protein exactly" if matches else "MISMATCH!"))
    print("  NdeI (CATATG) sites      : %d %s" % (n_nde, "(OK, exactly one)" if n_nde == 1 else "(PROBLEM)"))
    print("  XhoI (CTCGAG) sites      : %d %s" % (n_xho, "(OK, exactly one)" if n_xho == 1 else "(PROBLEM)"))
    print("  Note: output is the INSERT only - no plasmid backbone is invented.")
    return cds_start, cds_end, matches and n_nde == 1 and n_xho == 1


# ---------------------------------------------------------------------------
# SECTION 5: Writing files
# ---------------------------------------------------------------------------

def write_outputs(outdir, name, final_dna, tagged_protein, cds_start, cds_end):
    os.makedirs(outdir, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in name)

    # Plain FASTA of the DNA.
    fasta_path = os.path.join(outdir, "%s.dna.fasta" % safe)
    SeqIO.write(
        SeqRecord(Seq(final_dna), id=safe,
                  description="codon-optimized E. coli insert, NdeI/XhoI, C-term %s" % HIS_TAG),
        fasta_path, "fasta")

    # Annotated GenBank file.
    record = SeqRecord(Seq(final_dna), id=safe[:16], name=safe[:16],
                       description="Codon-optimized insert for E. coli expression")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"

    # Where the His-tag sits inside the DNA (last 8 codons before the stop).
    tag_dna_len = len(HIS_TAG) * 3
    tag_start = cds_end - 3 - tag_dna_len   # -3 skips the TAA stop

    record.features = [
        SeqFeature(FeatureLocation(0, len(NDEI)), type="misc_feature",
                   qualifiers={"label": ["NdeI site"], "note": ["CATATG, ATG = start codon"]}),
        SeqFeature(FeatureLocation(cds_start, cds_end), type="CDS",
                   qualifiers={"label": ["%s-%s" % (safe, HIS_TAG)],
                               "translation": [tagged_protein],
                               "codon_start": ["1"],
                               "transl_table": ["11"],
                               "product": ["%s with C-terminal His6 tag" % safe]}),
        SeqFeature(FeatureLocation(tag_start, cds_end - 3), type="misc_feature",
                   qualifiers={"label": ["%s (His6) tag" % HIS_TAG]}),
        SeqFeature(FeatureLocation(cds_end - 3, cds_end), type="misc_feature",
                   qualifiers={"label": ["TAA stop codon"]}),
        SeqFeature(FeatureLocation(len(final_dna) - len(XHOI), len(final_dna)),
                   type="misc_feature", qualifiers={"label": ["XhoI site"]}),
    ]
    gb_path = os.path.join(outdir, "%s.gb" % safe)
    SeqIO.write(record, gb_path, "genbank")
    return fasta_path, gb_path


# ---------------------------------------------------------------------------
# SECTION 6: Main program
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Design codon-optimized, His-tagged E. coli inserts from a protein FASTA.")
    ap.add_argument("fasta", help="Path to a FASTA file of protein sequences")
    ap.add_argument("--outdir", default="output", help="Output folder (default: output)")
    args = ap.parse_args()

    if not os.path.isfile(args.fasta):
        raise SystemExit("ERROR: file not found: %s" % args.fasta)

    records = list(SeqIO.parse(args.fasta, "fasta"))
    if not records:
        raise SystemExit("ERROR: no FASTA sequences found in %s" % args.fasta)

    print("Read %d sequence(s) from %s\n" % (len(records), args.fasta))
    ok_count = 0

    for record in records:
        name = record.id
        print("-" * 68)
        print("Sequence: %s" % name)
        valid, protein = validate_protein(name, record.seq)
        if not valid:
            print("")
            continue

        tagged, added_met = add_met_and_tag(protein)
        try:
            coding_dna, _problem = optimise_gene(tagged)
        except Exception as exc:
            print("  ERROR: DNA Chisel could not satisfy the constraints: %s. SKIPPED." % exc)
            print("")
            continue

        final_dna = add_flanks(coding_dna)
        cds_start, cds_end, all_good = report(name, protein, tagged, final_dna, added_met)
        fasta_path, gb_path = write_outputs(args.outdir, name, final_dna, tagged,
                                            cds_start, cds_end)
        print("  Wrote %s" % fasta_path)
        print("  Wrote %s" % gb_path)
        print("  Status                   : %s" % ("ALL CHECKS PASSED" if all_good else "CHECK FAILED"))
        print("")
        ok_count += 1 if all_good else 0

    print("-" * 68)
    print("Done: %d of %d sequence(s) designed successfully." % (ok_count, len(records)))


if __name__ == "__main__":
    main()
