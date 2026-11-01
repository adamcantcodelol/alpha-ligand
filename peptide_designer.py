#!/usr/bin/env python3
"""
peptide_designer.py
===================

THE PEPTIDE BRANCH of alpha-ligand.

The rest of this project designs a SMALL MOLECULE repressor: profile a pocket,
build molecules, dock them, and you end up with something you buy from a
chemical vendor. This script is the other road. It takes the SAME pharmacophore
JSON from pocket_profiler.py and designs short PEPTIDE (mini-protein) binder
candidates instead - sequences you can express in E. coli from the plasmid that
plasmid_builder.py designs for you, or order as synthetic peptides.

HOW IT WORKS (plain English)
----------------------------
1. pocket_profiler.py already worked out, for every residue lining the pocket,
   what a binder should present back at it ("this wall is positively charged,
   so bring something negative"). Each of those wishes has a priority score.
2. This script translates each of those wishes into the AMINO ACIDS that would
   satisfy it. A pocket feature asking for anionic chemistry asks the peptide
   for Asp/Glu; a feature asking for cationic chemistry asks for Arg/Lys; an
   aromatic feature asks for Phe/Trp/Tyr; and so on.
   The feature list and its scores are read straight from the profiler's JSON
   (which is produced by COMPLEMENT_RULES in pocket_profiler.py), so there is
   exactly one complementarity table in this project and this script does not
   invent a second, contradictory one.
3. It picks a sensible peptide length from the pocket's size and how many
   distinct chemical wishes it has (8-20 residues; a pocket binder peptide
   longer than that is a different design problem).
4. It generates many candidate sequences by sampling from that weighted
   amino-acid demand, scores each on how well its composition covers the
   weighted wish list, applies sanity filters (net charge, GRAVY, no silly
   runs, Cys avoided, no Met at position 2, etc.), and ranks what survives.
5. It reports real, computable properties per candidate with Biopython's
   ProtParam: length, net charge at pH 7, GRAVY, predicted pI, molecular
   weight.
6. It writes a FASTA that plasmid_builder.py accepts directly, and can run
   plasmid_builder.py for you (--plasmid) so you finish with an annotated,
   expression-ready GenBank construct.

WHAT THIS IS NOT - PLEASE READ
------------------------------
This is COMPOSITION-BASED COMPLEMENTARITY. It is a souped-up way of saying
"this pocket is greasy and negative, so a binder probably wants greasy and
positive residues, roughly this many of them."

It does NOT:
  * fold the peptide or predict its 3D structure,
  * dock it, or place any residue anywhere near the pocket,
  * predict binding affinity, Kd, IC50 or whether it binds AT ALL,
  * account for the peptide's backbone conformation, entropy, solubility in
    your buffer, protease stability, or cell permeability.

A peptide from this script is a starting hypothesis to test at the bench, in
the same spirit as the small-molecule shortlist. Real structure-based
mini-protein binder design uses tools like RFdiffusion, ProteinMPNN, AlphaFold2
/ AF-Multimer or Rosetta, which actually build and evaluate a 3D interface.
Nothing here replaces those. Expect most designs from a composition-only method
to do nothing; that is the honest prior.

USAGE
-----
  python peptide_designer.py output/3ERT_pharmacophore.json
  python peptide_designer.py output/3ERT_pharmacophore.json --n 8 --plasmid

Only Biopython is needed (plus DNA Chisel if you use --plasmid, because that
is plasmid_builder.py's dependency, not this script's).
"""

import argparse
import json
import os
import random
import subprocess
import sys
from collections import Counter

from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ---------------------------------------------------------------------------
# SECTION 1: Feature -> amino acid table
#
# The KEYS here are exactly the ligand_feature_type strings that
# pocket_profiler.py's COMPLEMENT_RULES emits. We are not re-deciding what
# complements what - the profiler already did that - we are only saying which
# amino acids can physically present each requested chemistry.
# ---------------------------------------------------------------------------

FEATURE_TO_AA = {
    # profiler asks the binder for a NEGATIVE group (pocket wall is Arg/Lys)
    "anionic": {"D": 1.0, "E": 1.0},
    # profiler asks for a POSITIVE group (pocket wall is Asp/Glu)
    "cationic": {"R": 1.0, "K": 1.0},
    # His wall: either a negative group or a ring to stack with
    "anionic_or_aromatic": {"D": 0.5, "E": 0.5, "F": 0.5, "Y": 0.5, "W": 0.4},
    # aromatic wall wants a ring back (Arg/Lys can also do cation-pi)
    "aromatic": {"F": 1.0, "W": 0.9, "Y": 0.9, "R": 0.3, "K": 0.2},
    # Tyr wall: ring plus an H-bond partner
    "aromatic_and_hbond": {"Y": 1.0, "W": 0.8, "F": 0.6, "H": 0.6,
                           "N": 0.4, "Q": 0.4, "S": 0.3, "T": 0.3},
    # the wall DONATES an H-bond, so we want an acceptor back
    "hbond_acceptor": {"D": 0.8, "E": 0.8, "N": 0.9, "Q": 0.9,
                       "S": 0.7, "T": 0.7, "Y": 0.5},
    # amide wall: donates and accepts, so anything that does both
    "hbond_mixed": {"N": 1.0, "Q": 1.0, "S": 0.8, "T": 0.8, "H": 0.7,
                    "Y": 0.6, "D": 0.4, "R": 0.4},
    # Cys wall. We answer with a plain H-bond acceptor. We deliberately do NOT
    # suggest a Cys in the peptide: see the Cys filter below.
    "hbond_acceptor_or_covalent": {"N": 0.8, "Q": 0.8, "S": 0.7, "T": 0.7,
                                   "D": 0.6, "E": 0.6},
    # greasy wall wants grease back
    "hydrophobic": {"L": 1.0, "I": 1.0, "V": 0.9, "M": 0.7, "A": 0.5, "F": 0.6},
    # Gly wall: there is room here. Small/flexible residues exploit it, and a
    # Gly in the peptide buys backbone freedom to reach in.
    "shape": {"G": 0.8, "A": 0.8, "S": 0.6, "P": 0.4},
}

# Residues we will never place by sampling (reasons given in the report).
BANNED = set("C")            # free Cys: disulfide scrambling, oxidation
# Pro is allowed only sparingly - it wrecks regular backbone geometry.
PRO_MAX_FRACTION = 0.12

KD = {  # Kyte-Doolittle, same numbers pocket_profiler.py uses
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Sanity windows for a soluble, synthesisable, expressible peptide.
MAX_RUN = 3               # no more than 3 identical residues in a row
GRAVY_WINDOW = (-2.0, 1.2)
NET_CHARGE_ABS_MAX = 6    # |net charge at pH 7|
MAX_SINGLE_AA_FRACTION = 0.40


# ---------------------------------------------------------------------------
# SECTION 2: Turning the pharmacophore into an amino-acid "shopping list"
# ---------------------------------------------------------------------------

def load_pharmacophore(path):
    with open(path) as fh:
        return json.load(fh)


def aa_demand(features):
    """Weighted demand for each amino acid, and the weight of each feature
    type, straight from the profiler's priority scores."""
    demand = Counter()
    by_type = Counter()
    unknown = Counter()
    for f in features:
        ftype = f.get("ligand_feature_type")
        w = float(f.get("priority_score", 1.0))
        table = FEATURE_TO_AA.get(ftype)
        if table is None:
            unknown[ftype] += 1
            continue
        by_type[ftype] += w
        for aa, share in table.items():
            if aa in BANNED:
                continue
            demand[aa] += w * share
    return demand, by_type, unknown


def choose_length(data, by_type, user_length=None):
    """Pick a peptide length from pocket span and how many distinct chemical
    demands the pocket makes. Clamped to 8-20 residues."""
    if user_length:
        return max(4, min(60, user_length)), "set by you with --length"

    span = float(data.get("geometry", {}).get("max_atom_span_A", 20.0) or 20.0)
    # An extended peptide advances ~3.3 A per residue; a helix ~1.5 A. Pocket
    # binders are rarely fully extended, so we use ~2.4 A/residue and only
    # cover the half of the span that is actually the binding face.
    by_span = (span * 0.5) / 2.4
    # Each distinct chemical demand wants at least one residue to answer it.
    by_demand = len(by_type) * 1.6
    raw = (by_span + by_demand) / 2.0
    length = int(round(max(8, min(20, raw))))
    why = ("pocket span %.1f A -> ~%.0f residues of binding face; %d distinct "
           "feature types -> ~%.0f residues; averaged and clamped to 8-20"
           % (span, by_span, len(by_type), by_demand))
    return length, why


def target_counts(demand, length):
    """How many of each amino acid the pocket 'wants', scaled to the length."""
    total = sum(demand.values()) or 1.0
    return {aa: (w / total) * length for aa, w in demand.items()}


# ---------------------------------------------------------------------------
# SECTION 3: Generating and scoring candidate sequences
# ---------------------------------------------------------------------------

def coverage_score(seq, targets):
    """0-100: how well this sequence's composition covers the weighted demand.

    For each amino acid the pocket wants n copies of, we credit min(have, n).
    Over-supplying one residue therefore earns nothing extra, which stops the
    optimiser collapsing to 'AAAAAAAA' for an all-hydrophobic pocket.
    """
    have = Counter(seq)
    wanted = sum(targets.values()) or 1.0
    got = sum(min(have.get(aa, 0), n) for aa, n in targets.items())
    return 100.0 * got / wanted


def sanity_check(seq):
    """Return a list of reasons this sequence is unusable (empty = fine)."""
    problems = []
    pa = ProteinAnalysis(seq)

    if any(c in BANNED for c in seq):
        problems.append("contains Cys")

    # runs
    run, prev = 1, ""
    for c in seq:
        run = run + 1 if c == prev else 1
        prev = c
        if run > MAX_RUN:
            problems.append("run of %d identical residues (%s)" % (run, c))
            break

    counts = Counter(seq)
    top_aa, top_n = counts.most_common(1)[0]
    if top_n / len(seq) > MAX_SINGLE_AA_FRACTION:
        problems.append("%.0f%% of the peptide is %s" % (100 * top_n / len(seq), top_aa))

    if counts.get("P", 0) / len(seq) > PRO_MAX_FRACTION:
        problems.append("too many prolines")

    gravy = pa.gravy()
    if not (GRAVY_WINDOW[0] <= gravy <= GRAVY_WINDOW[1]):
        problems.append("GRAVY %.2f outside %.1f..%.1f (solubility risk)"
                        % (gravy, *GRAVY_WINDOW))

    charge = pa.charge_at_pH(7.0)
    if abs(charge) > NET_CHARGE_ABS_MAX:
        problems.append("net charge %+.1f at pH 7 is extreme" % charge)

    # Met at position 2: after the Met that plasmid_builder.py adds at the N
    # terminus, a second Met right behind it invites internal-start and
    # Met-aminopeptidase ambiguity. Cheap to avoid, so we avoid it.
    if seq.startswith("M"):
        problems.append("starts with Met (clashes with the added start Met)")

    # A charged or Pro residue at position 1 after the start Met also governs
    # N-end-rule stability in E. coli; Pro at the N terminus is the worst case.
    if seq.startswith("P"):
        problems.append("starts with Pro")

    return problems


def generate(demand, targets, length, n_keep, tries, rng):
    """Sample many sequences from the weighted demand, keep the best legal
    ones. Sampling (not greedy packing) gives genuinely different candidates
    rather than permutations of one answer."""
    alphabet = [aa for aa in demand if aa not in BANNED]
    weights = [demand[aa] for aa in alphabet]
    if not alphabet:
        raise SystemExit("ERROR: the pharmacophore produced no usable amino-acid demand.")

    seen = {}
    rejected = Counter()
    for _ in range(tries):
        seq = "".join(rng.choices(alphabet, weights=weights, k=length))
        if seq in seen:
            continue
        problems = sanity_check(seq)
        if problems:
            rejected[problems[0].split(" (")[0]] += 1
            continue
        seen[seq] = coverage_score(seq, targets)

    ranked = sorted(seen.items(), key=lambda kv: -kv[1])
    return ranked[:n_keep], len(seen), rejected


def properties(seq):
    pa = ProteinAnalysis(seq)
    return {
        "length": len(seq),
        "mw_Da": round(pa.molecular_weight(), 1),
        "net_charge_pH7": round(pa.charge_at_pH(7.0), 2),
        "pI": round(pa.isoelectric_point(), 2),
        "gravy": round(pa.gravy(), 3),
        "aromatic_fraction": round(sum(seq.count(a) for a in "FWY") / len(seq), 3),
        "hydrophobic_fraction": round(sum(seq.count(a) for a in "AVLIMFW") / len(seq), 3),
        "instability_index": round(pa.instability_index(), 1),
    }


# ---------------------------------------------------------------------------
# SECTION 4: The handoff to plasmid_builder.py
#
# plasmid_builder.py refuses anything shorter than 20 aa (MIN_LEN). Our
# peptides are deliberately 8-20 aa, so most of them would be silently
# skipped. Rather than let that happen we pad the construct up to 20 aa with a
# flexible GS linker on the C-terminal side and SAY SO, both here and in the
# FASTA header. The binder residues are untouched; the linker is inert spacer
# that also keeps the designed face away from the His tag.
# ---------------------------------------------------------------------------

PLASMID_MIN_LEN = 20
LINKER_UNIT = "GSGGS"


def pad_for_expression(seq, min_len=PLASMID_MIN_LEN):
    """Return (padded_sequence, linker_added). Empty linker if already long."""
    if len(seq) >= min_len:
        return seq, ""
    need = min_len - len(seq)
    linker = (LINKER_UNIT * (need // len(LINKER_UNIT) + 1))[:need]
    return seq + linker, linker


# ---------------------------------------------------------------------------
# SECTION 5: Reporting
# ---------------------------------------------------------------------------

HONESTY = """\
WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
These sequences were designed by COMPOSITION-BASED COMPLEMENTARITY: the pocket
profile says which chemistry a binder should present, and this script picked
amino acids that present it, in roughly the right proportions, at a sensible
length.

That is ALL it did. In particular it did NOT:
  * fold the peptide or model its structure,
  * dock it into the pocket or place any residue at any coordinate,
  * compute, estimate or imply a binding affinity, Kd, IC50 or percentage
    chance of working,
  * check protease stability, aggregation, solubility in your buffer, or
    whether the peptide is even capable of reaching the site in a folded
    protein.

Composition is a real signal - binders to greasy, negatively charged pockets
really are enriched in greasy, positively charged residues - but it is a weak
one on its own. Real structure-based mini-protein binder design uses
RFdiffusion, ProteinMPNN, Rosetta and AlphaFold2 / AF-Multimer, which build an
actual 3D interface and score it. If this project ever claims otherwise,
the claim is wrong.

Treat the list as a cheap first panel to express and screen, and expect most
of it to do nothing. That is the normal outcome for a composition-only method,
not a failure of your bench work.
"""


def render_report(data, length, length_why, by_type, demand, targets,
                  designs, n_sampled, rejected, pharm_path):
    L = []
    add = L.append
    name = data.get("structure", "?")
    add("=" * 72)
    add("PEPTIDE BINDER CANDIDATES (composition-complementarity design)")
    add("Structure      : %s" % name)
    add("Pharmacophore  : %s" % pharm_path)
    add("=" * 72)
    add("")
    add("1. WHAT THE POCKET IS ASKING FOR")
    prof = data.get("profile", {})
    add("   Pocket residues      : %s" % prof.get("n_residues"))
    add("   Mean KD hydropathy   : %+.2f (%s)" % (
        prof.get("mean_kd_hydropathy", 0.0),
        "greasy" if prof.get("mean_kd_hydropathy", 0) > 0 else "polar"))
    add("   Net pocket charge    : %+d at pH 7" % prof.get("net_formal_charge_pH7", 0))
    add("   Pocket span          : %s A" % data.get("geometry", {}).get("max_atom_span_A"))
    add("")
    add("   Weighted demand by feature type (from pocket_profiler.py):")
    total_w = sum(by_type.values()) or 1.0
    for ftype, w in by_type.most_common():
        add("     %-28s %6.1f  (%4.1f%%)" % (ftype, w, 100 * w / total_w))
    add("")
    add("   Translated into amino acids the peptide should contain:")
    tot_d = sum(demand.values()) or 1.0
    for aa, w in sorted(demand.items(), key=lambda kv: -kv[1]):
        add("     %s  weight %7.1f  (%4.1f%%)  -> want about %.1f copies"
            % (aa, w, 100 * w / tot_d, targets.get(aa, 0.0)))
    add("")
    add("2. PEPTIDE LENGTH")
    add("   Chosen length        : %d residues" % length)
    add("   Why                  : %s" % length_why)
    add("")
    add("3. CANDIDATES (ranked by composition coverage of the demand above)")
    add("   %-4s %-24s %6s %5s %7s %7s %6s %9s" % (
        "RANK", "SEQUENCE", "COVER", "LEN", "CHARGE", "GRAVY", "pI", "MW (Da)"))
    for i, d in enumerate(designs, 1):
        p = d["properties"]
        add("   %-4d %-24s %5.1f%% %5d %+7.2f %7.2f %6.2f %9.1f" % (
            i, d["sequence"], d["coverage_score"], p["length"],
            p["net_charge_pH7"], p["gravy"], p["pI"], p["mw_Da"]))
    add("")
    add("   COVER  = percent of the pocket's weighted amino-acid demand that")
    add("            this composition covers. It is NOT a binding score.")
    add("   CHARGE = net charge at pH 7 (ProtParam). GRAVY = mean hydropathy;")
    add("            positive is greasy, negative is water-loving.")
    add("")
    add("4. SEARCH SUMMARY")
    add("   Sequences sampled    : %d legal unique (from the weighted demand)" % n_sampled)
    if rejected:
        add("   Rejected by filters  :")
        for reason, n in rejected.most_common(8):
            add("     %-52s %d" % (reason, n))
    add("   Filters applied      : no Cys; no run of >%d identical residues;" % MAX_RUN)
    add("                          no single residue >%.0f%% of the peptide;"
        % (100 * MAX_SINGLE_AA_FRACTION))
    add("                          GRAVY within %.1f..%.1f; |net charge| <= %d;"
        % (GRAVY_WINDOW[0], GRAVY_WINDOW[1], NET_CHARGE_ABS_MAX))
    add("                          no leading Met (a start Met gets added later)")
    add("                          and no leading Pro.")
    add("")
    add("5. EXPRESSION HANDOFF")
    add("   plasmid_builder.py only accepts proteins of %d aa or more, so any" % PLASMID_MIN_LEN)
    add("   design shorter than that is padded to %d aa with a flexible" % PLASMID_MIN_LEN)
    add("   %s linker at the C terminus before the handoff. The designed" % LINKER_UNIT)
    add("   residues are never changed; the linker is inert spacer that also")
    add("   keeps the binder face away from the His tag. Padding is recorded in")
    add("   every FASTA header.")
    add("")
    add(HONESTY)
    return "\n".join(L)


# ---------------------------------------------------------------------------
# SECTION 6: Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Design candidate peptide binders from a pocket_profiler.py "
                    "pharmacophore JSON, and hand them to plasmid_builder.py.")
    ap.add_argument("pharmacophore_json",
                    help="output/<NAME>_pharmacophore.json from pocket_profiler.py")
    ap.add_argument("--n", type=int, default=10, help="how many candidates to keep (default 10)")
    ap.add_argument("--length", type=int, default=None,
                    help="force a peptide length (default: chosen from pocket size)")
    ap.add_argument("--tries", type=int, default=40000,
                    help="how many sequences to sample (default 40000)")
    ap.add_argument("--seed", type=int, default=0, help="random seed (default 0, reproducible)")
    ap.add_argument("--outdir", default="output", help="output folder (default: output)")
    ap.add_argument("--plasmid", action="store_true",
                    help="also run plasmid_builder.py on the designed peptides")
    ap.add_argument("--python", default=sys.executable,
                    help="python used to run plasmid_builder.py with --plasmid")
    args = ap.parse_args()

    if not os.path.isfile(args.pharmacophore_json):
        raise SystemExit("ERROR: cannot find %s\nRun pocket_profiler.py first."
                         % args.pharmacophore_json)

    data = load_pharmacophore(args.pharmacophore_json)
    features = data.get("pharmacophore_features") or []
    if not features:
        raise SystemExit("ERROR: no pharmacophore_features in %s" % args.pharmacophore_json)

    demand, by_type, unknown = aa_demand(features)
    if unknown:
        for ftype, n in unknown.items():
            print("WARNING: %d feature(s) of unknown type '%s' were ignored. "
                  "If pocket_profiler.py gained a new feature type, add it to "
                  "FEATURE_TO_AA in peptide_designer.py." % (n, ftype))

    length, length_why = choose_length(data, by_type, args.length)
    targets = target_counts(demand, length)

    rng = random.Random(args.seed)
    ranked, n_sampled, rejected = generate(demand, targets, length, args.n, args.tries, rng)
    if not ranked:
        raise SystemExit("ERROR: no sequence survived the sanity filters. Try "
                         "--tries 200000, or a different --length.")

    name = data.get("structure", "peptide")
    designs = []
    for i, (seq, score) in enumerate(ranked, 1):
        padded, linker = pad_for_expression(seq)
        designs.append({
            "rank": i,
            "id": "%s_pep%02d" % (name, i),
            "sequence": seq,
            "coverage_score": round(score, 1),
            "properties": properties(seq),
            "expression_sequence": padded,
            "linker_added": linker,
            "expression_properties": properties(padded),
        })

    report = render_report(data, length, length_why, by_type, demand, targets,
                           designs, n_sampled, rejected, args.pharmacophore_json)
    print(report)

    os.makedirs(args.outdir, exist_ok=True)
    txt_path = os.path.join(args.outdir, "%s_peptide_report.txt" % name)
    json_path = os.path.join(args.outdir, "%s_peptides.json" % name)
    fasta_path = os.path.join(args.outdir, "%s_peptides.fasta" % name)

    with open(txt_path, "w") as fh:
        fh.write(report + "\n")

    with open(json_path, "w") as fh:
        json.dump({
            "structure": name,
            "source_pharmacophore": args.pharmacophore_json,
            "design_length": length,
            "length_rationale": length_why,
            "feature_type_weights": dict(by_type),
            "amino_acid_demand": {k: round(v, 2) for k, v in demand.items()},
            "target_counts": {k: round(v, 2) for k, v in targets.items()},
            "plasmid_min_length": PLASMID_MIN_LEN,
            "linker_unit": LINKER_UNIT,
            "designs": designs,
            "method": "composition-based complementarity; no folding, no docking, "
                      "no affinity prediction",
            "disclaimer": " ".join(HONESTY.split()),
        }, fh, indent=2)

    with open(fasta_path, "w") as fh:
        for d in designs:
            note = ("binder %s + %s linker pad to %d aa"
                    % (d["sequence"], d["linker_added"], PLASMID_MIN_LEN)
                    ) if d["linker_added"] else ("binder %s, no padding needed" % d["sequence"])
            fh.write(">%s rank=%d coverage=%.1f%% %s\n%s\n"
                     % (d["id"], d["rank"], d["coverage_score"], note,
                        d["expression_sequence"]))

    print("Wrote %s" % txt_path)
    print("Wrote %s" % json_path)
    print("Wrote %s  (ready for plasmid_builder.py)" % fasta_path)
    print("")
    print("NEXT STEP - turn these into expression constructs:")
    print("  python plasmid_builder.py %s" % fasta_path)

    if args.plasmid:
        print("")
        print("=" * 72)
        print("RUNNING plasmid_builder.py ON %s" % fasta_path)
        print("=" * 72)
        here = os.path.dirname(os.path.abspath(__file__))
        cmd = [args.python, os.path.join(here, "plasmid_builder.py"),
               fasta_path, "--outdir", args.outdir]
        # Capture and re-print, so the child output lands in the right place
        # even when stdout is redirected to a file or a notebook cell.
        proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True)
        print(proc.stdout, end="")
        rc = proc.returncode
        if rc != 0:
            print("plasmid_builder.py exited with code %d" % rc)
            return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())
