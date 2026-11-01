#!/usr/bin/env python3
"""
make_peptide_report.py
======================

Turns peptide_designer.py's output into ONE self-contained HTML page you can
double-click, read without any software, and email to a labmate.

It is the peptide-branch twin of make_report.py (which does the same job for
the small-molecule branch). It is a separate file on purpose: the two branches
produce different data, and keeping them apart means neither can break the
other.

  python make_peptide_report.py output/3ERT_peptides.json

Optional: --constructs output  (folder where plasmid_builder.py wrote the
<id>.gb / <id>.dna.fasta files) so the page can confirm, per peptide, that an
expression construct exists and that its DNA really translates back to the
designed sequence.

Writes output/<NAME>_peptide_results.html and output/<NAME>_peptides.csv.
Nothing in the page is fetched from the internet.
"""

import argparse
import csv
import html
import json
import os
import sys

CSS = """
:root{
  --bg:#0f1115; --card:#171a21; --line:#272c37; --ink:#e7eaf0; --dim:#9aa3b2;
  --accent:#5fd3a6; --accent2:#7aa7ff; --warn:#f0b254; --bad:#ef6a6a;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
header{padding:40px 24px 4px;max-width:1000px;margin:0 auto}
h1{margin:0 0 6px;font-size:32px;letter-spacing:-.02em}
.tagline{margin:0;color:var(--dim);max-width:70ch}
main{max-width:1000px;margin:0 auto;padding:20px 24px 60px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;
  padding:20px 22px;margin:0 0 18px}
h2{margin:0 0 12px;font-size:19px}
h3{margin:18px 0 6px;font-size:13px;text-transform:uppercase;letter-spacing:.06em;color:var(--dim)}
p{max-width:78ch}
code,.smiles{background:#0f1218;border:1px solid var(--line);border-radius:5px;
  padding:2px 6px;font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:12.5px;
  word-break:break-all}
.meta{display:flex;flex-wrap:wrap;gap:10px 28px;margin:10px 0 0;color:var(--dim);font-size:13.5px}
.meta b{color:var(--ink);font-weight:600}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{text-align:left;padding:8px 10px;border-bottom:1px solid var(--line);vertical-align:top}
th{color:var(--dim);font-weight:600;font-size:12.5px;text-transform:uppercase;letter-spacing:.04em}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.hit{display:grid;grid-template-columns:320px 1fr;gap:18px;border:1px solid var(--line);
  border-radius:11px;padding:16px;margin:0 0 14px;background:#141720;break-inside:avoid}
.pic{background:#fff;border-radius:9px;padding:6px;display:flex;align-items:center;justify-content:center}
img.mol{max-width:100%;height:auto;display:block}
.molsvg{background:#fff;border-radius:9px}
.nopic{color:#888;font-size:13px;padding:40px 0;text-align:center}
.hitname{font-size:16px;font-weight:650;margin:0 0 2px}
.badges{display:flex;flex-wrap:wrap;gap:6px;margin:8px 0}
.badge{font-size:11.5px;padding:2px 8px;border-radius:999px;border:1px solid var(--line);color:var(--dim)}
.badge.rank{background:var(--accent);color:#0f1115;border-color:var(--accent);font-weight:700}
.badge.tie{background:#3a2f16;border-color:var(--warn);color:var(--warn)}
.badge.ctrl{background:#1b2a3f;border-color:var(--accent2);color:var(--accent2)}
.scorerow{display:flex;flex-wrap:wrap;gap:22px;margin:10px 0 4px}
.score .v{font-size:24px;font-weight:700;font-variant-numeric:tabular-nums}
.score .k{display:block;font-size:11.5px;color:var(--dim);text-transform:uppercase;letter-spacing:.05em}
.good{color:var(--accent)} .mid{color:var(--warn)}
.props{display:flex;flex-wrap:wrap;gap:6px 10px;margin:10px 0;font-size:12.5px;color:var(--dim)}
.props span b{color:var(--ink)}
.note{color:var(--dim);font-size:13.5px}
.callout{border-left:3px solid var(--accent2);background:#141b26;padding:12px 16px;border-radius:0 8px 8px 0;margin:12px 0}
.callout.warn{border-left-color:var(--warn);background:#1e1a12}
.callout.bad{border-left-color:var(--bad);background:#1f1414}
ul{max-width:78ch} li{margin:5px 0}
footer{max-width:1000px;margin:0 auto;padding:0 24px 50px;color:var(--dim);font-size:12.5px}
@media print{
  :root{--bg:#fff;--card:#fff;--line:#c8ccd4;--ink:#111;--dim:#555;
        --accent:#0a7a55;--accent2:#1a4ea8;--warn:#8a5a00;--bad:#9a2020}
  body{background:#fff;color:#111;font-size:11pt}
  .card,.hit{background:#fff;border:1px solid #c8ccd4;box-shadow:none}
  .hit{break-inside:avoid;page-break-inside:avoid}
  .badge.rank{background:#eee;color:#111;border-color:#999}
  .callout{background:#f6f7f9}
  header{padding-top:10px}
}

.seqbox{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:15px;
  letter-spacing:.06em;background:#0f1218;border:1px solid var(--line);border-radius:8px;
  padding:10px 12px;word-break:break-all;margin:8px 0;user-select:all}
.linker{color:var(--dim)}
.pep{border:1px solid var(--line);border-radius:11px;padding:16px;margin:0 0 14px;
  background:#141720;break-inside:avoid}
.ok{color:var(--accent)} .no{color:var(--bad)}
"""


AA_CLASS = {
    "D": "neg", "E": "neg", "K": "pos", "R": "pos", "H": "pos",
    "F": "aro", "W": "aro", "Y": "aro",
    "A": "pho", "V": "pho", "L": "pho", "I": "pho", "M": "pho", "P": "pho",
    "S": "pol", "T": "pol", "N": "pol", "Q": "pol", "C": "pol", "G": "gly",
}
CLASS_COLOUR = {"neg": "#ef6a6a", "pos": "#7aa7ff", "aro": "#c792ea",
                "pho": "#f0b254", "pol": "#5fd3a6", "gly": "#9aa3b2"}
CLASS_NAME = {"neg": "negative (Asp/Glu)", "pos": "positive (Lys/Arg/His)",
              "aro": "aromatic (Phe/Trp/Tyr)", "pho": "greasy (Ala/Val/Leu/Ile/Met/Pro)",
              "pol": "polar (Ser/Thr/Asn/Gln/Cys)", "gly": "glycine (flexible)"}

FEATURE_PLAIN = {
    "anionic": "the pocket is positively charged here, so the peptide brings Asp or Glu",
    "cationic": "the pocket is negatively charged here, so the peptide brings Lys or Arg",
    "anionic_or_aromatic": "a histidine wall: answer with a negative charge or a ring",
    "aromatic": "a flat ring in the pocket wall to stack a Phe/Trp/Tyr against",
    "aromatic_and_hbond": "a tyrosine wall: wants a ring plus a hydrogen bond",
    "hbond_acceptor": "the pocket donates a hydrogen bond, so the peptide accepts one",
    "hbond_mixed": "an amide wall that both donates and accepts hydrogen bonds",
    "hbond_acceptor_or_covalent": "a cysteine wall; answered with a plain hydrogen-bond acceptor",
    "hydrophobic": "a greasy wall, answered with greasy residues",
    "shape": "a glycine in the wall leaves space here",
}


def esc(x):
    return html.escape(str(x), quote=True)


def seq_html(seq, linker=""):
    """Colour a sequence by residue class; grey out any appended linker."""
    core = seq[:len(seq) - len(linker)] if linker else seq
    out = []
    for c in core:
        colour = CLASS_COLOUR.get(AA_CLASS.get(c, "gly"), "#e7eaf0")
        out.append('<span style="color:%s">%s</span>' % (colour, esc(c)))
    if linker:
        out.append('<span class="linker">%s</span>' % esc(linker))
    return "".join(out)


def read_constructs(folder, designs):
    """If plasmid_builder.py has been run, check each construct for real."""
    results = {}
    if not folder or not os.path.isdir(folder):
        return results
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
    except ImportError:
        return results
    for d in designs:
        gb = os.path.join(folder, "%s.gb" % d["id"])
        fa = os.path.join(folder, "%s.dna.fasta" % d["id"])
        if not os.path.isfile(gb):
            continue
        entry = {"gb": gb, "fasta": fa if os.path.isfile(fa) else None,
                 "translation_ok": False, "bp": None, "protein": None}
        try:
            rec = SeqIO.read(gb, "genbank")
            entry["bp"] = len(rec.seq)
            cds = [f for f in rec.features if f.type == "CDS"]
            if cds:
                protein = str(Seq(str(cds[0].extract(rec.seq))).translate()).rstrip("*")
                entry["protein"] = protein
                expected = "M" + d["expression_sequence"] + "LEHHHHHH"
                entry["translation_ok"] = (protein == expected)
                entry["expected"] = expected
        except Exception as exc:          # a broken file should not kill the report
            entry["error"] = str(exc)
        results[d["id"]] = entry
    return results


def build_html(data, constructs, args):
    name = data.get("structure", "?")
    designs = data.get("designs", [])
    weights = data.get("feature_type_weights", {})
    total_w = sum(weights.values()) or 1.0

    H = []
    a = H.append
    a("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    a("<meta name='viewport' content='width=device-width,initial-scale=1'>")
    a("<title>%s - peptide binder candidates</title>" % esc(name))
    a("<style>%s</style></head><body>" % CSS)

    a("<header><h1>Peptide binder candidates &mdash; %s</h1>" % esc(name))
    a("<p class='tagline'>Short peptide sequences whose residue chemistry complements "
      "this pocket. Designed from the pocket profile only: nothing here has been folded, "
      "docked or measured.</p>")
    a("<div class='meta'><span>Design length <b>%s aa</b></span>"
      "<span>Candidates <b>%d</b></span>"
      "<span>Source <b>%s</b></span></div></header><main>"
      % (esc(data.get("design_length")), len(designs),
         esc(os.path.basename(data.get("source_pharmacophore", "")))))

    # --- what this is -------------------------------------------------------
    a("<div class='card'><h2>Read this first</h2>")
    a("<div class='callout warn'><p><b>This is a composition-based design, not a "
      "structure-based one.</b> The pocket profile said which chemistry a binder should "
      "present; these sequences present it, in roughly the right proportions, at a "
      "sensible length. The method never folded a peptide, never placed a residue in the "
      "pocket, and produced no affinity estimate. There is no Kd here because none was "
      "calculated.</p></div>")
    a("<p>Real structure-based mini-protein binder design uses RFdiffusion, ProteinMPNN, "
      "Rosetta and AlphaFold2/AF-Multimer, which build an actual 3D interface and score it. "
      "Treat this list as a cheap first panel to express and screen, and expect most of it "
      "to do nothing &mdash; that is the normal outcome for a composition-only method, not "
      "a failure of your bench work.</p>")
    a("<h3>Small molecule or peptide?</h3>")
    a("<ul><li><b>Small molecule</b> (the other branch of this project): you search a "
      "vendor catalogue and buy it. Fast, but you can only buy what exists.</li>"
      "<li><b>Peptide</b> (this page): you express it in <i>E. coli</i> from the plasmid "
      "insert <code>plasmid_builder.py</code> designs, or order it synthesised. Slower, "
      "but any sequence is makeable.</li></ul></div>")

    # --- what the pocket asked for -----------------------------------------
    a("<div class='card'><h2>What the pocket asked for</h2>")
    a("<p class='note'>%s</p>" % esc(data.get("length_rationale", "")))
    a("<table><tr><th>Pocket feature</th><th class='num'>Weight</th>"
      "<th class='num'>Share</th><th>In plain English</th></tr>")
    for ftype, w in sorted(weights.items(), key=lambda kv: -kv[1]):
        a("<tr><td><code>%s</code></td><td class='num'>%.1f</td>"
          "<td class='num'>%.1f%%</td><td class='note'>%s</td></tr>"
          % (esc(ftype), w, 100 * w / total_w,
             esc(FEATURE_PLAIN.get(ftype, "&mdash;"))))
    a("</table>")
    targets = data.get("target_counts", {})
    if targets:
        top = sorted(targets.items(), key=lambda kv: -kv[1])[:10]
        a("<h3>Amino acids the pocket wants most</h3><p>")
        a(", ".join("<b>%s</b> (~%.1f copies)" % (esc(k), v) for k, v in top))
        a("</p>")
    a("</div>")

    # --- candidates ---------------------------------------------------------
    a("<div class='card'><h2>Candidates</h2>")
    a("<p class='note'>Ranked by <b>coverage</b>: the percentage of the pocket's weighted "
      "amino-acid demand this composition satisfies. Coverage is a checklist score, "
      "<b>not</b> a predicted affinity, and candidates within a few points of each other "
      "should be treated as equivalent.</p>")
    a("<p class='note'>Colours: " + " &nbsp;".join(
        "<span style='color:%s'>&#9632;</span> %s" % (CLASS_COLOUR[k], esc(v))
        for k, v in CLASS_NAME.items()) + "</p>")

    for d in designs:
        p = d["properties"]
        a("<div class='pep'>")
        a("<p class='hitname'>#%d &nbsp;<span class='badge rank'>coverage %.1f%%</span> "
          "&nbsp;<span class='badge'>%s</span></p>" % (d["rank"], d["coverage_score"],
                                                       esc(d["id"])))
        a("<h3>Designed binder sequence (click to select, then copy)</h3>")
        a("<div class='seqbox'>%s</div>" % seq_html(d["sequence"]))
        a("<div class='props'>"
          "<span>Length <b>%d aa</b></span>"
          "<span>Net charge pH 7 <b>%+.2f</b></span>"
          "<span>GRAVY <b>%+.2f</b></span>"
          "<span>pI <b>%.2f</b></span>"
          "<span>MW <b>%.1f Da</b></span>"
          "<span>Aromatic <b>%.0f%%</b></span>"
          "<span>Greasy <b>%.0f%%</b></span>"
          "<span>Instability index <b>%.1f</b></span></div>"
          % (p["length"], p["net_charge_pH7"], p["gravy"], p["pI"], p["mw_Da"],
             100 * p["aromatic_fraction"], 100 * p["hydrophobic_fraction"],
             p["instability_index"]))
        if d.get("linker_added"):
            a("<h3>Sequence sent to the plasmid designer</h3>")
            a("<div class='seqbox'>%s</div>" % seq_html(d["expression_sequence"],
                                                        d["linker_added"]))
            a("<p class='note'>The grey tail is an inert <code>%s</code> linker added only "
              "because <code>plasmid_builder.py</code> will not accept a protein shorter "
              "than %s aa. It is flexible spacer and also keeps the designed face away from "
              "the His tag. The designed residues are unchanged.</p>"
              % (esc(d["linker_added"]), esc(data.get("plasmid_min_length", 20))))
        c = constructs.get(d["id"])
        if c:
            if c.get("translation_ok"):
                a("<p class='note'><span class='ok'>&#10003; Expression construct built"
                  "</span> &mdash; %s bp insert in <code>%s</code>; the DNA translates back "
                  "to <code>%s</code> exactly (start Met + designed peptide + LEHHHHHH "
                  "His tag).</p>" % (esc(c["bp"]), esc(os.path.basename(c["gb"])),
                                     esc(c.get("protein", ""))))
            else:
                a("<p class='note'><span class='no'>&#10007; Construct found but its "
                  "translation does not match the design.</span> Do not order it; rerun "
                  "<code>plasmid_builder.py</code> and check the console output.</p>")
        a("</div>")
    a("</div>")

    # --- how to use it ------------------------------------------------------
    a("<div class='card'><h2>What to do with this at the bench</h2><ol>")
    a("<li><b>Pick a few, not one.</b> Coverage differences of a point or two mean "
      "nothing. Taking the top 3&ndash;5 across different charge/GRAVY values gives you a "
      "better panel than taking the single top row.</li>")
    a("<li><b>Decide how to make them.</b> Either order them as synthetic peptides (fast, "
      "costs money per peptide, no cloning) or express them from the insert in the "
      "<code>.gb</code> files (cheap per peptide once cloned, needs an NdeI/XhoI vector "
      "such as a pET backbone, and purification off the His tag).</li>")
    a("<li><b>Screen for binding, do not assume it.</b> A pull-down on the His tag, "
      "fluorescence polarisation, BLI/SPR or a simple thermal shift of your target are all "
      "reasonable first assays. Include a scrambled version of the same composition as your "
      "negative control &mdash; that is the control that tests whether the <i>order</i> of "
      "residues matters, which is exactly what this method does not model.</li>")
    a("<li><b>If nothing binds, that is information.</b> It is the expected result often "
      "enough that a negative panel should push you toward structure-based design "
      "(RFdiffusion/ProteinMPNN) rather than toward more composition sampling.</li>")
    a("</ol></div>")

    a("<div class='card'><h2>Honest limitations</h2><ul>")
    for line in [
        "No folding, no docking, no 3D structure, no energy calculation.",
        "No affinity, Kd, IC50 or probability of binding is computed or implied anywhere.",
        "Residue ORDER is essentially arbitrary: the score only sees composition, so two "
        "sequences with the same residues score identically however they are arranged.",
        "A peptide may not be able to reach a buried pocket at all. Pockets that fit a "
        "small molecule are often too enclosed for a peptide backbone.",
        "Solubility, aggregation, protease sensitivity and cell permeability are not "
        "modelled. GRAVY and the instability index are crude proxies only.",
        "Cys is excluded by default, so disulfide-stapled or cyclic designs are out of "
        "scope here.",
        "The pocket profile itself came from one crystal structure and the residues you "
        "named. Wrong residues in, wrong peptides out.",
    ]:
        a("<li>%s</li>" % esc(line))
    a("</ul></div>")

    a("</main><footer>Generated by make_peptide_report.py (alpha-ligand peptide branch) "
      "from %s. Self-contained: no internet needed.</footer></body></html>"
      % esc(os.path.basename(args.peptides_json)))
    return "\n".join(H)


CSV_COLUMNS = ["rank", "id", "sequence", "coverage_percent", "length",
               "net_charge_pH7", "gravy", "pI", "mw_Da", "instability_index",
               "expression_sequence", "linker_added", "construct_translation_ok"]


def write_csv(path, designs, constructs):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLUMNS)
        for d in designs:
            p = d["properties"]
            c = constructs.get(d["id"], {})
            w.writerow([d["rank"], d["id"], d["sequence"], d["coverage_score"],
                        p["length"], p["net_charge_pH7"], p["gravy"], p["pI"],
                        p["mw_Da"], p["instability_index"],
                        d["expression_sequence"], d["linker_added"] or "",
                        "" if not c else ("yes" if c.get("translation_ok") else "no")])


def main():
    ap = argparse.ArgumentParser(
        description="Readable self-contained HTML page for the peptide branch.")
    ap.add_argument("peptides_json",
                    help="output/<NAME>_peptides.json from peptide_designer.py")
    ap.add_argument("--constructs", default=None,
                    help="folder holding <id>.gb files from plasmid_builder.py "
                         "(default: same folder as the JSON)")
    ap.add_argument("--outdir", default=None,
                    help="where to write the report (default: same folder as the JSON)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.peptides_json):
        sys.exit("Cannot find %s\nRun peptide_designer.py first." % args.peptides_json)

    with open(args.peptides_json) as fh:
        data = json.load(fh)
    designs = data.get("designs", [])
    if not designs:
        sys.exit("No designs in %s" % args.peptides_json)

    folder = os.path.dirname(os.path.abspath(args.peptides_json))
    outdir = args.outdir or folder
    constructs = read_constructs(args.constructs or folder, designs)

    os.makedirs(outdir, exist_ok=True)
    name = data.get("structure", "peptides")
    html_path = os.path.join(outdir, "%s_peptide_results.html" % name)
    csv_path = os.path.join(outdir, "%s_peptides.csv" % name)

    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(build_html(data, constructs, args))
    write_csv(csv_path, designs, constructs)

    if not args.quiet:
        built = sum(1 for c in constructs.values() if c.get("translation_ok"))
        print("peptides in report  : %d" % len(designs))
        print("expression constructs: %d verified (DNA translates back to the design)"
              % built)
        print("HTML report         : %s  (%.0f KB, self-contained)"
              % (html_path, os.path.getsize(html_path) / 1024.0))
        print("CSV                 : %s" % csv_path)
        print()
        print("Double-click the HTML file to open it in any browser.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
