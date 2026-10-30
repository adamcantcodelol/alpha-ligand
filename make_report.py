#!/usr/bin/env python3
"""
make_report.py - step 4 of the alpha-ligand pipeline: the human-readable export.

This script does NOT do any new science. It takes the files the earlier steps
already wrote and turns them into two things a person can actually read:

    output/<NAME>_results.html    one self-contained web page (double-click it)
    output/<NAME>_top_hits.csv    the same table, for Excel / Google Sheets

The HTML file is deliberately SELF-CONTAINED: the styling is embedded in the
file, and every molecule picture is embedded as a base64 PNG (or an inline SVG
if the PNG drawing backend is missing). There are no links to external CSS, JS,
images or CDNs, so the page renders identically with no internet connection and
survives being emailed as an attachment.

Inputs it reads
---------------
  * the pharmacophore JSON from pocket_profiler.py   (required)
  * the candidate molecule SDF from ligand_builder.py (required)
  * the docked pose SDF from dock_ligands.py          (OPTIONAL)

If the docking file is missing the report is still written - it just becomes a
"chemistry only" report, and says so clearly, because without docking there is
no fit score at all.

Usage
-----
    python make_report.py output/3ERT_pharmacophore.json
    python make_report.py output/3ERT_pharmacophore.json \
        --ligands output/3ERT_ligands.sdf \
        --docked  output/3ERT_docked.sdf \
        --outdir  output
"""

import argparse
import base64
import csv
import datetime
import html
import io
import json
import os
import sys

# ---------------------------------------------------------------------------
# RDKit imports. RDKit is the chemistry toolkit; it reads the SDF files, works
# out molecular properties, and draws the pictures.
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem
    from rdkit import RDLogger
    from rdkit.Chem import Descriptors, Draw, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import AllChem
except ImportError:  # pragma: no cover - only hit when RDKit is not installed
    sys.exit(
        "RDKit is not installed in this Python environment.\n"
        "Install it with:  pip install rdkit\n"
        "(In the Colab notebook, step 0 does this for you.)"
    )

RDLogger.DisableLog("rdApp.*")          # keep RDKit's chatter out of the output

# How close two docking scores have to be before we call them a statistical tie.
TIE_WINDOW_KCAL = 1.0


# ===========================================================================
# SECTION 1 - reading the files the earlier steps produced
# ===========================================================================

def load_pharmacophore(path):
    """Read the JSON blueprint written by pocket_profiler.py."""
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def read_sdf(path):
    """
    Read an SDF file into a list of RDKit molecules, skipping anything RDKit
    cannot parse. Returns an empty list if the file does not exist.
    """
    if not path or not os.path.exists(path):
        return []
    mols = []
    supplier = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
    for mol in supplier:
        if mol is not None:
            mols.append(mol)
    return mols


def prop(mol, name, default=None):
    """Safely read a named property off an RDKit molecule."""
    return mol.GetProp(name) if mol.HasProp(name) else default


def float_prop(mol, name, default=None):
    """Same as prop(), but converted to a float where possible."""
    raw = prop(mol, name)
    if raw is None:
        return default
    try:
        return float(raw)
    except (TypeError, ValueError):
        return default


# ===========================================================================
# SECTION 2 - collecting one record per molecule
# ===========================================================================

def descriptor_block(mol):
    """
    Work out the handful of properties a medicinal chemist glances at first.

    MW   - molecular weight in daltons (how big it is)
    cLogP- calculated greasiness; high = oily, low = water-loving
    HBD  - hydrogen-bond donors
    HBA  - hydrogen-bond acceptors
    RotB - rotatable bonds (how floppy it is)
    TPSA - topological polar surface area
    """
    try:
        plain = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    except Exception:
        plain = None
    target = plain if plain is not None else mol
    try:
        return {
            "mw": round(Descriptors.MolWt(target), 1),
            "clogp": round(Descriptors.MolLogP(target), 2),
            "hbd": rdMolDescriptors.CalcNumHBD(target),
            "hba": rdMolDescriptors.CalcNumHBA(target),
            "rotb": rdMolDescriptors.CalcNumRotatableBonds(target),
            "tpsa": round(rdMolDescriptors.CalcTPSA(target), 1),
        }
    except Exception:
        return {"mw": None, "clogp": None, "hbd": None,
                "hba": None, "rotb": None, "tpsa": None}


def collect_molecules(ligand_mols, docked_mols):
    """
    Merge the two SDF files into one list of records, keyed on molecule_name.

    ligand_builder.py writes several conformers of the same molecule, and
    dock_ligands.py writes several poses of the same molecule, so both files
    are de-duplicated down to one entry per distinct molecule.
    """
    records = {}
    order = []

    # --- the candidates, from the builder --------------------------------
    for mol in ligand_mols:
        name = prop(mol, "molecule_name") or prop(mol, "_Name") or "unnamed"
        if name in records:
            continue
        smiles = prop(mol, "smiles")
        if not smiles:
            smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        records[name] = {
            "name": name,
            "smiles": smiles,
            "source": prop(mol, "source", "generated"),
            "checklist": float_prop(mol, "checklist_score_0_100"),
            "chem_fraction": float_prop(mol, "chemistry_match_fraction"),
            "window_fraction": float_prop(mol, "property_window_fraction"),
            "vina": None,
            "poses": [],
            "is_control": False,
            "props": descriptor_block(mol),
        }
        order.append(name)

    # --- the docking results ---------------------------------------------
    for mol in docked_mols:
        name = prop(mol, "molecule_name") or prop(mol, "_Name") or "unnamed"
        score = float_prop(mol, "vina_score_kcal_per_mol")
        best = float_prop(mol, "vina_best_score_kcal_per_mol", score)
        is_control = str(prop(mol, "is_control", "no")).lower() in ("yes", "true", "1")

        if name not in records:
            # A docked molecule with no builder entry - almost always the
            # control ligand, which is supplied directly as SMILES.
            smiles = prop(mol, "smiles")
            if not smiles:
                try:
                    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
                except Exception:
                    smiles = ""
            records[name] = {
                "name": name,
                "smiles": smiles,
                "source": "control" if is_control else "docked only",
                "checklist": None,
                "chem_fraction": None,
                "window_fraction": None,
                "vina": None,
                "poses": [],
                "is_control": is_control,
                "props": descriptor_block(mol),
            }
            order.append(name)

        rec = records[name]
        rec["is_control"] = rec["is_control"] or is_control
        if score is not None:
            rec["poses"].append(score)
        if best is not None:
            rec["vina"] = best if rec["vina"] is None else min(rec["vina"], best)

    for rec in records.values():
        if rec["vina"] is None and rec["poses"]:
            rec["vina"] = min(rec["poses"])
        rec["poses"] = sorted(rec["poses"])

    return [records[n] for n in order]


def rank_molecules(records):
    """
    Split out the control, then rank the candidates.

    Ranking is by docking score when we have one (more negative = better fit),
    and falls back to the chemistry checklist score when docking was not run.
    """
    candidates = [r for r in records if not r["is_control"]]
    controls = [r for r in records if r["is_control"]]
    have_docking = any(r["vina"] is not None for r in candidates)

    if have_docking:
        candidates.sort(key=lambda r: (r["vina"] is None,
                                       r["vina"] if r["vina"] is not None else 0))
    else:
        candidates.sort(key=lambda r: -(r["checklist"] or 0))

    for i, rec in enumerate(candidates, start=1):
        rec["rank"] = i

    # Flag the statistical ties: anything within TIE_WINDOW_KCAL of the best.
    if have_docking and candidates and candidates[0]["vina"] is not None:
        best = candidates[0]["vina"]
        for rec in candidates:
            rec["tied"] = (rec["vina"] is not None
                           and abs(rec["vina"] - best) <= TIE_WINDOW_KCAL)
    else:
        for rec in candidates:
            rec["tied"] = False

    return candidates, controls, have_docking


# ===========================================================================
# SECTION 3 - drawing the molecules, embedded straight into the page
# ===========================================================================

def depiction(smiles, width=300, height=200):
    """
    Draw a 2D picture of a molecule and return it as a snippet of HTML that
    carries the image INSIDE itself - either a base64 <img> or an inline <svg>.
    Nothing is written to a separate image file, which is what keeps the
    report a single shareable document.
    """
    if not smiles:
        return '<div class="nopic">no structure</div>'
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return '<div class="nopic">could not draw</div>'
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass

    # Preferred: a PNG, base64-encoded directly into the src attribute.
    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        opts.bondLineWidth = 2
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        png = drawer.GetDrawingText()
        b64 = base64.b64encode(png).decode("ascii")
        return ('<img class="mol" alt="structure" '
                'src="data:image/png;base64,%s">' % b64)
    except Exception:
        pass

    # Fallback: inline SVG (no Cairo build of RDKit available).
    try:
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.drawOptions().clearBackground = True
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("<?xml version='1.0' encoding='iso-8859-1'?>", "")
        return '<div class="mol molsvg">%s</div>' % svg
    except Exception:
        return '<div class="nopic">could not draw</div>'


# ===========================================================================
# SECTION 4 - summarising what the pocket asked for
# ===========================================================================

FEATURE_PLAIN_ENGLISH = {
    "cationic": ("a positively charged group",
                 "the pocket wall is lined with negative residues (Asp/Glu), "
                 "so a positive group on the ligand is pulled in"),
    "anionic": ("a negatively charged group",
                "the pocket has positive residues (Arg/Lys) looking for a partner"),
    "aromatic": ("a flat ring",
                 "aromatic residues (Phe/Trp/Tyr) stack face-to-face with flat rings"),
    "hydrophobic": ("a greasy, non-polar part",
                    "greasy residues (Leu/Ile/Val/Met) prefer to touch something greasy"),
    "hbond_donor": ("a hydrogen-bond donor (an -OH or -NH)",
                    "there are acceptor atoms on the wall waiting for a hydrogen"),
    "hbond_acceptor": ("a hydrogen-bond acceptor (an O or N lone pair)",
                       "there are -OH/-NH groups on the wall offering a hydrogen"),
    "anionic_or_aromatic": ("either a negative group or a flat ring",
                            "histidine can do either, depending on its protonation"),
    "open_space": ("room to spare",
                   "glycine leaves a gap the ligand may grow into"),
}


def pocket_demand(data):
    """
    Add up the priority weights pocket_profiler.py gave each kind of chemistry,
    so the report can say "this pocket mostly wants X" with a number behind it.
    """
    totals = {}
    examples = {}
    for feat in data.get("pharmacophore_features", []):
        kind = feat.get("ligand_feature_type", "other")
        totals[kind] = totals.get(kind, 0.0) + float(feat.get("priority_score", 0) or 0)
        if kind not in examples:
            examples[kind] = feat.get("suggested_ligand_chemistry", "")
    ordered = sorted(totals.items(), key=lambda kv: -kv[1])
    return [(k, v, examples.get(k, "")) for k, v in ordered]


# ===========================================================================
# SECTION 5 - the HTML itself
# ===========================================================================

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
"""


def esc(x):
    return html.escape(str(x), quote=True)


def fmt(v, nd=2, dash="—"):
    return dash if v is None else ("%.*f" % (nd, v))


def score_class(vina, control_score):
    """Colour a docking number: green if it is near the control, amber otherwise."""
    if vina is None:
        return ""
    if control_score is not None and vina <= control_score + 1.0:
        return "good"
    if vina <= -8.0:
        return "good"
    return "mid"


def hit_card(rec, control_score, have_docking):
    badges = []
    if rec.get("is_control"):
        badges.append('<span class="badge ctrl">benchmark — real known ligand</span>')
    else:
        badges.append('<span class="badge rank">#%d</span>' % rec.get("rank", 0))
        if rec.get("tied"):
            badges.append('<span class="badge tie">statistically tied for top</span>')
    if rec.get("source"):
        badges.append('<span class="badge">%s</span>' % esc(rec["source"]))

    p = rec["props"]
    scores = []
    if have_docking:
        scores.append(
            '<div class="score"><span class="k">docking score</span>'
            '<span class="v %s">%s</span>'
            '<span class="k">kcal/mol &middot; more negative = better fit</span></div>'
            % (score_class(rec["vina"], control_score), fmt(rec["vina"], 2))
        )
    if rec["checklist"] is not None:
        scores.append(
            '<div class="score"><span class="k">chemistry checklist</span>'
            '<span class="v">%s</span>'
            '<span class="k">out of 100 &middot; how much of the pocket\'s wish list it covers</span></div>'
            % fmt(rec["checklist"], 1)
        )
    if rec["poses"]:
        poses = ", ".join("%.2f" % s for s in rec["poses"][:5])
    else:
        poses = None

    return """
<div class="hit">
  <div class="pic">%(pic)s</div>
  <div>
    <p class="hitname">%(name)s</p>
    <div class="badges">%(badges)s</div>
    <div class="scorerow">%(scores)s</div>
    <div class="props">
      <span>MW <b>%(mw)s</b> Da</span><span>cLogP <b>%(clogp)s</b></span>
      <span>H-bond donors <b>%(hbd)s</b></span><span>H-bond acceptors <b>%(hba)s</b></span>
      <span>rotatable bonds <b>%(rotb)s</b></span><span>TPSA <b>%(tpsa)s</b></span>
    </div>
    <h3>SMILES — paste this into a vendor catalogue search</h3>
    <div class="smiles">%(smiles)s</div>
    %(posesline)s
  </div>
</div>""" % {
        "pic": depiction(rec["smiles"]),
        "name": esc(rec["name"]),
        "badges": "".join(badges),
        "scores": "".join(scores) or '<span class="note">no scores available</span>',
        "mw": fmt(p["mw"], 1), "clogp": fmt(p["clogp"], 2),
        "hbd": p["hbd"] if p["hbd"] is not None else "—",
        "hba": p["hba"] if p["hba"] is not None else "—",
        "rotb": p["rotb"] if p["rotb"] is not None else "—",
        "tpsa": fmt(p["tpsa"], 1),
        "smiles": esc(rec["smiles"]),
        "posesline": ('<p class="note">All poses kept for this molecule: %s kcal/mol.</p>' % poses)
                     if poses else "",
    }


def build_html(data, candidates, controls, have_docking, args):
    code = data.get("structure", "structure")
    residues = ", ".join(data.get("active_site_residues", [])) or "not recorded"
    cutoff = data.get("cutoff_A", "?")
    profile = data.get("profile", {})
    geom = data.get("geometry", {})
    windows = data.get("suggested_ligand_properties", {})
    now = datetime.datetime.now().strftime("%d %B %Y, %H:%M")
    control = controls[0] if controls else None
    control_score = control["vina"] if control else None

    parts = []
    A = parts.append

    A("<!DOCTYPE html>\n<html lang='en'><head><meta charset='utf-8'>")
    A("<meta name='viewport' content='width=device-width,initial-scale=1'>")
    A("<title>Alpha-Ligand results — %s</title>" % esc(code))
    A("<style>%s</style></head><body>" % CSS)

    # ---- header ----------------------------------------------------------
    A("<header><h1>Alpha-Ligand results — %s</h1>" % esc(code))
    A("<p class='tagline'>Candidate repressor molecules designed from the shape and "
      "chemistry of one pocket in this protein, then scored for how well they fit it.</p>")
    A("<div class='meta'>"
      "<span>Structure <b>%s</b></span>"
      "<span>Active-site residues <b>%s</b></span>"
      "<span>Pocket cutoff <b>%s &Aring;</b></span>"
      "<span>Pocket residues <b>%s</b></span>"
      "<span>Generated <b>%s</b></span>"
      "</div></header><main>"
      % (esc(code), esc(residues), esc(cutoff),
         esc(profile.get("n_residues", "?")), esc(now)))

    # ---- what this is ----------------------------------------------------
    A("<div class='card'><h2>What this is</h2>")
    A("<p>A protein does its job using a small pocket on its surface — the active site. "
      "If you can plug that pocket with a small molecule, you switch the protein off. "
      "That molecule is a <b>repressor</b>, and finding one is what this pipeline is for.</p>")
    A("<p>The computer did three things. First it measured the pocket of <b>%s</b> around the "
      "residues you named (%s) and wrote down what chemistry the pocket is asking for — "
      "charges, flat rings, greasy patches. Then it assembled small molecules that carry that "
      "chemistry. "
      % (esc(code), esc(residues)))
    if have_docking:
        A("Finally it used AutoDock Vina to physically place each molecule inside the pocket "
          "and score how well it fits. Those scores are the ranking below.</p>")
    else:
        A("Docking was <b>not</b> run, so nothing below has been physically placed inside the "
          "pocket. The ranking is chemistry-only: it says the molecule carries the right groups, "
          "not that it fits.</p>")
    A("<div class='callout'><b>The one-line summary:</b> this is a shortlist worth testing at "
      "the bench, not a discovery. Most molecules that score well here do not bind in real life. "
      "That is normal, and it is why you test a handful rather than one.</div>")
    A("</div>")

    # ---- how to read this ------------------------------------------------
    A("<div class='card'><h2>How to read this</h2><ul>")
    if have_docking:
        A("<li><b>Docking score (kcal/mol).</b> How well Vina thinks the molecule sits in the "
          "pocket. <b>More negative is better.</b> &minus;9 is a good fit, &minus;5 is a poor one. "
          "Differences smaller than about 1 kcal/mol are inside the noise — treat those molecules "
          "as tied, not ranked. Anything tied for the top is flagged below.</li>")
    A("<li><b>Chemistry checklist (0–100).</b> Purely how much of the pocket's wish list the "
      "molecule's chemistry covers. It is <i>not</i> an affinity and <i>not</i> a docking score. "
      "A molecule can score 95 here and still be far too bulky to get in.</li>")
    A("<li><b>MW</b> is molecular weight — how big it is. <b>cLogP</b> is greasiness: high means "
      "oily, low means water-loving; roughly 0–3 is comfortable for a drug-like molecule. "
      "<b>H-bond donors/acceptors</b> count the places it can hydrogen-bond to the protein. "
      "<b>Rotatable bonds</b> is how floppy it is; floppier molecules pay a bigger penalty on "
      "binding, so fewer is usually better.</li>")
    A("<li><b>SMILES</b> is the molecule written as a line of text. Copy it and paste it into the "
      "search box at a chemical vendor (Enamine, Mcule, Sigma-Aldrich) to see whether that "
      "molecule, or something close to it, is buyable.</li>")
    A("</ul>")
    A("<p class='note'>Want to look at the 3D poses? Open the <code>_docked.sdf</code> file "
      "alongside the protein structure in a free viewer such as ChimeraX or PyMOL.</p></div>")

    # ---- the benchmark ---------------------------------------------------
    A("<div class='card'><h2>Benchmark</h2>")
    if control and control["vina"] is not None:
        A("<p><b>%s</b> is a molecule already known to bind this pocket, so it was docked too, "
          "as a control. It scored <b>%s kcal/mol</b>. That is the yardstick: it shows the setup "
          "is self-consistent, and it tells you how far the generated molecules still have to go.</p>"
          % (esc(control["name"]), fmt(control["vina"], 2)))
        if candidates and candidates[0]["vina"] is not None:
            gap = candidates[0]["vina"] - control["vina"]
            A("<p>The best generated molecule is <b>%.2f kcal/mol behind</b> the known binder. "
              "A gap like that usually means the generated molecules are smaller and do not fill "
              "the pocket as completely.</p>" % gap)
        A(hit_card(control, control_score, have_docking))
    else:
        A("<div class='callout warn'><b>There is no benchmark for this run.</b> No known binder "
          "was supplied as a control, so there is nothing to calibrate against. "
          "<b>The scores below are only comparable with each other</b> — a &minus;8 here does not "
          "mean the same thing as a &minus;8 in a run on a different protein, and it does not tell "
          "you the molecule binds. If a ligand is known for this protein, re-run with it as the "
          "control to get a yardstick.</div>")
    A("</div>")

    # ---- the ranking -----------------------------------------------------
    A("<div class='card'><h2>Ranked candidates</h2>")
    if have_docking:
        tied = [r for r in candidates if r.get("tied")]
        if len(tied) > 1:
            A("<div class='callout warn'><b>%d molecules are within %.1f kcal/mol of the best "
              "score.</b> That is inside Vina's margin of error, so they are statistically tied. "
              "Do not read anything into their exact order — if you can test several, test several "
              "of these.</div>" % (len(tied), TIE_WINDOW_KCAL))
        A("<p class='note'>Ranked by docking score, best fit first.</p>")
    else:
        A("<div class='callout warn'><b>Chemistry-only report.</b> No docking results were found, "
          "so these are ranked by the chemistry checklist alone. Nothing here has been placed "
          "inside the pocket, and none of these numbers is a fit score.</div>")

    # summary table first, cards after
    A("<h3>At a glance</h3><table><thead><tr><th class='num'>#</th><th>Molecule</th>")
    if have_docking:
        A("<th class='num'>Docking<br>kcal/mol</th>")
    A("<th class='num'>Checklist<br>/100</th><th class='num'>MW</th><th class='num'>cLogP</th>"
      "<th class='num'>HBD</th><th class='num'>HBA</th><th class='num'>RotB</th></tr></thead><tbody>")
    for r in candidates:
        p = r["props"]
        A("<tr><td class='num'>%d</td><td>%s%s</td>" % (
            r["rank"], esc(r["name"]),
            " <span class='badge tie'>tied</span>" if r.get("tied") else ""))
        if have_docking:
            A("<td class='num'>%s</td>" % fmt(r["vina"], 2))
        A("<td class='num'>%s</td><td class='num'>%s</td><td class='num'>%s</td>"
          "<td class='num'>%s</td><td class='num'>%s</td><td class='num'>%s</td></tr>"
          % (fmt(r["checklist"], 1), fmt(p["mw"], 1), fmt(p["clogp"], 2),
             p["hbd"] if p["hbd"] is not None else "—",
             p["hba"] if p["hba"] is not None else "—",
             p["rotb"] if p["rotb"] is not None else "—"))
    A("</tbody></table>")

    A("<h3>Each molecule in detail</h3>")
    for r in candidates:
        A(hit_card(r, control_score, have_docking))
    A("</div>")

    # ---- pocket profile --------------------------------------------------
    A("<div class='card'><h2>Why these molecules — the pocket's shopping list</h2>")
    A("<p>Every feature below was read off the pocket wall. The weight is how insistently the "
      "pocket asks for that kind of chemistry: residues you named as the active site count "
      "extra, and so do residues closest to the pocket centre.</p>")
    A("<table><thead><tr><th>What the pocket wants</th><th class='num'>Weight</th>"
      "<th>Why</th><th>Chemistry that provides it</th></tr></thead><tbody>")
    for kind, weight, chem in pocket_demand(data):
        plain, why = FEATURE_PLAIN_ENGLISH.get(kind, (kind.replace("_", " "), ""))
        A("<tr><td><b>%s</b></td><td class='num'>%.1f</td><td>%s</td><td>%s</td></tr>"
          % (esc(plain), weight, esc(why), esc(chem)))
    A("</tbody></table>")

    A("<h3>The pocket in numbers</h3><div class='meta'>"
      "<span>Residues lining it <b>%s</b></span>"
      "<span>Greasy (hydrophobic) <b>%s%%</b></span>"
      "<span>Net charge at pH 7 <b>%s</b></span>"
      "<span>Mean hydropathy <b>%s</b></span>"
      "<span>Widest span <b>%s &Aring;</b></span>"
      "</div>" % (
          esc(profile.get("n_residues", "?")),
          esc(round(100 * profile.get("hydrophobic", {}).get("fraction", 0))),
          esc(profile.get("net_formal_charge_pH7", "?")),
          esc(profile.get("mean_kd_hydropathy", "?")),
          esc(geom.get("max_atom_span_A", "?"))))
    A("<p class='note'>A negative pocket charge is why the shopping list is topped by positive "
      "groups, and a wide span is why the suggested molecules are not tiny.</p>")

    if windows:
        A("<h3>Size and property window the profiler suggested</h3><div class='meta'>")
        mw = windows.get("molecular_weight_Da") or ["?", "?"]
        lp = windows.get("clogp") or ["?", "?"]
        A("<span>MW <b>%s–%s</b> Da</span><span>cLogP <b>%s–%s</b></span>"
          "<span>Rotatable bonds <b>&le; %s</b></span>"
          % (esc(mw[0]), esc(mw[-1]), esc(lp[0]), esc(lp[-1]),
             esc(windows.get("rotatable_bonds_max", "?"))))
        A("</div>")
        if windows.get("basis"):
            A("<p class='note'>%s</p>" % esc(windows["basis"]))
    A("</div>")

    # ---- honesty ---------------------------------------------------------
    A("<div class='card'><h2>What these numbers are not</h2>")
    A("<div class='callout bad'>No measured affinity appears anywhere in this report. "
      "There is no K<sub>d</sub>, K<sub>i</sub> or IC<sub>50</sub> here, because none has been "
      "measured. Anyone who quotes you one from a docking run is guessing.</div><ul>")
    A("<li><b>A docking score is an empirical estimate, not a measurement.</b> Vina's scoring "
      "function is fitted to known structures and correlates only loosely with real binding "
      "strength. Treat it as accurate to roughly &plusmn;2–3 kcal/mol.</li>")
    A("<li><b>A good score is not a validated binder.</b> It means the molecule can be placed in "
      "the pocket without clashing, in a way the scoring function likes. Plenty of molecules that "
      "do that never bind.</li>")
    A("<li><b>The protein was held rigid.</b> Real side chains move to accommodate a ligand, and "
      "whole loops can shift. A molecule rejected here might fit a protein that is allowed to "
      "flex, and vice versa.</li>")
    A("<li><b>Water was stripped out.</b> Crystal waters often bridge a ligand to the protein; "
      "removing them is standard practice but it throws away real interactions.</li>")
    A("<li><b>One protonation state was assumed.</b> Charges on both the ligand and the protein "
      "were assigned by rule at pH 7, and getting one wrong changes the score.</li>")
    A("<li><b>Only the pocket you named was searched.</b> If the residues were wrong, everything "
      "here describes the wrong hole in the protein.</li>")
    A("<li><b>Small differences mean nothing.</b> Inside about 1 kcal/mol the ordering is noise.</li>")
    A("</ul>")
    A("<h3>So what do I actually do with this?</h3><ol>")
    A("<li>Read the ranked list and pick the top few — including the tied ones, not just #1.</li>")
    A("<li>Copy each SMILES into a vendor search (Enamine, Mcule, Sigma-Aldrich) and see what is "
      "buyable, or close enough to buyable.</li>")
    A("<li>Open the poses in ChimeraX or PyMOL and sanity-check them by eye: is the charged group "
      "actually near the residue it is supposed to pair with?</li>")
    A("<li>Order a handful and assay them against your protein. The bench is the only thing here "
      "that decides.</li>")
    A("<li>If nothing is buyable, the pocket's shopping list above is still useful on its own — "
      "it tells you what chemistry to look for in a compound library you already have.</li>")
    A("</ol></div>")

    A("</main><footer>Generated by make_report.py, part of the alpha-ligand pipeline "
      "(pocket_profiler.py &rarr; ligand_builder.py &rarr; dock_ligands.py &rarr; make_report.py). "
      "This file is self-contained: the styling and every picture are embedded inside it, so it "
      "works offline and can be emailed as-is. Structure %s, generated %s."
      "</footer></body></html>" % (esc(code), esc(now)))

    return "\n".join(parts)


# ===========================================================================
# SECTION 6 - the CSV
# ===========================================================================

CSV_COLUMNS = [
    "rank", "molecule_name", "smiles", "docking_score_kcal_per_mol",
    "checklist_score_0_100", "molecular_weight_Da", "clogp",
    "hbond_donors", "hbond_acceptors", "rotatable_bonds", "tpsa",
    "source", "is_control", "tied_with_best",
]


def write_csv(path, candidates, controls):
    """One row per molecule, newest spreadsheet-friendly layout."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(CSV_COLUMNS)
        for rec in list(candidates) + list(controls):
            p = rec["props"]
            writer.writerow([
                rec.get("rank", ""),
                rec["name"],
                rec["smiles"],
                "" if rec["vina"] is None else "%.2f" % rec["vina"],
                "" if rec["checklist"] is None else "%.1f" % rec["checklist"],
                "" if p["mw"] is None else p["mw"],
                "" if p["clogp"] is None else p["clogp"],
                "" if p["hbd"] is None else p["hbd"],
                "" if p["hba"] is None else p["hba"],
                "" if p["rotb"] is None else p["rotb"],
                "" if p["tpsa"] is None else p["tpsa"],
                rec.get("source", ""),
                "yes" if rec.get("is_control") else "no",
                "yes" if rec.get("tied") else "no",
            ])


# ===========================================================================
# SECTION 7 - putting it together
# ===========================================================================

def default_sibling(pharm_path, suffix):
    """Guess output/<NAME>_<suffix> from output/<NAME>_pharmacophore.json."""
    base = os.path.basename(pharm_path)
    name = base[:-len("_pharmacophore.json")] if base.endswith("_pharmacophore.json") \
        else os.path.splitext(base)[0]
    return os.path.join(os.path.dirname(pharm_path) or ".", name + suffix), name


def main():
    ap = argparse.ArgumentParser(
        description="Turn alpha-ligand pipeline output into one readable, "
                    "self-contained HTML page plus a CSV of the top hits.")
    ap.add_argument("pharmacophore_json",
                    help="output/<NAME>_pharmacophore.json from pocket_profiler.py")
    ap.add_argument("--ligands", default=None,
                    help="candidate SDF from ligand_builder.py "
                         "(default: output/<NAME>_ligands.sdf)")
    ap.add_argument("--docked", default=None,
                    help="docked pose SDF from dock_ligands.py "
                         "(default: output/<NAME>_docked.sdf; optional)")
    ap.add_argument("--outdir", default="output",
                    help="where to write the report (default: output)")
    ap.add_argument("--quiet", action="store_true", help="print less")
    args = ap.parse_args()

    if not os.path.exists(args.pharmacophore_json):
        sys.exit("Cannot find pharmacophore file: %s\n"
                 "Run pocket_profiler.py first." % args.pharmacophore_json)

    guess_lig, name = default_sibling(args.pharmacophore_json, "_ligands.sdf")
    guess_dock, _ = default_sibling(args.pharmacophore_json, "_docked.sdf")
    ligand_path = args.ligands or guess_lig
    docked_path = args.docked or guess_dock

    data = load_pharmacophore(args.pharmacophore_json)

    ligand_mols = read_sdf(ligand_path)
    if not ligand_mols:
        sys.exit("No molecules could be read from %s\n"
                 "Run ligand_builder.py first." % ligand_path)

    docked_mols = read_sdf(docked_path)
    if not args.quiet:
        print("pharmacophore : %s" % args.pharmacophore_json)
        print("candidates    : %s  (%d SDF entries)" % (ligand_path, len(ligand_mols)))
        if docked_mols:
            print("docked poses  : %s  (%d SDF entries)" % (docked_path, len(docked_mols)))
        else:
            print("docked poses  : none found at %s - writing a chemistry-only report"
                  % docked_path)

    records = collect_molecules(ligand_mols, docked_mols)
    candidates, controls, have_docking = rank_molecules(records)

    os.makedirs(args.outdir, exist_ok=True)
    html_path = os.path.join(args.outdir, "%s_results.html" % name)
    csv_path = os.path.join(args.outdir, "%s_top_hits.csv" % name)

    page = build_html(data, candidates, controls, have_docking, args)
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(page)
    write_csv(csv_path, candidates, controls)

    if not args.quiet:
        print()
        print("molecules in report : %d%s" % (
            len(candidates),
            " (+1 control: %s)" % controls[0]["name"] if controls else " (no control)"))
        print("ranked by           : %s" % ("docking score" if have_docking
                                            else "chemistry checklist (no docking found)"))
        print("HTML report         : %s  (%.0f KB, self-contained)"
              % (html_path, os.path.getsize(html_path) / 1024.0))
        print("CSV of hits         : %s" % csv_path)
        print()
        print("Double-click the HTML file to open it in any browser. It needs no internet")
        print("and no software, so you can email it to anyone.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
