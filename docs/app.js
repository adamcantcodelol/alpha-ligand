/* =====================================================================
 * app.js - the page. All the science lives in core.js; this file only
 * wires up the buttons, draws the results, and makes the downloads.
 * ===================================================================== */
(function () {
  "use strict";

  var Core = window.PocketCore;

  // Everything we currently have loaded.
  var state = { residues: null, name: null, result: null, picked: [] };

  function $(id) { return document.getElementById(id); }

  function setStatus(el, message, kind) {
    el.textContent = message || "";
    el.className = "status" + (kind ? " " + kind : "");
  }

  /* ---------------- Step 1: loading a structure ---------------- */

  // RCSB serves these with permissive CORS headers, so the browser can read
  // them directly. We try the PDB format first and fall back to mmCIF, which
  // is the only option for very large structures.
  var SOURCES = [
    { url: function (c) { return "https://files.rcsb.org/download/" + c + ".pdb"; }, ext: ".pdb" },
    { url: function (c) { return "https://files.rcsb.org/download/" + c + ".cif"; }, ext: ".cif" },
  ];

  function loadFromCode() {
    var code = $("pdbcode").value.trim().toUpperCase();
    var status = $("loadStatus");
    if (!/^[0-9][A-Z0-9]{3}$/.test(code)) {
      setStatus(status, "A PDB code is four characters starting with a digit, like 3ERT or 1LMB.", "error");
      return;
    }
    setStatus(status, "Fetching " + code + " from RCSB\u2026", "working");
    $("fetchBtn").disabled = true;

    var attempt = 0;
    function tryNext() {
      if (attempt >= SOURCES.length) {
        $("fetchBtn").disabled = false;
        setStatus(status,
          "Could not download " + code + ".\nCheck the code is real, and that you are online. " +
          "If your network blocks rcsb.org, download the file from rcsb.org yourself and use the " +
          "file box above instead.", "error");
        return;
      }
      var src = SOURCES[attempt++];
      fetch(src.url(code))
        .then(function (r) {
          if (!r.ok) throw new Error("HTTP " + r.status);
          return r.text();
        })
        .then(function (text) { accept(text, code + src.ext, code); })
        .catch(function () { tryNext(); });
    }
    tryNext();
  }

  function loadFromFile(file) {
    var status = $("loadStatus");
    setStatus(status, "Reading " + file.name + "\u2026", "working");
    var reader = new FileReader();
    reader.onerror = function () { setStatus(status, "Could not read that file.", "error"); };
    reader.onload = function () {
      accept(String(reader.result), file.name, file.name.replace(/\.[^.]+$/, ""));
    };
    reader.readAsText(file);
  }

  // Parse whatever we just got and open up step 2.
  function accept(text, filename, name) {
    var status = $("loadStatus");
    $("fetchBtn").disabled = false;
    var residues;
    try {
      residues = Core.parseStructure(text, filename);
    } catch (err) {
      setStatus(status, "Could not read that structure: " + err.message, "error");
      return;
    }
    var aminoAcids = residues.filter(function (r) { return r.letter; });
    if (aminoAcids.length === 0) {
      setStatus(status, "That file parsed, but it contains no standard amino acids.", "error");
      return;
    }
    state.residues = residues;
    state.name = name;
    state.result = null;
    state.picked = [];
    $("results").hidden = true;

    var chains = Core.chainsPresent(aminoAcids);
    setStatus(status, "Loaded " + name + ": " + aminoAcids.length + " amino-acid residues across " +
      chains.length + " chain" + (chains.length === 1 ? "" : "s") + " (" + chains.join(", ") + ").", "ok");

    var sel = $("chainSelect");
    sel.innerHTML = "";
    chains.forEach(function (c) {
      var opt = document.createElement("option");
      opt.value = c; opt.textContent = "Chain " + c;
      sel.appendChild(opt);
    });
    $("siteCard").hidden = false;
    drawResidueList();
    $("siteCard").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  /* ---------------- Step 2: the residue picker ---------------- */

  function drawResidueList() {
    var chain = $("chainSelect").value;
    var filter = $("resFilter").value.trim().toUpperCase();
    var list = $("resList");
    list.innerHTML = "";
    state.residues.forEach(function (r) {
      if (!r.letter || r.chain !== chain) return;
      var key = Core.residueKey(r);
      var label = r.resname + " " + r.seq;
      if (filter && label.indexOf(filter) === -1 && key.toUpperCase().indexOf(filter) === -1) return;
      var b = document.createElement("button");
      b.type = "button";
      b.textContent = r.letter + r.seq;
      b.title = label + " (" + key + ")";
      if (state.picked.indexOf(key) !== -1) b.className = "on";
      b.onclick = function () { togglePick(key); };
      list.appendChild(b);
    });
    if (!list.children.length) {
      var p = document.createElement("span");
      p.className = "hint";
      p.textContent = "Nothing matches that filter.";
      list.appendChild(p);
    }
  }

  function togglePick(key) {
    var i = state.picked.indexOf(key);
    if (i === -1) state.picked.push(key); else state.picked.splice(i, 1);
    $("residues").value = state.picked.join(",");
    drawResidueList();
  }

  // Typing in the box keeps the clickable list in sync.
  function syncPickedFromBox() {
    state.picked = $("residues").value.split(",")
      .map(function (s) { return s.trim(); })
      .filter(function (s) { return s.length > 0; });
    drawResidueList();
  }

  /* ---------------- Step 3: run the analysis ---------------- */

  function run() {
    var status = $("runStatus");
    var spec = $("residues").value.trim();
    if (!spec) {
      setStatus(status, "Name at least one active-site residue first.", "error");
      return;
    }
    var cutoff = parseFloat($("cutoff").value);
    setStatus(status, "Profiling\u2026", "working");

    // Let the browser paint the "Profiling..." message before we block it.
    setTimeout(function () {
      var result;
      try {
        result = Core.analyse(state.residues, spec, cutoff);
      } catch (err) {
        setStatus(status, err.message, "error");
        $("results").hidden = true;
        return;
      }
      state.result = result;
      setStatus(status, "Done. " + result.profile.n_residues + " residues in the pocket shell.", "ok");
      draw(result);
      $("results").hidden = false;
      $("results").scrollIntoView({ behavior: "smooth", block: "start" });
    }, 20);
  }

  /* ---------------- Drawing the results ---------------- */

  function cell(key, value, sub) {
    var d = document.createElement("div");
    d.className = "cell";
    var k = document.createElement("div"); k.className = "k"; k.textContent = key;
    var v = document.createElement("div"); v.className = "v"; v.textContent = value;
    d.appendChild(k); d.appendChild(v);
    if (sub) { var s = document.createElement("div"); s.className = "s"; s.textContent = sub; d.appendChild(s); }
    return d;
  }

  function pct(d) { return Math.round(100 * d.fraction) + "%"; }

  function draw(result) {
    var p = result.profile, g = result.geometry, props = result.properties;

    var grid = $("summaryGrid");
    grid.innerHTML = "";
    grid.appendChild(cell("Pocket residues", p.n_residues, "within " + result.cutoff.toFixed(1) + " \u00C5"));
    grid.appendChild(cell("Hydrophobic", pct(p.hydrophobic), p.hydrophobic.count + " residues"));
    grid.appendChild(cell("Aromatic", pct(p.aromatic), p.aromatic.count + " residues"));
    grid.appendChild(cell("Positive (K/R)", pct(p.positive), p.positive.count + " residues"));
    grid.appendChild(cell("Negative (D/E)", pct(p.negative), p.negative.count + " residues"));
    grid.appendChild(cell("Polar neutral", pct(p.polar_neutral), p.polar_neutral.count + " residues"));
    grid.appendChild(cell("Mean hydropathy",
      (p.mean_kd_hydropathy >= 0 ? "+" : "") + p.mean_kd_hydropathy.toFixed(2),
      p.mean_kd_hydropathy > 0 ? "greasy pocket" : "polar pocket"));
    grid.appendChild(cell("Net charge at pH 7",
      (p.net_formal_charge_pH7 >= 0 ? "+" : "") + p.net_formal_charge_pH7, "formal charge"));
    grid.appendChild(cell("Side-chain H-bonds",
      p.sidechain_hbond_donors + " / " + p.sidechain_hbond_acceptors, "donors / acceptors"));
    grid.appendChild(cell("Composition",
      p.composition.slice(0, 5).map(function (kv) { return kv[0] + kv[1]; }).join("  "),
      p.composition.map(function (kv) { return kv[0] + "=" + kv[1]; }).join(", ")));

    $("shellList").textContent = result.shell.map(Core.residueKey).sort().join(", ");

    var gg = $("geomGrid");
    gg.innerHTML = "";
    gg.appendChild(cell("Centroid", g.centroid.join(", "), "x, y, z"));
    gg.appendChild(cell("Bounding box", g.bounding_box_xyz.join(" \u00D7 "), "\u00C5"));
    gg.appendChild(cell("Longest span", g.max_atom_span_A + " \u00C5", "atom to atom"));
    gg.appendChild(cell("Pocket atoms", g.n_pocket_atoms, "counted"));

    var tbody = $("featureTable").querySelector("tbody");
    tbody.innerHTML = "";
    result.features.forEach(function (f) {
      var tr = document.createElement("tr");
      if (f.is_named_active_site) tr.className = "named";
      var star = document.createElement("td"); star.className = "star";
      star.textContent = f.is_named_active_site ? "\u2605" : "";
      var res = document.createElement("td");
      res.innerHTML = "<strong>" + f.pocket_residue + "</strong>";
      res.title = f.rationale;
      var type = document.createElement("td");
      type.innerHTML = '<span class="tag ' + f.ligand_feature_type + '">' +
        f.ligand_feature_type.replace(/_/g, " ") + "</span>";
      var chem = document.createElement("td"); chem.textContent = f.suggested_ligand_chemistry;
      var xyz = document.createElement("td"); xyz.className = "xyz"; xyz.textContent = f.approx_xyz.join(", ");
      var sc = document.createElement("td"); sc.className = "score"; sc.textContent = f.priority_score;
      [star, res, type, chem, xyz, sc].forEach(function (td) { tr.appendChild(td); });
      tbody.appendChild(tr);
    });

    var pg = $("propGrid");
    pg.innerHTML = "";
    pg.appendChild(cell("Molecular weight", props.molecular_weight_Da.join(" \u2013 "), "Da"));
    pg.appendChild(cell("cLogP", props.clogp.join(" \u2013 "), "lipophilicity"));
    pg.appendChild(cell("Rotatable bonds", "\u2264 " + props.rotatable_bonds_max, "flexibility"));
    pg.appendChild(cell("H-bond donors", props.hbond_donors.join(" \u2013 "), ""));
    pg.appendChild(cell("H-bond acceptors", props.hbond_acceptors.join(" \u2013 "), ""));
    pg.appendChild(cell("Net charge",
      (props.suggested_net_charge >= 0 ? "+" : "") + props.suggested_net_charge, "opposes the pocket"));
    $("propBasis").textContent = props.basis;
    $("disclaimer").textContent = Core.DISCLAIMER;
  }

  /* ---------------- Downloads, made in the browser ---------------- */

  function download(filename, text, mime) {
    var blob = new Blob([text], { type: mime });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  /* ---------------- Wiring ---------------- */

  window.addEventListener("DOMContentLoaded", function () {
    $("fetchBtn").onclick = loadFromCode;
    $("pdbcode").addEventListener("keydown", function (e) { if (e.key === "Enter") loadFromCode(); });
    $("fileinput").onchange = function (e) { if (e.target.files[0]) loadFromFile(e.target.files[0]); };

    $("chainSelect").onchange = drawResidueList;
    $("resFilter").oninput = drawResidueList;
    $("residues").oninput = syncPickedFromBox;

    var slider = $("cutoff");
    slider.oninput = function () { $("cutoffOut").textContent = parseFloat(slider.value).toFixed(1) + " \u00C5"; };
    slider.oninput();

    $("runBtn").onclick = run;

    $("dlTxt").onclick = function () {
      if (!state.result) return;
      download(state.name + "_pocket_report.txt",
        Core.renderReport(state.name, state.result, 1e9) + "\n", "text/plain");
    };
    $("dlJson").onclick = function () {
      if (!state.result) return;
      download(state.name + "_pharmacophore.json",
        JSON.stringify(Core.toJSON(state.name, state.result), null, 2), "application/json");
    };
  });
})();
