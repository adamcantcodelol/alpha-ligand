/* =====================================================================
 * app.js - the page. All the science lives in core.js; this file only
 * wires up the buttons, draws the results, and makes the downloads.
 * ===================================================================== */
(function () {
  "use strict";

  var Core = window.PocketCore;

  // Everything we currently have loaded. Shared by both modes, so
  // switching mode never loses the structure, pocket or result.
  var state = {
    residues: null, name: null, result: null, picked: [],
    simplePicked: [], auto: null, mode: "simple", simpleStep: 1, simpleDrawnFor: null,
  };

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

  function loadFromCode(opts) {
    opts = opts || {};
    var simple = !!opts.simple;
    var code = (opts.code !== undefined ? opts.code : $(simple ? "sxCode" : "pdbcode").value)
      .trim().toUpperCase();
    var status = $(simple ? "sxLoadStatus" : "loadStatus");
    var btn = $(simple ? "sxFetch" : "fetchBtn");
    if (!/^[0-9][A-Z0-9]{3}$/.test(code)) {
      setStatus(status, simple
        ? "A PDB code is four characters and starts with a number, like 3ERT or 1LMB."
        : "A PDB code is four characters starting with a digit, like 3ERT or 1LMB.", "error");
      return;
    }
    setStatus(status, simple ? "Looking up " + code + "\u2026" : "Fetching " + code + " from RCSB\u2026", "working");
    btn.disabled = true;

    var attempt = 0;
    function tryNext() {
      if (attempt >= SOURCES.length) {
        btn.disabled = false;
        setStatus(status, simple
          ? "We could not download " + code + ". Check the code is right and that you are online. " +
            "If it still fails, download the file from rcsb.org yourself and use the file box below."
          : "Could not download " + code + ".\nCheck the code is real, and that you are online. " +
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
        .then(function (text) { accept(text, code + src.ext, code, simple, opts.then); })
        .catch(function () { tryNext(); });
    }
    tryNext();
  }

  function loadFromFile(file, simple) {
    var status = $(simple ? "sxLoadStatus" : "loadStatus");
    setStatus(status, "Reading " + file.name + "\u2026", "working");
    var reader = new FileReader();
    reader.onerror = function () { setStatus(status, "Could not read that file.", "error"); };
    reader.onload = function () {
      accept(String(reader.result), file.name, file.name.replace(/\.[^.]+$/, ""), simple);
    };
    reader.readAsText(file);
  }

  // Parse whatever we just got and open up step 2.
  function accept(text, filename, name, simple, then) {
    var status = $(simple ? "sxLoadStatus" : "loadStatus");
    $("fetchBtn").disabled = false;
    $("sxFetch").disabled = false;
    var residues;
    try {
      residues = Core.parseStructure(text, filename);
    } catch (err) {
      setStatus(status, (simple ? "We could not read that structure. " : "Could not read that structure: ") +
        err.message, "error");
      return;
    }
    var aminoAcids = residues.filter(function (r) { return r.letter; });
    if (aminoAcids.length === 0) {
      setStatus(status, simple
        ? "That file opened, but there is no protein in it - no amino acids were found."
        : "That file parsed, but it contains no standard amino acids.", "error");
      return;
    }
    state.residues = residues;
    state.name = name;
    state.result = null;
    state.picked = [];
    state.simplePicked = [];
    state.auto = null;
    $("results").hidden = true;
    $("sxStep3").hidden = true;

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
    var sel2 = $("sxChain");
    sel2.innerHTML = "";
    chains.forEach(function (c) {
      var opt = document.createElement("option");
      opt.value = c; opt.textContent = "Chain " + c;
      sel2.appendChild(opt);
    });

    $("siteCard").hidden = false;
    drawResidueList();
    drawSimpleResidueList();

    // Look for a bound ligand so Simple mode can offer a pocket.
    try {
      state.auto = Core.findLikelyPocket(residues, 5.0, 8);
    } catch (err) {
      state.auto = null;
    }

    if (simple) {
      setStatus(status, "Loaded " + name + ".", "ok");
      showSimpleStep(2);
      if (then) then();
    } else {
      $("siteCard").scrollIntoView({ behavior: "smooth", block: "start" });
    }
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

  function suggestAdvanced() {
    var status = $("autoStatus");
    if (!state.residues) { setStatus(status, "Load a structure first.", "error"); return; }
    if (!state.auto) {
      setStatus(status, "Nothing suitable is bound in this structure, so there is nothing to " +
        "measure a pocket around. Enter the residues yourself.", "error");
      return;
    }
    $("residues").value = state.auto.keys.join(",");
    syncPickedFromBox();
    setStatus(status, "Residues within " + state.auto.radius.toFixed(1) + " \u00C5 of " +
      state.auto.ligand.resname + " " + state.auto.ligand.key + " (" +
      state.auto.ligand.n_atoms + " atoms): " + state.auto.labels.join(", ") + ".", "ok");
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


  /* =====================================================================
   * SIMPLE MODE
   * A three-screen wizard over exactly the same science. It never
   * computes anything the Advanced view does not; it only says less,
   * and says it in words.
   * =================================================================== */

  var MODE_KEY = "alphaligand.mode";

  function setMode(mode, opts) {
    opts = opts || {};
    state.mode = (mode === "advanced") ? "advanced" : "simple";
    var simple = state.mode === "simple";
    $("simpleUI").hidden = !simple;
    $("advancedUI").hidden = simple;
    $("modeSimple").setAttribute("aria-pressed", simple ? "true" : "false");
    $("modeAdvanced").setAttribute("aria-pressed", simple ? "false" : "true");
    $("modeSimple").classList.toggle("on", simple);
    $("modeAdvanced").classList.toggle("on", !simple);
    $("modeHint").textContent = simple
      ? "Step-by-step, in plain English. Switch any time \u2014 you won't lose your place."
      : "Every number, table and download. Switch back any time.";
    $("tagline").textContent = simple
      ? "Point it at a protein and it tells you, in plain English, what kind of molecule would fit the pocket you care about."
      : "Give it a protein structure and its active-site residues. It reads the chemistry of the pocket and works backwards to the blueprint of a ligand that would fill it \u2014 a candidate repressor.";
    try { localStorage.setItem(MODE_KEY, state.mode); } catch (err) { /* private mode */ }

    // Carry whatever is loaded across to the other view.
    if (!simple) {
      if (state.residues) $("siteCard").hidden = false;
      if (state.result) $("results").hidden = false;
    } else {
      if (state.result && state.simpleDrawnFor !== state.result) {
        drawSimple(state.result);
        state.simpleDrawnFor = state.result;
      }
      var dataStep = state.result ? 3 : (state.residues ? 2 : 1);
      showSimpleStep(Math.max(state.simpleStep, dataStep), true);
    }
    if (!opts.silent) window.scrollTo({ top: 0, behavior: "smooth" });
  }

  function showSimpleStep(n, silent) {
    // Never show a step we have no data for.
    if (n >= 2 && !state.residues) n = 1;
    if (n >= 3 && !state.result) n = state.residues ? 2 : 1;
    state.simpleStep = n;
    $("sxStep1").hidden = n !== 1;
    $("sxStep2").hidden = n !== 2;
    $("sxStep3").hidden = n !== 3;
    var items = $("progress").children;
    for (var i = 0; i < items.length; i++) {
      items[i].classList.toggle("on", (i + 1) === n);
      items[i].classList.toggle("done", (i + 1) < n);
    }
    if (n === 2) drawSimpleAuto();
    if (!silent && state.mode === "simple") {
      window.scrollTo({ top: 0, behavior: "smooth" });
    }
  }

  /* ---- step 2: the automatic pocket suggestion ---- */

  function drawSimpleAuto() {
    if (!state.residues) return;
    var aa = state.residues.filter(function (r) { return r.letter; });
    var chains = Core.chainsPresent(aa);
    $("sxLoadedLine").textContent = "Loaded " + state.name + " \u2014 " + aa.length +
      " building blocks in " + chains.length + " chain" + (chains.length === 1 ? "" : "s") +
      " (" + chains.join(", ") + ").";

    var auto = state.auto;
    if (auto) {
      $("sxAutoBox").hidden = false;
      $("sxNoAutoBox").hidden = true;
      $("sxAutoText").textContent =
        "This structure has a molecule called " + auto.ligand.resname + " already sitting in it, " +
        "which almost always marks the spot that matters. The " + auto.keys.length +
        " building blocks closest to it are " + auto.friendly.join(", ") +
        ". We can use those as your pocket.";
      // Pre-fill but do not commit until they press the button.
      if (state.simplePicked.length === 0) $("sxResidues").value = auto.keys.join(",");
      $("sxManualBox").hidden = state.simplePicked.length === 0;
      // Only one obvious button at a time: the big green "use this pocket"
      // is enough until they choose to type their own.
      $("sxRunRow").hidden = $("sxManualBox").hidden;
    } else {
      $("sxAutoBox").hidden = true;
      $("sxNoAutoBox").hidden = false;
      $("sxManualBox").hidden = false;
      $("sxRunRow").hidden = false;
    }
    drawSimpleResidueList();
  }

  function drawSimpleResidueList() {
    if (!state.residues) return;
    var chain = $("sxChain").value;
    var filter = $("sxFilter").value.trim().toUpperCase();
    var list = $("sxResList");
    list.innerHTML = "";
    var picked = $("sxResidues").value.split(",").map(function (t) { return t.trim(); });
    state.residues.forEach(function (r) {
      if (!r.letter || r.chain !== chain) return;
      var key = Core.residueKey(r);
      var label = r.resname + " " + r.seq;
      if (filter && label.indexOf(filter) === -1 && key.toUpperCase().indexOf(filter) === -1) return;
      var b = document.createElement("button");
      b.type = "button";
      b.textContent = r.letter + r.seq;
      b.title = label + " (" + key + ")";
      if (picked.indexOf(key) !== -1) b.className = "on";
      b.onclick = function () { toggleSimplePick(key); };
      list.appendChild(b);
    });
    if (!list.children.length) {
      var p = document.createElement("span");
      p.className = "hint";
      p.textContent = "Nothing matches that filter.";
      list.appendChild(p);
    }
  }

  function toggleSimplePick(key) {
    var picked = $("sxResidues").value.split(",")
      .map(function (t) { return t.trim(); })
      .filter(function (t) { return t.length; });
    var i = picked.indexOf(key);
    if (i === -1) picked.push(key); else picked.splice(i, 1);
    $("sxResidues").value = picked.join(",");
    state.simplePicked = picked;
    drawSimpleResidueList();
  }

  /* ---- step 2 -> 3: run it ---- */

  function runSimple(spec, then) {
    var status = $("sxRunStatus");
    spec = (spec || $("sxResidues").value).trim();
    if (!spec) {
      setStatus(status, "Pick at least one building block first \u2014 either accept the pocket we " +
        "found, or click a few in the list.", "error");
      return;
    }
    setStatus(status, "Reading the pocket\u2026", "working");
    setTimeout(function () {
      var result;
      try {
        result = Core.analyse(state.residues, spec, parseFloat($("cutoff").value));
      } catch (err) {
        setStatus(status, plainError(err.message), "error");
        return;
      }
      state.result = result;
      // Keep Advanced mode in step with Simple mode.
      $("residues").value = spec;
      syncPickedFromBox();
      draw(result);
      setStatus($("runStatus"), "Done. " + result.profile.n_residues + " residues in the pocket shell.", "ok");
      setStatus(status, "", "");
      drawSimple(result);
      state.simpleDrawnFor = result;
      showSimpleStep(3);
      if (then) then();
    }, 20);
  }

  // Core's errors are already friendly-ish; soften the jargon a little.
  function plainError(msg) {
    return msg
      .replace("standard amino acid", "normal protein building block")
      .replace("Try a larger cutoff.", "Try picking a different building block.");
  }

  /* ---- step 3: the plain-English answer ---- */

  function drawSimple(result) {
    var plain = Core.plainSummary(state.name, result);

    var box = $("sxSentences");
    box.innerHTML = "";
    plain.sentences.forEach(function (text, i) {
      var p = document.createElement("p");
      if (i === 0) p.className = "lede";
      p.textContent = text;
      box.appendChild(p);
    });

    var bars = $("sxBars");
    bars.innerHTML = "";
    var title = document.createElement("div");
    title.className = "barstitle";
    title.textContent = "What the pocket wall is made of";
    bars.appendChild(title);
    plain.bars.forEach(function (b) {
      var row = document.createElement("div");
      row.className = "bar";
      var lab = document.createElement("span"); lab.className = "barlabel"; lab.textContent = b.label;
      var track = document.createElement("span"); track.className = "bartrack";
      var fill = document.createElement("span"); fill.className = "barfill";
      fill.style.width = Math.max(1, b.pct) + "%";
      track.appendChild(fill);
      track.title = b.hint;
      var val = document.createElement("span"); val.className = "barval"; val.textContent = b.pct + "%";
      row.appendChild(lab); row.appendChild(track); row.appendChild(val);
      bars.appendChild(row);
    });

    $("sxCaveat").textContent = plain.caveat;

    // The details expander: the same numbers, just folded away.
    var p = result.profile, g = result.geometry, props = result.properties;
    var body = $("sxDetailBody");
    body.innerHTML = "";
    var grid = document.createElement("div");
    grid.className = "grid";
    grid.appendChild(cell("Pocket residues", p.n_residues, "within " + result.cutoff.toFixed(1) + " \u00C5"));
    grid.appendChild(cell("Hydrophobic", pct(p.hydrophobic), p.hydrophobic.count + " residues"));
    grid.appendChild(cell("Aromatic", pct(p.aromatic), p.aromatic.count + " residues"));
    grid.appendChild(cell("Net charge at pH 7",
      (p.net_formal_charge_pH7 >= 0 ? "+" : "") + p.net_formal_charge_pH7, "formal charge"));
    grid.appendChild(cell("Mean hydropathy",
      (p.mean_kd_hydropathy >= 0 ? "+" : "") + p.mean_kd_hydropathy.toFixed(2),
      p.mean_kd_hydropathy > 0 ? "greasy pocket" : "polar pocket"));
    grid.appendChild(cell("Longest span", g.max_atom_span_A + " \u00C5", "atom to atom"));
    grid.appendChild(cell("Molecular weight", props.molecular_weight_Da.join(" \u2013 "), "Da"));
    grid.appendChild(cell("cLogP", props.clogp.join(" \u2013 "), "lipophilicity"));
    grid.appendChild(cell("Suggested net charge",
      (props.suggested_net_charge >= 0 ? "+" : "") + props.suggested_net_charge, "opposes the pocket"));
    body.appendChild(grid);

    var h = document.createElement("h3");
    h.className = "plain";
    h.textContent = "The five strongest things the pocket asks for";
    body.appendChild(h);
    var ul = document.createElement("ul");
    ul.className = "plainlist";
    result.features.slice(0, 5).forEach(function (f) {
      var li = document.createElement("li");
      li.textContent = f.pocket_residue + " \u2192 " + f.suggested_ligand_chemistry +
        " (" + f.rationale + ")";
      ul.appendChild(li);
    });
    body.appendChild(ul);

    var site = result.activeSite.map(function (r) { return Core.residueKey(r); }).join(",");
    $("sxCarryOver").textContent = "When you get there, use PDB code " +
      String(state.name).toUpperCase() + " and residues " + site + " \u2014 the same ones you used here.";
  }

  /* ---------------- Wiring ---------------- */

  window.addEventListener("DOMContentLoaded", function () {
    $("fetchBtn").onclick = function () { loadFromCode({ simple: false }); };
    $("pdbcode").addEventListener("keydown", function (e) { if (e.key === "Enter") loadFromCode({ simple: false }); });
    $("fileinput").onchange = function (e) { if (e.target.files[0]) loadFromFile(e.target.files[0], false); };
    $("autoBtn").onclick = suggestAdvanced;

    $("chainSelect").onchange = drawResidueList;
    $("resFilter").oninput = drawResidueList;
    $("residues").oninput = syncPickedFromBox;

    var slider = $("cutoff");
    slider.oninput = function () { $("cutoffOut").textContent = parseFloat(slider.value).toFixed(1) + " \u00C5"; };
    slider.oninput();

    $("runBtn").onclick = run;

    /* ---- simple mode wiring ---- */
    $("modeSimple").onclick = function () { setMode("simple"); };
    $("modeAdvanced").onclick = function () { setMode("advanced"); };
    $("sxToAdvanced").onclick = function () { setMode("advanced"); };

    $("sxFetch").onclick = function () { loadFromCode({ simple: true }); };
    $("sxCode").addEventListener("keydown", function (e) {
      if (e.key === "Enter") loadFromCode({ simple: true });
    });
    $("sxFile").onchange = function (e) { if (e.target.files[0]) loadFromFile(e.target.files[0], true); };

    // The example: fill everything in and run the whole thing.
    $("sxExample").onclick = function () {
      $("sxCode").value = "3ERT";
      loadFromCode({
        simple: true, code: "3ERT",
        then: function () {
          var spec = state.auto ? state.auto.keys.join(",") : "A:351,A:383,A:524";
          $("sxResidues").value = spec;
          runSimple(spec);
        },
      });
    };

    $("sxUseAuto").onclick = function () {
      if (!state.auto) return;
      $("sxResidues").value = state.auto.keys.join(",");
      runSimple(state.auto.keys.join(","));
    };
    $("sxManualLink").onclick = function () {
      $("sxManualBox").hidden = false;
      $("sxRunRow").hidden = false;
      drawSimpleResidueList();
      $("sxManualBox").scrollIntoView({ behavior: "smooth", block: "center" });
    };
    $("sxChain").onchange = drawSimpleResidueList;
    $("sxFilter").oninput = drawSimpleResidueList;
    $("sxResidues").oninput = function () {
      state.simplePicked = $("sxResidues").value.split(",")
        .map(function (t) { return t.trim(); }).filter(function (t) { return t.length; });
      drawSimpleResidueList();
    };
    $("sxRun").onclick = function () { runSimple(); };
    $("sxBack1").onclick = function () { showSimpleStep(1); };
    $("sxBack2").onclick = function () { showSimpleStep(2); };
    $("sxRestart").onclick = function () {
      state.residues = null; state.result = null; state.auto = null; state.simpleDrawnFor = null;
      state.picked = []; state.simplePicked = [];
      $("sxCode").value = ""; $("sxResidues").value = "";
      $("residues").value = ""; $("results").hidden = true; $("siteCard").hidden = true;
      setStatus($("sxLoadStatus"), "", ""); setStatus($("sxRunStatus"), "", "");
      showSimpleStep(1);
    };
    $("sxDlTxt").onclick = function () {
      if (!state.result) return;
      download(state.name + "_pocket_report.txt",
        Core.renderReport(state.name, state.result, 1e9) + "\n", "text/plain");
    };
    $("sxDlJson").onclick = function () {
      if (!state.result) return;
      download(state.name + "_pharmacophore.json",
        JSON.stringify(Core.toJSON(state.name, state.result), null, 2), "application/json");
    };

    // Simple is the default for a first-time visitor.
    var saved = null;
    try { saved = localStorage.getItem(MODE_KEY); } catch (err) { saved = null; }
    setMode(saved === "advanced" ? "advanced" : "simple", { silent: true });
    showSimpleStep(1, true);

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
