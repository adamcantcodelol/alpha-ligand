/* =====================================================================
 * core.js - the science.
 * A direct JavaScript port of pocket_profiler.py. Every number this file
 * produces is meant to match the Python script exactly.
 *
 * It is written so it works BOTH in the browser (loaded with a <script>
 * tag, which puts everything on window.PocketCore) AND in Node, so the
 * numbers can be tested against the Python version.
 * ===================================================================== */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.PocketCore = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /* -------------------------------------------------------------------
   * SECTION 1: Chemistry lookup tables.
   * Plain facts about the 20 amino acids. Nothing clever here.
   * ----------------------------------------------------------------- */

  const THREE_TO_ONE = {
    ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C",
    GLN: "Q", GLU: "E", GLY: "G", HIS: "H", ILE: "I",
    LEU: "L", LYS: "K", MET: "M", PHE: "F", PRO: "P",
    SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V",
  };

  // Kyte-Doolittle hydropathy: positive = water-hating, negative = water-loving.
  const KD = {
    A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5,
    Q: -3.5, E: -3.5, G: -0.4, H: -3.2, I: 4.5,
    L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6,
    S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2,
  };

  const HYDROPHOBIC = new Set("AVLIMFWCP");
  const AROMATIC = new Set("FWYH");
  const POSITIVE = new Set("KR");       // formally + at pH 7
  const NEGATIVE = new Set("DE");       // formally - at pH 7
  const POLAR_NEUTRAL = new Set("STNQCYH");

  // How many H-bond donors / acceptors each SIDE CHAIN offers.
  const DONORS = { R: 3, K: 1, W: 1, N: 1, Q: 1, S: 1, T: 1, Y: 1, H: 1, C: 1 };
  const ACCEPTORS = { D: 2, E: 2, N: 1, Q: 1, S: 1, T: 1, Y: 1, H: 1 };

  // The chemically important atoms of each side chain, used to place the
  // pharmacophore anchor points.
  const KEY_ATOMS = {
    R: ["CZ", "NH1", "NH2"], K: ["NZ"], D: ["CG", "OD1", "OD2"],
    E: ["CD", "OE1", "OE2"], F: ["CZ"], Y: ["OH", "CZ"],
    W: ["CZ2", "NE1"], H: ["NE2", "ND1"], S: ["OG"], T: ["OG1"],
    N: ["OD1", "ND2"], Q: ["OE1", "NE2"], C: ["SG"], M: ["SD"],
    L: ["CD1"], I: ["CD1"], V: ["CG1"], A: ["CB"], P: ["CG"], G: ["CA"],
  };

  const BACKBONE = new Set(["N", "CA", "C", "O", "OXT"]);

  // For each pocket residue type: what the LIGAND should bring to face it.
  // [feature type, suggested chemistry, why, base priority score]
  const COMPLEMENT_RULES = {
    R: ["anionic", "carboxylate / phosphate / sulfonate", "Arg guanidinium: salt bridge + bidentate H-bond", 10],
    K: ["anionic", "carboxylate / phosphate", "Lys ammonium: salt bridge", 10],
    H: ["anionic_or_aromatic", "carboxylate, or aryl ring for stacking", "His can be + at pH 7 and also stacks", 7],
    D: ["cationic", "amine / amidine / guanidine", "Asp carboxylate: salt bridge", 10],
    E: ["cationic", "amine / amidine / guanidine", "Glu carboxylate: salt bridge", 10],
    F: ["aromatic", "aryl / heteroaryl ring", "Phe ring: pi-pi stacking", 8],
    W: ["aromatic", "aryl ring or cation for cation-pi", "Trp indole: stacking / cation-pi", 9],
    Y: ["aromatic_and_hbond", "aryl ring plus an H-bond acceptor", "Tyr: stacking ring + donor OH", 9],
    S: ["hbond_acceptor", "carbonyl / ester / ether / N of heteroaryl", "Ser OH donates an H-bond", 5],
    T: ["hbond_acceptor", "carbonyl / ether oxygen", "Thr OH donates an H-bond", 5],
    N: ["hbond_mixed", "H-bond acceptor facing the NH2, donor facing the C=O", "Asn amide donates and accepts", 6],
    Q: ["hbond_mixed", "H-bond acceptor facing the NH2, donor facing the C=O", "Gln amide donates and accepts", 6],
    C: ["hbond_acceptor_or_covalent", "H-bond acceptor; or a mild electrophile", "Cys thiol - only consider covalent chemistry deliberately", 4],
    M: ["hydrophobic", "alkyl / thioether-friendly greasy group", "Met side chain is soft and greasy", 4],
    L: ["hydrophobic", "alkyl group", "Leu greasy wall", 4],
    I: ["hydrophobic", "alkyl group", "Ile greasy wall", 4],
    V: ["hydrophobic", "small alkyl group", "Val greasy wall", 3],
    A: ["hydrophobic", "methyl-sized contact", "Ala small greasy contact", 2],
    P: ["hydrophobic", "rigid hydrophobic contact", "Pro rigid greasy contact", 2],
    G: ["shape", "room to grow - Gly leaves space", "Gly has no side chain; the pocket is open here", 2],
  };

  const DISCLAIMER =
    "NOTE: This is a HYPOTHESIS-GENERATING pharmacophore derived purely from " +
    "the geometry and chemistry of the pocket residues. Nothing here has been " +
    "docked, scored, simulated or experimentally validated. Treat it as a " +
    "starting blueprint for designing or searching for a repressor ligand.";

  /* -------------------------------------------------------------------
   * SECTION 2: Small helpers
   * ----------------------------------------------------------------- */

  function round2(v) { return Math.round(v * 100) / 100; }
  function round3(v) { return Math.round(v * 1000) / 1000; }

  // The side chain is everything that is not backbone N/CA/C/O.
  // Glycine has no side chain, so we fall back to its CA.
  function sidechainAtoms(residue) {
    let atoms = residue.atoms.filter(function (a) { return !BACKBONE.has(a.name); });
    if (atoms.length === 0) atoms = residue.atoms.filter(function (a) { return a.name === "CA"; });
    return atoms;
  }

  // A readable label like "A:351".
  function residueKey(residue) { return residue.chain + ":" + residue.seq + (residue.icode || ""); }

  // Turn "A:351,A:383,A:524" (or "351,383" with a default chain) into pairs.
  function parseResidueSpec(spec, defaultChain) {
    const wanted = [];
    const tokens = String(spec).split(",");
    for (let i = 0; i < tokens.length; i++) {
      const token = tokens[i].trim();
      if (!token) continue;
      let chain, number;
      if (token.indexOf(":") !== -1) {
        const bits = token.split(":");
        chain = bits[0]; number = bits.slice(1).join(":");
      } else {
        chain = defaultChain; number = token;
      }
      if (chain === null || chain === undefined || chain === "") {
        throw new Error("Residue '" + token + "' has no chain. Write it as A:" + token + " instead.");
      }
      const n = parseInt(String(number).trim(), 10);
      if (isNaN(n)) throw new Error("Could not read a residue number from '" + token + "'.");
      wanted.push([chain.trim(), n]);
    }
    if (wanted.length === 0) throw new Error("No residues given.");
    return wanted;
  }

  /* -------------------------------------------------------------------
   * SECTION 3: Reading a structure file
   * Handles PDB format (ATOM/HETATM lines) and a basic mmCIF atom_site loop.
   * ----------------------------------------------------------------- */

  // Build the residue list from parsed atom records, keeping only the first
  // model and only one alternate location per atom (matching Biopython's
  // default of using a single conformation).
  function assembleResidues(atomRecords) {
    const byResidue = new Map();
    const order = [];
    for (let i = 0; i < atomRecords.length; i++) {
      const rec = atomRecords[i];
      const key = rec.chain + "|" + rec.seq + "|" + (rec.icode || "");
      let res = byResidue.get(key);
      if (!res) {
        res = {
          chain: rec.chain, seq: rec.seq, icode: rec.icode || "",
          resname: rec.resname, atoms: [], seen: new Set(),
        };
        byResidue.set(key, res);
        order.push(res);
      }
      // Skip a second alternate location for an atom name we already have.
      if (res.seen.has(rec.name)) continue;
      res.seen.add(rec.name);
      res.atoms.push({ name: rec.name, xyz: rec.xyz });
    }
    order.forEach(function (r) {
      r.letter = THREE_TO_ONE[r.resname.toUpperCase()] || null;
      delete r.seen;
    });
    return order;
  }

  function parsePDB(text) {
    const records = [];
    const lines = text.split("\n");
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];
      // Stop at the end of the first model (NMR files hold several).
      if (line.substring(0, 6) === "ENDMDL") break;
      const tag = line.substring(0, 6);
      if (tag !== "ATOM  " && tag !== "HETATM") continue;
      const altloc = line.charAt(16);
      if (altloc !== " " && altloc !== "A") continue;
      const name = line.substring(12, 16).trim();
      const resname = line.substring(17, 20).trim();
      const chain = line.substring(21, 22).trim() || " ";
      const seq = parseInt(line.substring(22, 26).trim(), 10);
      const icode = line.substring(26, 27).trim();
      const x = parseFloat(line.substring(30, 38));
      const y = parseFloat(line.substring(38, 46));
      const z = parseFloat(line.substring(46, 54));
      if (isNaN(seq) || isNaN(x) || isNaN(y) || isNaN(z)) continue;
      records.push({ name: name, resname: resname, chain: chain, seq: seq, icode: icode, xyz: [x, y, z] });
    }
    if (records.length === 0) throw new Error("No ATOM or HETATM records found - is this really a PDB file?");
    return assembleResidues(records);
  }

  // A deliberately simple mmCIF reader: it only understands the atom_site
  // loop, which is all we need.
  function parseMMCIF(text) {
    const lines = text.split("\n");
    const records = [];
    let cols = null, colIndex = {};
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];
      const t = line.trim();
      if (t.indexOf("_atom_site.") === 0) {
        if (!cols) { cols = []; colIndex = {}; }
        const field = t.split(/\s+/)[0].substring("_atom_site.".length);
        colIndex[field] = cols.length;
        cols.push(field);
        continue;
      }
      if (!cols) continue;
      if (t === "" || t === "#" || t.indexOf("loop_") === 0) { if (records.length) break; else continue; }
      const parts = t.split(/\s+/);
      if (parts.length < cols.length) continue;
      const get = function (f) { const idx = colIndex[f]; return idx === undefined ? null : parts[idx]; };
      const group = get("group_PDB");
      if (group !== "ATOM" && group !== "HETATM") continue;
      const model = get("pdbx_PDB_model_num");
      if (model !== null && parseInt(model, 10) !== 1) continue;
      const alt = get("label_alt_id");
      if (alt && alt !== "." && alt !== "?" && alt !== "A") continue;
      let chain = get("auth_asym_id"); if (!chain || chain === ".") chain = get("label_asym_id");
      let seqRaw = get("auth_seq_id"); if (!seqRaw || seqRaw === ".") seqRaw = get("label_seq_id");
      let resname = get("auth_comp_id"); if (!resname || resname === ".") resname = get("label_comp_id");
      let name = get("auth_atom_id"); if (!name || name === ".") name = get("label_atom_id");
      let icode = get("pdbx_PDB_ins_code"); if (!icode || icode === "." || icode === "?") icode = "";
      const seq = parseInt(seqRaw, 10);
      const x = parseFloat(get("Cartn_x")), y = parseFloat(get("Cartn_y")), z = parseFloat(get("Cartn_z"));
      if (isNaN(seq) || isNaN(x) || isNaN(y) || isNaN(z)) continue;
      records.push({
        name: String(name).replace(/"/g, ""), resname: String(resname), chain: String(chain),
        seq: seq, icode: icode, xyz: [x, y, z],
      });
    }
    if (records.length === 0) throw new Error("No atom_site records found - is this really an mmCIF file?");
    return assembleResidues(records);
  }

  function parseStructure(text, filename) {
    const looksCif = /\.(cif|mmcif)$/i.test(filename || "") || text.indexOf("_atom_site.") !== -1;
    return looksCif ? parseMMCIF(text) : parsePDB(text);
  }

  function chainsPresent(residues) {
    const seen = [];
    residues.forEach(function (r) { if (seen.indexOf(r.chain) === -1) seen.push(r.chain); });
    return seen;
  }

  /* -------------------------------------------------------------------
   * SECTION 4: Selecting the active site and the pocket shell
   * ----------------------------------------------------------------- */

  function collectActiveSite(residues, wanted) {
    const chains = chainsPresent(residues);
    const found = [];
    for (let i = 0; i < wanted.length; i++) {
      const chainId = wanted[i][0], number = wanted[i][1];
      if (chains.indexOf(chainId) === -1) {
        throw new Error("Chain '" + chainId + "' is not in this structure. Chains present: " + chains.join(", "));
      }
      const hit = residues.find(function (r) { return r.chain === chainId && r.seq === number; });
      if (!hit) throw new Error("Residue " + number + " was not found in chain " + chainId + ".");
      if (!hit.letter) {
        throw new Error("Residue " + chainId + ":" + number + " is '" + hit.resname + "', not a standard amino acid.");
      }
      found.push(hit);
    }
    return found;
  }

  // Every standard amino acid with any atom within `cutoff` angstroms of an
  // active-site SIDE-CHAIN atom. That set is the wall of the pocket.
  function collectPocketShell(residues, activeSite, cutoff) {
    const probes = [];
    activeSite.forEach(function (res) {
      sidechainAtoms(res).forEach(function (a) { probes.push(a.xyz); });
    });
    const c2 = cutoff * cutoff;
    const shell = [];
    for (let i = 0; i < residues.length; i++) {
      const res = residues[i];
      if (!res.letter) continue;              // skip water, ions, bound ligands
      let near = false;
      for (let j = 0; j < res.atoms.length && !near; j++) {
        const a = res.atoms[j].xyz;
        for (let k = 0; k < probes.length; k++) {
          const p = probes[k];
          const dx = a[0] - p[0], dy = a[1] - p[1], dz = a[2] - p[2];
          if (dx * dx + dy * dy + dz * dz <= c2) { near = true; break; }
        }
      }
      if (near) shell.push(res);
    }
    return shell;
  }

  /* -------------------------------------------------------------------
   * SECTION 5: Profiling the pocket chemistry and geometry
   * ----------------------------------------------------------------- */

  function profilePocket(shell) {
    const letters = shell.map(function (r) { return r.letter; });
    const n = letters.length;
    function frac(group) {
      const hits = letters.filter(function (c) { return group.has(c); }).length;
      return { count: hits, fraction: n ? round3(hits / n) : 0 };
    }
    const hydrophobic = frac(HYDROPHOBIC), aromatic = frac(AROMATIC);
    const positive = frac(POSITIVE), negative = frac(NEGATIVE), polar = frac(POLAR_NEUTRAL);

    // Composition, most common first (ties keep first-seen order).
    const counts = new Map();
    letters.forEach(function (c) { counts.set(c, (counts.get(c) || 0) + 1); });
    const composition = Array.from(counts.entries()).sort(function (a, b) { return b[1] - a[1]; });

    let kdSum = 0;
    letters.forEach(function (c) { kdSum += KD[c]; });

    let donors = 0, acceptors = 0;
    letters.forEach(function (c) { donors += DONORS[c] || 0; acceptors += ACCEPTORS[c] || 0; });

    return {
      n_residues: n,
      composition: composition,
      hydrophobic: hydrophobic, aromatic: aromatic,
      positive: positive, negative: negative, polar_neutral: polar,
      mean_kd_hydropathy: n ? round3(kdSum / n) : 0,
      net_formal_charge_pH7: positive.count - negative.count,
      sidechain_hbond_donors: donors,
      sidechain_hbond_acceptors: acceptors,
    };
  }

  function pocketGeometry(shell) {
    const coords = [];
    shell.forEach(function (r) { r.atoms.forEach(function (a) { coords.push(a.xyz); }); });
    if (coords.length === 0) return {};

    let sx = 0, sy = 0, sz = 0;
    let minx = Infinity, miny = Infinity, minz = Infinity;
    let maxx = -Infinity, maxy = -Infinity, maxz = -Infinity;
    coords.forEach(function (c) {
      sx += c[0]; sy += c[1]; sz += c[2];
      if (c[0] < minx) minx = c[0]; if (c[0] > maxx) maxx = c[0];
      if (c[1] < miny) miny = c[1]; if (c[1] > maxy) maxy = c[1];
      if (c[2] < minz) minz = c[2]; if (c[2] > maxz) maxz = c[2];
    });
    const n = coords.length;

    // Longest span between any two pocket atoms. Big pockets are
    // subsampled so this stays fast, exactly as the Python does.
    let sample = coords;
    if (n > 400) {
      const step = Math.max(1, Math.floor(n / 400));
      sample = [];
      for (let i = 0; i < n; i += step) sample.push(coords[i]);
    }
    let span = 0;
    for (let i = 0; i < sample.length; i++) {
      for (let j = i + 1; j < sample.length; j++) {
        const a = sample[i], b = sample[j];
        const dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (d > span) span = d;
      }
    }

    return {
      centroid: [round2(sx / n), round2(sy / n), round2(sz / n)],
      bounding_box_xyz: [round2(maxx - minx), round2(maxy - miny), round2(maxz - minz)],
      max_atom_span_A: round2(span),
      n_pocket_atoms: n,
    };
  }

  /* -------------------------------------------------------------------
   * SECTION 6: Working backwards -> the complementary pharmacophore
   * ----------------------------------------------------------------- */

  // A representative 3D point for a side chain: the average position of its
  // chemically important atoms.
  function anchorPoint(residue) {
    const names = KEY_ATOMS[residue.letter] || [];
    let atoms = residue.atoms.filter(function (a) { return names.indexOf(a.name) !== -1; });
    if (atoms.length === 0) atoms = sidechainAtoms(residue);
    const n = atoms.length;
    const out = [0, 0, 0];
    atoms.forEach(function (a) { out[0] += a.xyz[0]; out[1] += a.xyz[1]; out[2] += a.xyz[2]; });
    return [round2(out[0] / n), round2(out[1] / n), round2(out[2] / n)];
  }

  function buildPharmacophore(activeSite, shell) {
    const activeKeys = new Set(activeSite.map(residueKey));
    const features = [];
    shell.forEach(function (res) {
      const rule = COMPLEMENT_RULES[res.letter];
      if (!rule) return;
      const key = residueKey(res);
      const isNamed = activeKeys.has(key);
      // The residues the user named matter most, so their score is boosted.
      const score = rule[3] * (isNamed ? 1.5 : 1.0);
      features.push({
        pocket_residue: res.letter + key,
        pocket_residue_name: res.resname,
        is_named_active_site: isNamed,
        ligand_feature_type: rule[0],
        suggested_ligand_chemistry: rule[1],
        rationale: rule[2],
        approx_xyz: anchorPoint(res),
        priority_score: round2(score),
      });
    });
    features.sort(function (a, b) { return b.priority_score - a.priority_score; });
    return features;
  }

  // Rough, honest ligand property windows from pocket size and chemistry.
  // These are rules of thumb, not measurements.
  function suggestProperties(profile, geometry) {
    const span = geometry.max_atom_span_A === undefined ? 15 : geometry.max_atom_span_A;
    const targetLen = Math.max(6.0, span * 0.55);          // ligand spans ~half the pocket
    let mwLow = Math.round(Math.max(150, targetLen * 16) / 10) * 10;   // ~16-30 Da per angstrom
    let mwHigh = Math.round(Math.min(650, targetLen * 30) / 10) * 10;
    if (mwHigh <= mwLow) mwHigh = mwLow + 150;

    const kd = profile.mean_kd_hydropathy;
    const logpMid = 1.5 + kd;                               // greasy pocket -> higher logP
    const logpLow = Math.round(Math.max(-1.0, logpMid - 1.5) * 10) / 10;
    const logpHigh = Math.round(Math.min(5.0, logpMid + 1.5) * 10) / 10;

    // The ligand charge should oppose the pocket's, but stay modest.
    const net = profile.net_formal_charge_pH7;
    const floorDiv2 = Math.floor(-net / 2);                 // Python's -net // 2
    const ligandCharge = net ? Math.max(-3, Math.min(3, floorDiv2)) : 0;

    return {
      molecular_weight_Da: [mwLow, mwHigh],
      clogp: [logpLow, logpHigh],
      rotatable_bonds_max: Math.max(3, Math.min(10, Math.trunc(targetLen / 2.5))),
      hbond_donors: [0, Math.max(1, Math.min(5, Math.floor(profile.sidechain_hbond_acceptors / 3)))],
      hbond_acceptors: [1, Math.max(2, Math.min(10, Math.floor(profile.sidechain_hbond_donors / 2)))],
      suggested_net_charge: ligandCharge,
      basis: "Derived from pocket span " + span.toFixed(1) + " A, mean Kyte-Doolittle " +
        kd.toFixed(2) + " and net pocket charge " + (net >= 0 ? "+" : "") + net +
        ". Rules of thumb only.",
    };
  }

  /* -------------------------------------------------------------------
   * SECTION 7: Run everything, and render the plain-text report
   * ----------------------------------------------------------------- */

  function analyse(residues, residueSpec, cutoff, defaultChain) {
    const wanted = parseResidueSpec(residueSpec, defaultChain || null);
    const activeSite = collectActiveSite(residues, wanted);
    const shell = collectPocketShell(residues, activeSite, cutoff);
    if (shell.length === 0) throw new Error("No residues found within the cutoff. Try a larger cutoff.");
    const profile = profilePocket(shell);
    const geometry = pocketGeometry(shell);
    const features = buildPharmacophore(activeSite, shell);
    const properties = suggestProperties(profile, geometry);
    return {
      activeSite: activeSite, shell: shell, profile: profile,
      geometry: geometry, features: features, properties: properties,
      cutoff: cutoff,
    };
  }

  function pad(s, width) { s = String(s); while (s.length < width) s += " "; return s; }
  function signed(n, places) {
    const v = Number(n).toFixed(places === undefined ? 0 : places);
    return (Number(n) >= 0 ? "+" : "") + v;
  }

  function renderReport(name, result, top) {
    top = top || 25;
    const p = result.profile, g = result.geometry, props = result.properties;
    const L = [];
    const bar = new Array(73).join("=");
    L.push(bar);
    L.push("POCKET PROFILE AND COMPLEMENTARY LIGAND BLUEPRINT");
    L.push("Structure : " + name);
    L.push("Cutoff    : " + result.cutoff.toFixed(1) + " A around active-site side chains");
    L.push(bar);
    L.push("");
    L.push("1. ACTIVE SITE RESIDUES YOU NAMED");
    result.activeSite.forEach(function (r) {
      L.push("   " + r.letter + residueKey(r) + " (" + r.resname + ")");
    });
    L.push("");
    L.push("2. POCKET SHELL (" + p.n_residues + " residues within " + result.cutoff.toFixed(1) + " A)");
    L.push("   " + result.shell.map(residueKey).sort().join(", "));
    L.push("");
    L.push("3. CHEMICAL CHARACTER OF THE POCKET");
    L.push("   Composition          : " + p.composition.map(function (kv) { return kv[0] + "=" + kv[1]; }).join(", "));
    [["Hydrophobic", "hydrophobic"], ["Aromatic", "aromatic"], ["Positive (K/R)", "positive"],
     ["Negative (D/E)", "negative"], ["Polar neutral", "polar_neutral"]].forEach(function (pair) {
      const d = p[pair[1]];
      L.push("   " + pad(pair[0], 21) + ": " + d.count + "  (" + Math.round(100 * d.fraction) + "%)");
    });
    L.push("   Mean KD hydropathy   : " + signed(p.mean_kd_hydropathy, 2) +
      "  (" + (p.mean_kd_hydropathy > 0 ? "greasy pocket" : "polar pocket") + ")");
    L.push("   Net formal charge    : " + signed(p.net_formal_charge_pH7) + " at pH 7");
    L.push("   Side-chain H-bond    : " + p.sidechain_hbond_donors + " donors, " + p.sidechain_hbond_acceptors + " acceptors");
    L.push("");
    L.push("4. POCKET SIZE AND SHAPE");
    L.push("   Centroid (x,y,z)     : [" + g.centroid.join(", ") + "]");
    L.push("   Bounding box (A)     : [" + g.bounding_box_xyz.join(", ") + "]");
    L.push("   Max atom span (A)    : " + g.max_atom_span_A);
    L.push("   Pocket atoms counted : " + g.n_pocket_atoms);
    L.push("   SASA                 : not computed in the browser version");
    L.push("");
    L.push("5. WORKING BACKWARDS: COMPLEMENTARY LIGAND FEATURES (ranked)");
    L.push("   " + pad("POCKET RES", 12) + " " + pad("LIGAND FEATURE", 22) + " " + pad("SUGGESTED CHEMISTRY", 34) + " APPROX XYZ");
    result.features.slice(0, top).forEach(function (f) {
      L.push("  " + (f.is_named_active_site ? "*" : " ") + pad(f.pocket_residue, 12) + " " +
        pad(f.ligand_feature_type, 22) + " " + pad(f.suggested_ligand_chemistry.substring(0, 34), 34) +
        " [" + f.approx_xyz.join(", ") + "]");
    });
    if (result.features.length > top) {
      L.push("   ... " + (result.features.length - top) + " more features in the JSON file.");
    }
    L.push("   (* = residue you named as active site; these are weighted highest.)");
    L.push("");
    L.push("6. SUGGESTED LIGAND PROPERTY WINDOWS");
    L.push("   Molecular weight (Da): " + props.molecular_weight_Da[0] + " - " + props.molecular_weight_Da[1]);
    L.push("   cLogP                : " + props.clogp[0] + " - " + props.clogp[1]);
    L.push("   Rotatable bonds      : <= " + props.rotatable_bonds_max);
    L.push("   H-bond donors        : " + props.hbond_donors[0] + " - " + props.hbond_donors[1]);
    L.push("   H-bond acceptors     : " + props.hbond_acceptors[0] + " - " + props.hbond_acceptors[1]);
    L.push("   Net charge           : " + signed(props.suggested_net_charge));
    L.push("   Basis                : " + props.basis);
    L.push("");
    L.push("7. HOW TO READ THIS");
    L.push("   Each line in section 5 says: 'the pocket offers X here, so your");
    L.push("   ligand should present the opposite of X at roughly these");
    L.push("   coordinates'. Build or search for a molecule that places those");
    L.push("   groups near those points, then dock and test it.");
    L.push("");
    L.push(DISCLAIMER);
    return L.join("\n");
  }

  function toJSON(name, result) {
    return {
      structure: name,
      cutoff_A: result.cutoff,
      active_site_residues: result.activeSite.map(function (r) { return r.letter + residueKey(r); }),
      pocket_shell_residues: result.shell.map(residueKey).sort(),
      profile: result.profile,
      geometry: result.geometry,
      solvent_accessibility: { available: false, note: "SASA is not computed in the browser version." },
      pharmacophore_features: result.features,
      suggested_ligand_properties: result.properties,
      disclaimer: DISCLAIMER,
    };
  }

  /* -------------------------------------------------------------------
   * SECTION 8: Finding a likely pocket automatically.
   * Same idea as the notebook's residue finder: take the biggest
   * non-water, non-buffer hetero group in the file (usually the bound
   * drug or substrate) and report the amino acids around it.
   * ----------------------------------------------------------------- */

  // Things that are in crystal structures but are not the interesting
  // ligand: solvent, cryoprotectants, buffers, ions, sugars.
  const IGNORED_HET = new Set([
    "HOH", "WAT", "DOD", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "MPD",
    "SO4", "PO4", "ACT", "ACY", "FMT", "CIT", "TRS", "MES", "EPE", "IMD", "DMS",
    "BME", "NO3", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "NI",
    "CD", "HG", "IOD", "BR", "F", "CO", "CS", "SR", "BA", "UNX", "UNL",
  ]);

  // The biggest hetero residue worth caring about, or null.
  function findBoundLigand(residues) {
    let best = null;
    residues.forEach(function (r) {
      if (r.letter) return;                              // standard amino acid
      const name = (r.resname || "").toUpperCase();
      if (IGNORED_HET.has(name)) return;
      if (r.atoms.length < 6) return;                    // ions and fragments
      if (!best || r.atoms.length > best.atoms.length) best = r;
    });
    return best;
  }

  // Amino acids with any atom within `radius` of the ligand, nearest first.
  function residuesNearLigand(residues, ligand, radius) {
    radius = radius || 5.0;
    const r2 = radius * radius;
    const near = [];
    residues.forEach(function (res) {
      if (!res.letter) return;
      let bestD2 = Infinity;
      for (let i = 0; i < res.atoms.length; i++) {
        const a = res.atoms[i].xyz;
        for (let j = 0; j < ligand.atoms.length; j++) {
          const b = ligand.atoms[j].xyz;
          const dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
          const d2 = dx * dx + dy * dy + dz * dz;
          if (d2 < bestD2) bestD2 = d2;
        }
      }
      if (bestD2 <= r2) near.push({ residue: res, dist: round2(Math.sqrt(bestD2)) });
    });
    near.sort(function (a, b) { return a.dist - b.dist; });
    return near;
  }

  // One call for the UI: {ligand, near, keys} or null if nothing is bound.
  function findLikelyPocket(residues, radius, maxResidues) {
    const ligand = findBoundLigand(residues);
    if (!ligand) return null;
    const near = residuesNearLigand(residues, ligand, radius || 5.0);
    if (near.length === 0) return null;
    const chosen = near.slice(0, maxResidues || 8);
    return {
      ligand: {
        resname: ligand.resname, chain: ligand.chain, seq: ligand.seq,
        n_atoms: ligand.atoms.length, key: residueKey(ligand),
      },
      radius: radius || 5.0,
      near: near,
      keys: chosen.map(function (h) { return residueKey(h.residue); }),
      labels: chosen.map(function (h) { return h.residue.letter + residueKey(h.residue); }),
      // Same residues, written for someone who has never read a PDB file.
      friendly: chosen.map(function (h) {
        return h.residue.resname + " " + h.residue.seq + " (" + residueKey(h.residue) + ")";
      }),
    };
  }

  /* -------------------------------------------------------------------
   * SECTION 9: Saying all of this in plain English (Simple mode).
   * No new science - it only re-words what analyse() already produced.
   * ----------------------------------------------------------------- */

  // What each pharmacophore feature type means, for someone who has never
  // taken a chemistry class.
  const PLAIN_FEATURE = {
    anionic: "a negatively charged part (the kind of group vinegar has)",
    cationic: "a positively charged part (the kind of group ammonia has)",
    anionic_or_aromatic: "either a negatively charged part or a flat ring",
    aromatic: "a flat ring of carbon atoms, like the ring in benzene",
    aromatic_and_hbond: "a flat ring that also has an oxygen or nitrogen on it",
    hydrophobic: "a greasy, oily part made of plain carbon",
    hbond_acceptor: "an oxygen or nitrogen that can grab onto a hydrogen",
    hbond_mixed: "oxygens and nitrogens that can both give and take hydrogens",
    hbond_acceptor_or_covalent: "an oxygen or nitrogen that can grab onto a hydrogen",
    shape: "a little extra bulk - there is spare room here",
  };

  // Everyday molecules to compare a weight against.
  const SIZE_REFS = [
    [46, "a molecule of alcohol"], [151, "a paracetamol tablet's active ingredient"],
    [180, "an aspirin molecule"], [194, "a caffeine molecule"],
    [206, "an ibuprofen molecule"], [300, "a typical small drug"],
    [388, "a tamoxifen molecule"], [480, "a large drug molecule"],
    [734, "a very large drug molecule"],
  ];

  function sizeComparison(mw) {
    let best = SIZE_REFS[0];
    for (let i = 0; i < SIZE_REFS.length; i++) {
      if (Math.abs(SIZE_REFS[i][0] - mw) < Math.abs(best[0] - mw)) best = SIZE_REFS[i];
    }
    return best[1];
  }

  // Turn an analyse() result into a handful of plain sentences.
  function plainSummary(name, result) {
    const p = result.profile, g = result.geometry, props = result.properties;
    const sentences = [];

    // 1. How big the pocket is and what it is made of.
    const greasy = Math.round(100 * p.hydrophobic.fraction);
    const rings = Math.round(100 * p.aromatic.fraction);
    const charged = Math.round(100 * (p.positive.fraction + p.negative.fraction));
    sentences.push(
      "The pocket you picked in " + name + " is lined by " + p.n_residues +
      " building blocks of the protein, and the wall of it measures about " +
      Math.round(g.max_atom_span_A) + " angstroms end to end (an angstrom is a ten-millionth of a " +
      "millimetre). The molecule that binds only has to fill the dent in the middle, so it is much " +
      "smaller than that.");

    let character;
    if (p.mean_kd_hydropathy > 0.5) character = "greasy and water-hating";
    else if (p.mean_kd_hydropathy > -0.5) character = "a mix of greasy and water-friendly";
    else character = "water-friendly rather than greasy";
    sentences.push(
      "Overall it is " + character + ". " + greasy + " out of every 100 of those building blocks are " +
      "oily ones, " + rings + " have flat rings, and " + charged + " carry an electrical charge.");

    // 2. The charge of the pocket, and so of the molecule that fits it.
    const net = p.net_formal_charge_pH7;
    if (net > 0) {
      sentences.push("The pocket is positively charged overall (" + (net > 0 ? "+" : "") + net +
        "), so a molecule that fits well should carry some negative charge - opposite charges pull together.");
    } else if (net < 0) {
      sentences.push("The pocket is negatively charged overall (" + net +
        "), so a molecule that fits well should carry some positive charge - opposite charges pull together.");
    } else {
      sentences.push("The pocket has no overall electrical charge, so the molecule that fits it " +
        "does not need one either.");
    }

    // 3. What chemistry the molecule should carry: top feature types by
    //    total priority, worded plainly.
    const weight = new Map();
    result.features.forEach(function (f) {
      weight.set(f.ligand_feature_type,
        (weight.get(f.ligand_feature_type) || 0) + f.priority_score);
    });
    const ranked = Array.from(weight.entries()).sort(function (a, b) { return b[1] - a[1]; });
    const wants = ranked.slice(0, 3).map(function (kv) {
      return PLAIN_FEATURE[kv[0]] || kv[0].replace(/_/g, " ");
    });
    let wantsText;
    if (wants.length === 1) wantsText = wants[0];
    else if (wants.length === 2) wantsText = wants[0] + ", and " + wants[1];
    else wantsText = wants[0] + ", " + wants[1] + ", and " + wants[2];
    sentences.push("So the molecule you are looking for should have " + wantsText + ".");

    // 4. Size window, in everyday terms.
    const mwLow = props.molecular_weight_Da[0], mwHigh = props.molecular_weight_Da[1];
    const mid = Math.round((mwLow + mwHigh) / 2);
    sentences.push("It should weigh somewhere between " + mwLow + " and " + mwHigh +
      " daltons - about the size of " + sizeComparison(mid) + ". Smaller than that rattles around " +
      "in the pocket; bigger than that will not fit in.");

    // The honesty line, in beginner words.
    const caveat =
      "Important: this is a shopping list, not a prediction. It says what kind of molecule the " +
      "pocket is shaped and charged for - nothing here has been tested, in a computer or in a lab, " +
      "and most molecules that look right on paper do not actually work. Treat it as a short list " +
      "worth trying at the bench, never as proof that something binds.";

    // The headline tags for the little bar chart.
    const bars = [
      { label: "Greasy", pct: greasy, hint: "oily, water-hating building blocks" },
      { label: "Flat rings", pct: rings, hint: "ring-shaped building blocks a molecule can stack against" },
      { label: "Positive", pct: Math.round(100 * p.positive.fraction), hint: "positively charged building blocks" },
      { label: "Negative", pct: Math.round(100 * p.negative.fraction), hint: "negatively charged building blocks" },
      { label: "Water-friendly", pct: Math.round(100 * p.polar_neutral.fraction), hint: "polar but uncharged building blocks" },
    ];

    return { sentences: sentences, caveat: caveat, bars: bars, wants: wants };
  }

  return {
    parseStructure: parseStructure, parsePDB: parsePDB, parseMMCIF: parseMMCIF,
    chainsPresent: chainsPresent, parseResidueSpec: parseResidueSpec,
    analyse: analyse, renderReport: renderReport, toJSON: toJSON,
    residueKey: residueKey, THREE_TO_ONE: THREE_TO_ONE, DISCLAIMER: DISCLAIMER,
    findBoundLigand: findBoundLigand, residuesNearLigand: residuesNearLigand,
    findLikelyPocket: findLikelyPocket, plainSummary: plainSummary,
  };
});
