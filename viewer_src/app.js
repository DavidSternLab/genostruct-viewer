/* =====================================================================
   Integrated genome / gene / protein / structure viewer  — front-end
   Three color-linked panels:
     (1) Mol* 3D structure        (top)
     (2) protein sequence track   (middle)
     (3) genome track (exons/CDS/introns + structural elements) (bottom)
   Each structural element (helix/strand) carries a unique color that is
   identical across all three panels. Hover/click on any panel highlights
   the same element everywhere.
   ===================================================================== */
(function () {
"use strict";

var STATE = {
  index: null, rec: null, plugin: null,
  structureRef: null,
  colorThemeName: "genostruct-elements",
  selected: null,
};

/* ---- embedded-asset decoding (gzip + base64, inflated with pako) ------- */
function b64ToU8(b64) {
  var bin = atob(b64), n = bin.length, u = new Uint8Array(n);
  for (var i = 0; i < n; i++) u[i] = bin.charCodeAt(i);
  return u;
}
function inflateStr(b64) { return pako.inflate(b64ToU8(b64), { to: "string" }); }
function decodeAsset(id) {
  var el = document.getElementById(id);
  return el ? inflateStr(el.textContent.trim()) : null;
}
function loadRecord(tid) { return JSON.parse(inflateStr(EMBED.data[tid])); }
function loadPDB(tid) { return inflateStr(EMBED.pdb[tid]); }

/* ---- color helpers ---------------------------------------------------- */
function hexToInt(hex) { return parseInt(hex.replace('#', ''), 16); }

/* =====================================================================
   MOL*  — structure panel + custom per-residue color theme
   ===================================================================== */
function registerElementTheme(plugin) {
  // Register exactly once. Mol*'s registry.add() THROWS if the name already
  // exists, and registry.has() expects a PROVIDER OBJECT (it reads e.name), so
  // has("string") is always false — using it as the guard would let add() run on
  // every transcript switch and throw on the 2nd, aborting the structure load.
  // We therefore track registration with our own flag. The provider reads
  // STATE.rec dynamically, so a single registration serves every transcript.
  if (STATE._themeRegistered) return;
  var reg = plugin.representation.structure.themes.colorThemeRegistry;
  var DEFAULT = 0xBFBFBF; // coil grey
  var provider = {
    name: STATE.colorThemeName,
    label: "Structural elements",
    category: "Custom",
    factory: function (ctx, props) {
      function colorFor(location) {
        var recColors = STATE.rec ? STATE.rec.residue_color : null;
        if (!recColors) return DEFAULT;
        var off = (STATE.rec.model && STATE.rec.model.offset) || 0;
        var unit = location.unit, el = location.element;
        if (!unit || el === undefined || el === null) return DEFAULT;
        var residueIndex;
        try { residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index[el]; }
        catch (e) { residueIndex = undefined; }
        if (residueIndex === undefined || residueIndex === null) return DEFAULT;
        // residueIndex is 0-based within model; map to full protein residue.
        var protRes = off + residueIndex;
        if (protRes >= 0 && protRes < recColors.length) {
          var c = recColors[protRes];
          if (c) return hexToInt(c);
        }
        return DEFAULT;
      }
      return {
        factory: provider.factory,
        granularity: "group",
        color: colorFor,
        props: props || {},
        description: "Genostruct structural elements",
      };
    },
    getParams: function () { return {}; },
    defaultValues: {},
    isApplicable: function () { return true; },
  };
  try {
    reg.add(provider);
    STATE._themeRegistered = true;
  } catch (e) {
    // Already registered from a prior call (e.g. hot reload) — treat as success
    // so we never retry and never let the throw propagate into loadStructure.
    STATE._themeRegistered = true;
  }
}

async function loadStructure(tid) {
  var plugin = STATE.plugin;
  var status = document.getElementById("structStatus");
  try {
    if (status) status.textContent = "";
    // register the theme BEFORE clearing so a theme error can't strand us with
    // an empty viewport; registration is idempotent (guarded).
    registerElementTheme(plugin);
    await plugin.clear();
    var pdbText = loadPDB(tid);
    if (!pdbText) throw new Error("no embedded PDB for " + tid);
    var data = await plugin.builders.data.rawData({ data: pdbText, label: tid });
    var traj = await plugin.builders.structure.parseTrajectory(data, "pdb");
    var model = await plugin.builders.structure.createModel(traj);
    var structure = await plugin.builders.structure.createStructure(model);
    // Try our custom element coloring; if the theme ever fails, fall back to a
    // built-in theme so the structure is ALWAYS visible.
    try {
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: "cartoon", color: STATE.colorThemeName,
      });
    } catch (themeErr) {
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: "cartoon", color: "chain-id",
      });
      if (status) status.textContent = "(structure shown with default coloring — custom theme unavailable)";
    }
    STATE.structureRef = structure;
    plugin.managers.camera.reset();
  } catch (err) {
    if (status) status.textContent = "Could not display structure for " + tid + ": " + (err && err.message ? err.message : err);
    STATE.structureRef = null;
    // rethrow-free: keep the app responsive so other panels/transcripts work
    if (window.console) console.error("loadStructure failed", err);
  }
}

/* =====================================================================
   SEQUENCE PANEL
   ===================================================================== */
function renderSequence() {
  var rec = STATE.rec;
  var host = document.getElementById("seqTrack");
  host.innerHTML = "";
  var seq = rec.protein_sequence, colors = rec.residue_color, elems = rec.residue_element;
  var wrap = document.createElement("div"); wrap.className = "seqwrap";
  if (!STATE.seqZoom) STATE.seqZoom = 1;
  wrap.style.fontSize = (13 * STATE.seqZoom) + "px";
  wrap.style.letterSpacing = (1 * STATE.seqZoom) + "px";
  for (var i = 0; i < seq.length; i++) {
    var span = document.createElement("span");
    span.className = "aa";
    span.textContent = seq[i];
    span.dataset.res = (i + 1);
    if (colors[i]) {
      span.style.background = colors[i];
      span.style.color = "#fff";
      span.dataset.elem = elems[i];
      span.classList.add("in-elem");
    }
    (function (span) {
      span.addEventListener("mouseenter", function () {
        if (span.dataset.elem) highlightElement(span.dataset.elem, true);
        showResidueTip(span);
      });
      span.addEventListener("mouseleave", function () {
        if (span.dataset.elem) highlightElement(span.dataset.elem, false);
        hideTip();
      });
      span.addEventListener("click", function () {
        if (span.dataset.elem) selectElement(span.dataset.elem);
      });
    })(span);
    wrap.appendChild(span);
  }
  host.appendChild(wrap);
}

/* =====================================================================
   GENOME PANEL (SVG)
   ===================================================================== */
function genomeExtent() {
  var exons = STATE.rec.exons.slice().sort(function (a, b) { return a[0] - b[0]; });
  return [exons[0][0], exons[exons.length - 1][1]];
}
function resetGenomeView() {
  var ext = genomeExtent();
  STATE.gview = { start: ext[0], end: ext[1] };
  renderGenome();
}
function renderGenome() {
  var rec = STATE.rec;
  var host = document.getElementById("genomeTrack");
  host.innerHTML = "";
  var exons = rec.exons.slice().sort(function (a, b) { return a[0] - b[0]; });
  var cds = rec.cds.slice().sort(function (a, b) { return a[0] - b[0]; });
  if (!exons.length) { host.textContent = "No exon features."; return; }
  var ext = genomeExtent();
  if (!STATE.gview) STATE.gview = { start: ext[0], end: ext[1] };
  // clamp view to the locus
  var vs = Math.max(ext[0], STATE.gview.start), ve = Math.min(ext[1], STATE.gview.end);
  if (ve - vs < 3) { ve = vs + 3; }
  STATE.gview = { start: vs, end: ve };
  var W = host.clientWidth || 1000, PAD = 50, H = 150;
  var innerW = W - 2 * PAD, vspan = (ve - vs) || 1;
  function X(g) { return PAD + (g - vs) / vspan * innerW; }
  function clampX(x) { return Math.max(PAD, Math.min(W - PAD, x)); }
  var NS = "http://www.w3.org/2000/svg";
  var svg = document.createElementNS(NS, "svg");
  svg.setAttribute("width", W); svg.setAttribute("height", H);
  svg.setAttribute("class", "genome-svg");
  svg.setAttribute("viewBox", "0 0 " + W + " " + H);
  var yLine = 78, yExon = 68, exonH = 20, yCDS = 100, cdsH = 28;

  // clip rect so features don't spill into the padding when zoomed
  var defs = document.createElementNS(NS, "defs");
  var clip = document.createElementNS(NS, "clipPath"); clip.setAttribute("id", "plotClip");
  var cr = document.createElementNS(NS, "rect");
  cr.setAttribute("x", PAD); cr.setAttribute("y", 0); cr.setAttribute("width", innerW); cr.setAttribute("height", H);
  clip.appendChild(cr); defs.appendChild(clip); svg.appendChild(defs);
  var plot = document.createElementNS(NS, "g"); plot.setAttribute("clip-path", "url(#plotClip)");
  svg.appendChild(plot);

  var line = document.createElementNS(NS, "line");
  line.setAttribute("x1", clampX(X(ext[0]))); line.setAttribute("x2", clampX(X(ext[1])));
  line.setAttribute("y1", yLine); line.setAttribute("y2", yLine);
  line.setAttribute("class", "intron-line");
  plot.appendChild(line);

  var step = innerW / 30, dir = rec.strand === "+" ? 1 : -1;
  for (var xa = PAD + 6; xa < W - PAD - 6; xa += step) {
    var ar = document.createElementNS(NS, "path");
    ar.setAttribute("d", "M" + xa + "," + (yLine - 4) + " L" + (xa + dir * 6) + "," + yLine + " L" + xa + "," + (yLine + 4));
    ar.setAttribute("class", "strand-arrow");
    plot.appendChild(ar);
  }
  function drawBox(a, b, y, h, cls, fill, elId) {
    if (b < vs || a > ve) return null;          // fully outside view
    var x0 = X(a), x1 = X(b);
    var r = document.createElementNS(NS, "rect");
    r.setAttribute("x", x0); r.setAttribute("y", y);
    r.setAttribute("width", Math.max(cls === "sse-genome" ? 1.5 : 1, x1 - x0));
    r.setAttribute("height", h); r.setAttribute("class", cls);
    if (fill) r.setAttribute("fill", fill);
    if (elId) r.dataset.elem = elId;
    plot.appendChild(r);
    return r;
  }
  exons.forEach(function (e) { drawBox(e[0], e[1], yExon, exonH, "exon-box"); });
  cds.forEach(function (c) { drawBox(c[0], c[1], yCDS, cdsH, "cds-box"); });
  rec.elements.forEach(function (el) {
    el.genome_ranges.forEach(function (gr) {
      var r = drawBox(gr[0], gr[1], yCDS, cdsH, "sse-genome", el.color, el.id);
      if (!r) return;
      r.addEventListener("mouseenter", function () { highlightElement(el.id, true); showElemTip(el, r); });
      r.addEventListener("mouseleave", function () { highlightElement(el.id, false); hideTip(); });
      r.addEventListener("click", function () { selectElement(el.id); });
    });
  });
  // coordinate labels reflect the CURRENT view window
  [[vs, "start"], [ve, "end"]].forEach(function (t) {
    var tx = document.createElementNS(NS, "text");
    tx.setAttribute("x", t[1] === "start" ? PAD : W - PAD); tx.setAttribute("y", yExon - 10);
    tx.setAttribute("class", "coord-label");
    if (t[1] === "end") tx.setAttribute("text-anchor", "end");
    tx.textContent = rec.scaffold + ":" + Math.round(t[0]).toLocaleString();
    svg.appendChild(tx);
  });
  [["exon", yExon + exonH / 2 + 4], ["CDS", yCDS + cdsH / 2 + 4]].forEach(function (l) {
    var tx = document.createElementNS(NS, "text");
    tx.setAttribute("x", 6); tx.setAttribute("y", l[1]);
    tx.setAttribute("class", "row-label"); tx.textContent = l[0];
    svg.appendChild(tx);
  });
  // zoom indicator
  var full = ext[1] - ext[0], z = full / vspan;
  var zt = document.createElementNS(NS, "text");
  zt.setAttribute("x", W - PAD); zt.setAttribute("y", H - 6);
  zt.setAttribute("text-anchor", "end"); zt.setAttribute("class", "row-label");
  zt.textContent = (z > 1.01 ? z.toFixed(1) + "\u00d7 zoom \u2022 " : "") +
    Math.round(vspan).toLocaleString() + " bp \u2022 scroll to zoom, drag to pan, dbl-click reset";
  svg.appendChild(zt);

  // per-render geometry the (once-installed) interaction handlers read
  STATE.gGeom = { svg: svg, X: X, vs: vs, ve: ve, innerW: innerW, PAD: PAD, W: W };
  // wheel + dblclick live on the freshly-created svg (replaced each render)
  svg.addEventListener("wheel", onGenomeWheel, { passive: false });
  svg.addEventListener("mousedown", onGenomeMouseDown);
  svg.addEventListener("dblclick", function (ev) { ev.preventDefault(); resetGenomeView(); });
  host.appendChild(svg);
}

function gAtClientX(clientX) {
  var g = STATE.gGeom;
  var rect = g.svg.getBoundingClientRect();
  var frac = (clientX - rect.left - g.PAD) / g.innerW;
  return g.vs + frac * (g.ve - g.vs);
}
function onGenomeWheel(ev) {
  ev.preventDefault();
  var g = STATE.gGeom, ext = genomeExtent(), vspan = g.ve - g.vs;
  var focus = Math.max(g.vs, Math.min(g.ve, gAtClientX(ev.clientX)));
  var factor = ev.deltaY < 0 ? 0.8 : 1.25;
  var newSpan = Math.min(ext[1] - ext[0], Math.max(6, vspan * factor));
  var leftFrac = (focus - g.vs) / vspan;
  var ns = focus - leftFrac * newSpan;
  STATE.gview = { start: ns, end: ns + newSpan };
  renderGenome();
}
function onGenomeMouseDown(ev) {
  if (ev.target.classList && ev.target.classList.contains("sse-genome")) return; // allow select
  STATE._drag = { on: true, lastX: ev.clientX };
  STATE.gGeom.svg.style.cursor = "grabbing";
}
// window-level drag handlers installed ONCE (see init)
function onWindowMouseMove(ev) {
  if (!STATE._drag || !STATE._drag.on) return;
  var g = STATE.gGeom, vspan = g.ve - g.vs, ext = genomeExtent();
  var dpx = ev.clientX - STATE._drag.lastX; STATE._drag.lastX = ev.clientX;
  var dg = -(dpx / g.innerW) * vspan;
  var ns = g.vs + dg, ne = g.ve + dg;
  if (ns < ext[0]) { ne += (ext[0] - ns); ns = ext[0]; }
  if (ne > ext[1]) { ns -= (ne - ext[1]); ne = ext[1]; }
  STATE.gview = { start: ns, end: ne };
  renderGenome();
}
function onWindowMouseUp() {
  if (STATE._drag) STATE._drag.on = false;
  if (STATE.gGeom && STATE.gGeom.svg) STATE.gGeom.svg.style.cursor = "";
}

/* =====================================================================
   LEGEND
   ===================================================================== */
function renderLegend() {
  var rec = STATE.rec, host = document.getElementById("legend");
  host.innerHTML = "";
  if (!rec.elements.length) { host.textContent = "No structural elements detected in this model."; return; }
  rec.elements.forEach(function (el) {
    var chip = document.createElement("span");
    chip.className = "chip"; chip.dataset.elem = el.id;
    var sw = document.createElement("span"); sw.className = "sw"; sw.style.background = el.color;
    chip.appendChild(sw);
    var lab = document.createElement("span");
    lab.textContent = el.id + " " + el.kind + " " + el.prot_start + "\u2013" + el.prot_end;
    chip.appendChild(lab);
    chip.addEventListener("mouseenter", function () { highlightElement(el.id, true); });
    chip.addEventListener("mouseleave", function () { highlightElement(el.id, false); });
    chip.addEventListener("click", function () { selectElement(el.id); });
    host.appendChild(chip);
  });
}

/* =====================================================================
   CROSS-PANEL LINKING
   ===================================================================== */
function highlightElement(elemId, on) {
  document.querySelectorAll('#seqTrack .aa[data-elem="' + elemId + '"]').forEach(function (s) { s.classList.toggle("hl", on); });
  document.querySelectorAll('#genomeTrack .sse-genome[data-elem="' + elemId + '"]').forEach(function (r) { r.classList.toggle("hl", on); });
  document.querySelectorAll('#legend .chip[data-elem="' + elemId + '"]').forEach(function (c) { c.classList.toggle("hl", on); });
  highlight3D(elemId, on);
}
function selectElement(elemId) {
  STATE.selected = (STATE.selected === elemId) ? null : elemId;
  document.querySelectorAll(".sel").forEach(function (n) { n.classList.remove("sel"); });
  if (STATE.selected) {
    document.querySelectorAll('[data-elem="' + elemId + '"]').forEach(function (n) { n.classList.add("sel"); });
    focus3D(elemId);
  }
}
function elementModelRange(elemId) {
  var rec = STATE.rec;
  var el = rec.elements.filter(function (e) { return e.id === elemId; })[0];
  if (!el) return null;
  var off = (rec.model && rec.model.offset) || 0;
  // full-protein 1-based prot_start/end -> model auth_seq_id (1-based) = protRes - off
  return [el.prot_start - off, el.prot_end - off];
}
function residueRangeLoci(startAuth, endAuth) {
  var plugin = STATE.plugin;
  var s = plugin.managers.structure.hierarchy.current.structures[0];
  if (!s || !s.cell || !s.cell.obj) return null;
  var structure = s.cell.obj.data;
  var MS = window.molstar;
  try {
    var sel = MS.Script.getStructureSelection(function (Q) {
      return Q.struct.generator.atomGroups({
        "residue-test": Q.core.rel.inRange([
          Q.struct.atomProperty.macromolecular.auth_seq_id(), startAuth, endAuth]),
      });
    }, structure);
    return MS.StructureSelection.toLociWithSourceUnits(sel);
  } catch (e) { return null; }
}
function highlight3D(elemId, on) {
  var plugin = STATE.plugin;
  if (!plugin || !STATE.structureRef) return;
  try {
    if (!on) { plugin.managers.interactivity.lociHighlights.clearHighlights(); return; }
    var rng = elementModelRange(elemId); if (!rng) return;
    var loci = residueRangeLoci(rng[0], rng[1]);
    if (loci) plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: loci });
  } catch (e) {}
}
function focus3D(elemId) {
  var plugin = STATE.plugin;
  if (!plugin || !STATE.structureRef) return;
  try {
    var rng = elementModelRange(elemId); if (!rng) return;
    var loci = residueRangeLoci(rng[0], rng[1]);
    if (loci) plugin.managers.camera.focusLoci(loci);
  } catch (e) {}
}

/* =====================================================================
   TOOLTIPS
   ===================================================================== */
function showResidueTip(span) {
  var rec = STATE.rec, i = parseInt(span.dataset.res, 10) - 1;
  var g = rec.residue_genome[i];
  var txt = "Residue " + (i + 1) + " " + rec.protein_sequence[i];
  if (span.dataset.elem) txt += " \u2022 " + span.dataset.elem;
  if (g) txt += " \u2022 " + rec.scaffold + ":" + Math.min(g[0], g[2]).toLocaleString() + "-" + Math.max(g[0], g[2]).toLocaleString();
  tip(txt, span);
}
function showElemTip(el, node) {
  var ranges = el.genome_ranges.map(function (r) { return r[0].toLocaleString() + "-" + r[1].toLocaleString(); }).join(", ");
  tip(el.id + " " + el.kind + " \u2022 protein " + el.prot_start + "-" + el.prot_end +
      " \u2022 " + STATE.rec.scaffold + ":" + ranges +
      " (" + el.genome_ranges.length + " exon segment" + (el.genome_ranges.length > 1 ? "s" : "") + ")", node);
}
function tip(text, node) {
  var t = document.getElementById("tooltip");
  t.textContent = text; t.style.display = "block";
  var r = node.getBoundingClientRect();
  t.style.left = Math.min(window.innerWidth - t.offsetWidth - 10, r.left) + "px";
  t.style.top = (r.bottom + 6 + window.scrollY) + "px";
}
function hideTip() { document.getElementById("tooltip").style.display = "none"; }

/* =====================================================================
   TRANSCRIPT SWITCHING + INIT
   ===================================================================== */
async function selectTranscript(tid) {
  STATE.rec = loadRecord(tid);
  STATE.gview = null;              // reset zoom/pan for new locus
  document.getElementById("meta").innerHTML = metaHTML(STATE.rec);
  renderSequence(); renderGenome(); renderLegend();
  await loadStructure(tid);
}
function metaHTML(rec) {
  return "<b>" + rec.transcript_id + "</b> &middot; gene " + rec.gene_id +
    " &middot; " + rec.scaffold + ":" + rec.span[0].toLocaleString() + "-" + rec.span[1].toLocaleString() +
    " (" + rec.strand + ") &middot; " + rec.protein_length + " aa &middot; " +
    rec.n_helix + " helices, " + rec.n_strand + " strands" +
    (rec.cds_translation_ok ? " &middot; <span class='ok'>CDS&#10003;</span>" : " &middot; <span class='warn'>CDS mismatch</span>");
}
function populateSelector() {
  var sel = document.getElementById("transcriptSelect");
  STATE.index.transcripts.forEach(function (t) {
    var o = document.createElement("option");
    o.value = t.transcript_id;
    o.textContent = t.transcript_id + "  [" + t.scaffold + ", " + t.n_helix + "H/" + t.n_strand + "S, " + t.protein_length + "aa]";
    sel.appendChild(o);
  });
  sel.addEventListener("change", function () { selectTranscript(sel.value); });
}
async function init() {
  STATE.index = JSON.parse(decodeAsset("index-json"));
  populateSelector();
  var viewer = await molstar.Viewer.create("molstar", {
    layoutIsExpanded: false, layoutShowControls: false, layoutShowSequence: false,
    layoutShowLog: false, layoutShowLeftPanel: false,
    viewportShowExpand: true, viewportShowSelectionMode: false, viewportShowAnimation: false,
  });
  STATE.plugin = viewer.plugin ? viewer.plugin : viewer;
  document.getElementById("exportBtn").addEventListener("click", exportGenBankUI);
  window.addEventListener("resize", function () { if (STATE.rec) renderGenome(); });
  // once-installed genome pan handlers
  window.addEventListener("mousemove", onWindowMouseMove);
  window.addEventListener("mouseup", onWindowMouseUp);
  // sequence zoom controls
  var sz = document.getElementById("seqZoomIn"), so = document.getElementById("seqZoomOut"), sr = document.getElementById("seqZoomReset");
  if (sz) sz.addEventListener("click", function () { STATE.seqZoom = Math.min(4, (STATE.seqZoom || 1) * 1.25); renderSequence(); });
  if (so) so.addEventListener("click", function () { STATE.seqZoom = Math.max(0.5, (STATE.seqZoom || 1) / 1.25); renderSequence(); });
  if (sr) sr.addEventListener("click", function () { STATE.seqZoom = 1; renderSequence(); });
  // Ctrl/Cmd + wheel over the sequence zooms it
  document.getElementById("seqTrack").addEventListener("wheel", function (ev) {
    if (!(ev.ctrlKey || ev.metaKey)) return;
    ev.preventDefault();
    STATE.seqZoom = Math.min(4, Math.max(0.5, (STATE.seqZoom || 1) * (ev.deltaY < 0 ? 1.1 : 0.9)));
    renderSequence();
  }, { passive: false });
  var first = STATE.index.transcripts[0].transcript_id;
  document.getElementById("transcriptSelect").value = first;
  await selectTranscript(first);
}

/* =====================================================================
   GENBANK EXPORT (client-side feature table for the current transcript)
   ===================================================================== */
function exportGenBankUI() {
  var rec = STATE.rec; if (!rec) return;
  var gb = buildGenBank(rec);
  var blob = new Blob([gb], { type: "text/plain" });
  var a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = rec.transcript_id.replace(/[^A-Za-z0-9._-]/g, "_") + ".gb";
  a.click();
}
function buildGenBank(rec) {
  var exons = rec.exons.slice().sort(function (a, b) { return a[0] - b[0]; });
  var gmin = exons[0][0], gmax = exons[exons.length - 1][1], L = gmax - gmin + 1;
  var lines = [];
  lines.push("LOCUS       " + (rec.transcript_id.slice(0, 16)).padEnd(16) + " " + L + " bp    DNA     linear   UNK");
  lines.push("DEFINITION  " + rec.transcript_id + " structural-feature annotation (" +
             rec.scaffold + ":" + gmin + "-" + gmax + ").");
  lines.push("ACCESSION   " + rec.scaffold);
  lines.push("KEYWORDS    genostruct; AlphaFold; secondary structure.");
  lines.push("FEATURES             Location/Qualifiers");
  lines.push("     source          1.." + L);
  lines.push("                     /note=\"" + rec.scaffold + ":" + gmin + "-" + gmax + "\"");
  lines.push("     gene            " + (gmin - gmin + 1) + ".." + (gmax - gmin + 1));
  lines.push("                     /gene=\"" + rec.gene_id + "\"");
  var cds = rec.cds.slice().sort(function (a, b) { return a[0] - b[0]; });
  var parts = cds.map(function (c) { return (c[0] - gmin + 1) + ".." + (c[1] - gmin + 1); });
  var jn = parts.length > 1 ? "join(" + parts.join(",") + ")" : parts[0];
  if (rec.strand === "-") jn = "complement(" + jn + ")";
  lines.push("     CDS             " + jn);
  lines.push("                     /transcript_id=\"" + rec.transcript_id + "\"");
  rec.elements.forEach(function (el) {
    var ep = el.genome_ranges.map(function (gr) { return (gr[0] - gmin + 1) + ".." + (gr[1] - gmin + 1); });
    var ej = ep.length > 1 ? "join(" + ep.join(",") + ")" : ep[0];
    if (rec.strand === "-") ej = "complement(" + ej + ")";
    lines.push("     misc_feature    " + ej);
    lines.push("                     /note=\"" + el.kind + " " + el.id + " (protein " + el.prot_start + "-" + el.prot_end + ")\"");
    lines.push("                     /label=\"" + el.id + "\"");
    lines.push("                     /color=\"" + el.color + "\"");
  });
  lines.push("ORIGIN");
  lines.push("//");
  lines.push("");
  lines.push("; NOTE: genomic sequence omitted from the in-browser export to keep");
  lines.push("; the viewer lightweight. For a complete GenBank record WITH sequence,");
  lines.push("; use genostruct_export.py (see README).");
  return lines.join("\n");
}

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
else init();

})();
