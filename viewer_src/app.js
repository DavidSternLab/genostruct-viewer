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
  plddtThemeName: "genostruct-plddt",
  selected: null,
  dimElem: null,     // element to focus in 3D (rest dimmed)
  hoverElem: null,   // element hovered (brief 3D emphasis)
  focusRange: null,  // {start,end} 1-based protein residue drag-selection
  plddtMode: false,  // false = element colors, true = AF2 pLDDT confidence
  flipGenome: true,  // true = show genome in protein N->C orientation (default)
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
// blend an integer color toward light grey by factor f in [0,1] (1 = full grey)
function dimInt(intColor, f) {
  var r = (intColor >> 16) & 255, g = (intColor >> 8) & 255, b = intColor & 255;
  var G = 224; // light grey target
  r = Math.round(r + (G - r) * f); g = Math.round(g + (G - g) * f); b = Math.round(b + (G - b) * f);
  return (r << 16) | (g << 8) | b;
}

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
        var rec = STATE.rec;
        var recColors = rec ? rec.residue_color : null;
        if (!recColors) return DEFAULT;
        var off = (rec.model && rec.model.offset) || 0;
        var unit = location.unit, el = location.element;
        if (!unit || el === undefined || el === null) return DEFAULT;
        var residueIndex;
        try { residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index[el]; }
        catch (e) { residueIndex = undefined; }
        if (residueIndex === undefined || residueIndex === null) return DEFAULT;
        // residueIndex is 0-based within model; map to full protein residue.
        var protRes = off + residueIndex;
        var base = DEFAULT;
        if (protRes >= 0 && protRes < recColors.length) {
          var c = recColors[protRes];
          if (c) base = hexToInt(c);
        }
        // focus/dim: when an element is focused, keep its residues vivid and
        // fade everything else toward light grey. A hovered element (no focus)
        // gets a lighter emphasis without dimming the rest.
        // (a) residue-range focus (drag-select on sequence): 1-based inclusive
        var fr = STATE.focusRange;
        if (fr) {
          var res1 = protRes + 1;
          if (res1 >= fr.start && res1 <= fr.end) return base;
          return dimInt(base, base === DEFAULT ? 0.55 : 0.82);
        }
        // (b) element focus (click)
        var focus = STATE.dimElem;
        if (focus) {
          var reElem = rec.residue_element;
          var thisElem = (reElem && protRes >= 0 && protRes < reElem.length) ? reElem[protRes] : null;
          if (thisElem === focus) return base;              // focused element: full color
          return dimInt(base, base === DEFAULT ? 0.55 : 0.82); // everything else: dimmed
        }
        return base;
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
  // Second theme: AF2 pLDDT confidence, colored from the per-residue B-factors we
  // stored in the pipeline (model.plddt). Self-contained so it works on plain PDB
  // input, where Mol*'s built-in plddt-confidence theme is not applicable.
  var plddtProvider = {
    name: STATE.plddtThemeName,
    label: "pLDDT confidence",
    category: "Custom",
    factory: function (ctx, props) {
      function colorFor(location) {
        var rec = STATE.rec;
        var plddt = rec && rec.model ? rec.model.plddt : null;
        if (!plddt) return DEFAULT;
        var unit = location.unit, el = location.element;
        if (!unit || el == null) return DEFAULT;
        var ri;
        try { ri = unit.model.atomicHierarchy.residueAtomSegments.index[el]; }
        catch (e) { ri = undefined; }
        if (ri == null || ri < 0 || ri >= plddt.length) return DEFAULT;
        return plddtColor(plddt[ri]);
      }
      return { factory: plddtProvider.factory, granularity: "group", color: colorFor,
               props: props || {}, description: "AF2 pLDDT confidence" };
    },
    getParams: function () { return {}; },
    defaultValues: {},
    isApplicable: function () { return true; },
  };
  try {
    reg.add(provider);
    reg.add(plddtProvider);
    STATE._themeRegistered = true;
  } catch (e) {
    // Already registered from a prior call (e.g. hot reload) — treat as success
    // so we never retry and never let the throw propagate into loadStructure.
    STATE._themeRegistered = true;
  }
}
// AlphaFold pLDDT color ramp -> branded int
function plddtColor(v) {
  if (v == null) return 0xBFBFBF;
  if (v >= 90) return 0x0053D6;   // very high (dark blue)
  if (v >= 70) return 0x65CBF3;   // confident (cyan)
  if (v >= 50) return 0xFFDB13;   // low (yellow)
  return 0xFF7D45;                // very low (orange)
}

function activeThemeName() { return STATE.plddtMode ? STATE.plddtThemeName : STATE.colorThemeName; }
async function loadStructure(tid) {
  var plugin = STATE.plugin;
  var status = document.getElementById("structStatus");
  // token guards against a stale async load (rapid transcript switching): if a
  // newer load started, this one stops touching the UI. This also prevents the
  // spurious "Could not find node" error that fired from a late op on a cleared
  // node even though the current structure rendered fine.
  var myToken = (STATE._loadToken = (STATE._loadToken || 0) + 1);
  var rendered = false;
  if (status) status.textContent = "";
  try {
    registerElementTheme(plugin);
    await plugin.clear();
    var pdbText = loadPDB(tid);
    if (!pdbText) throw new Error("no embedded PDB for " + tid);
    var data = await plugin.builders.data.rawData({ data: pdbText, label: tid });
    var traj = await plugin.builders.structure.parseTrajectory(data, "pdb");
    var model = await plugin.builders.structure.createModel(traj);
    var structure = await plugin.builders.structure.createStructure(model);
    try {
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: "cartoon", color: activeThemeName(),
      });
      rendered = true;
    } catch (themeErr) {
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: "cartoon", color: "chain-id",
      });
      rendered = true;
      if (myToken === STATE._loadToken && status) status.textContent = "(default coloring — custom theme unavailable)";
    }
    STATE.structureRef = structure;
    try { plugin.managers.camera.reset(); } catch (camErr) { /* late/stale camera op — harmless */ }
  } catch (err) {
    // Only surface an error if the structure did NOT render AND this is still the
    // current load. Stale-node errors from a superseded load are ignored.
    if (!rendered && myToken === STATE._loadToken) {
      if (status) status.textContent = "Could not display structure for " + tid + ": " + (err && err.message ? err.message : err);
      STATE.structureRef = null;
    }
    if (window.console) console.warn("loadStructure note", err && err.message);
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
        if (STATE._seqDrag && STATE._seqDrag.on) { extendSeqDrag(+span.dataset.res); return; }
        if (span.dataset.elem) highlightElement(span.dataset.elem, true);
        showResidueTip(span);
      });
      span.addEventListener("mouseleave", function () {
        if (STATE._seqDrag && STATE._seqDrag.on) return;
        if (span.dataset.elem) highlightElement(span.dataset.elem, false);
        hideTip();
      });
      span.addEventListener("mousedown", function (ev) {
        ev.preventDefault();
        STATE._seqDrag = { on: true, anchor: +span.dataset.res, moved: false };
      });
    })(span);
    wrap.appendChild(span);
  }
  host.appendChild(wrap);
  // apply any existing drag-selection highlight
  paintSeqSelection();
}
// --- protein sequence drag-select: highlights a residue range, zooms the genome
//     to the encoding bases, and focuses/dims the 3D structure to that range ---
function extendSeqDrag(res) {
  var d = STATE._seqDrag; if (!d || !d.on) return;
  d.moved = true; d.cur = res;
  paintSeqSelection();
}
function paintSeqSelection() {
  var d = STATE._seqDrag, fr = STATE.focusRange;
  var lo, hi;
  if (d && d.on && d.cur != null) { lo = Math.min(d.anchor, d.cur); hi = Math.max(d.anchor, d.cur); }
  else if (fr) { lo = fr.start; hi = fr.end; }
  document.querySelectorAll('#seqTrack .aa').forEach(function (s) {
    var r = +s.dataset.res;
    s.classList.toggle("rangesel", lo != null && r >= lo && r <= hi);
  });
}
function finishSeqDrag() {
  var d = STATE._seqDrag; if (!d || !d.on) return;
  d.on = false;
  if (!d.moved || d.cur == null) {
    // treat as a click: select the element under the residue, if any
    var rec = STATE.rec, idx = d.anchor - 1;
    var el = rec.residue_element && rec.residue_element[idx];
    STATE._seqDrag = null;
    if (el) selectElement(el);
    return;
  }
  var lo = Math.min(d.anchor, d.cur), hi = Math.max(d.anchor, d.cur);
  STATE._seqDrag = null;
  applyResidueRangeFocus(lo, hi);
}
// Focus a protein residue range [lo,hi] (1-based): dim 3D outside it, zoom the
// genome to the genomic bases encoding those residues, and mark the sequence.
function applyResidueRangeFocus(lo, hi) {
  STATE.selected = null; STATE.dimElem = null;      // range focus supersedes element focus
  STATE.focusRange = { start: lo, end: hi };
  paintSeqSelection();
  recolor3D();
  zoomGenomeToResidues(lo, hi);
  var el = document.getElementById("focusInfo");
  if (el) el.textContent = "Focused protein residues " + lo + "\u2013" + hi + " (click empty area or Reset to clear)";
}
function clearRangeFocus() {
  STATE.focusRange = null;
  paintSeqSelection();
  recolor3D();
  resetGenomeView();
  var el = document.getElementById("focusInfo"); if (el) el.textContent = "";
}
// map a protein residue range to genomic coordinates via residue_genome, then
// zoom the genome track to the enclosing window (handles minus strand + splits).
function zoomGenomeToResidues(lo, hi) {
  var rec = STATE.rec, rg = rec.residue_genome;
  if (!rg) return;
  var mn = Infinity, mx = -Infinity;
  for (var r = lo; r <= hi; r++) {
    var codon = rg[r - 1];
    if (!codon) continue;
    codon.forEach(function (p) { if (p < mn) mn = p; if (p > mx) mx = p; });
  }
  if (mn === Infinity) return;
  var ext = genomeExtent();
  var span = mx - mn, pad = Math.max(30, span * 0.4);
  STATE.gview = { start: Math.max(ext[0], mn - pad), end: Math.min(ext[1], mx + pad) };
  renderGenome();
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
// Zoom the genome track to a window enclosing the selected element's genomic
// ranges, with padding, clamped to the locus extent.
function zoomGenomeToElement(elemId) {
  var rec = STATE.rec;
  var el = rec.elements.filter(function (e) { return e.id === elemId; })[0];
  if (!el || !el.genome_ranges || !el.genome_ranges.length) return;
  var lo = Infinity, hi = -Infinity;
  el.genome_ranges.forEach(function (gr) { if (gr[0] < lo) lo = gr[0]; if (gr[1] > hi) hi = gr[1]; });
  var ext = genomeExtent();
  var span = hi - lo, pad = Math.max(30, span * 0.5);   // at least 30 bp of context
  var s = Math.max(ext[0], lo - pad), e = Math.min(ext[1], hi + pad);
  if (e - s < 12) { var mid = (s + e) / 2; s = Math.max(ext[0], mid - 6); e = Math.min(ext[1], mid + 6); }
  STATE.gview = { start: s, end: e };
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
  // orientation: when flipGenome is on and the gene is on the minus strand, draw
  // the axis high->low so the genome reads in the same N->C direction as the
  // protein. `flip` also drives the interaction handlers via STATE.gGeom.
  var flip = !!(STATE.flipGenome && rec.strand === "-");
  function X(g) {
    var f = (g - vs) / vspan;
    if (flip) f = 1 - f;
    return PAD + f * innerW;
  }
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

  // strand arrows point in the direction of transcription AS DRAWN. With flip on,
  // a minus-strand gene reads left->right (N->C), so arrows point right (+1).
  var txnDir = rec.strand === "+" ? 1 : -1;
  var dir = flip ? -txnDir : txnDir;
  var step = innerW / 30;
  for (var xa = PAD + 6; xa < W - PAD - 6; xa += step) {
    var ar = document.createElementNS(NS, "path");
    ar.setAttribute("d", "M" + xa + "," + (yLine - 4) + " L" + (xa + dir * 6) + "," + yLine + " L" + xa + "," + (yLine + 4));
    ar.setAttribute("class", "strand-arrow");
    plot.appendChild(ar);
  }
  function drawBox(a, b, y, h, cls, fill, elId) {
    if (b < vs || a > ve) return null;          // fully outside view
    var xa2 = X(a), xb2 = X(b), x0 = Math.min(xa2, xb2), x1 = Math.max(xa2, xb2);
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
  // coordinate labels reflect the CURRENT view window. With flip on, the LEFT
  // edge of the drawing is the higher coordinate (ve) and the right edge is vs.
  var leftCoord = flip ? ve : vs, rightCoord = flip ? vs : ve;
  [[leftCoord, "start"], [rightCoord, "end"]].forEach(function (t) {
    var tx = document.createElementNS(NS, "text");
    tx.setAttribute("x", t[1] === "start" ? PAD : W - PAD); tx.setAttribute("y", yExon - 10);
    tx.setAttribute("class", "coord-label");
    if (t[1] === "end") tx.setAttribute("text-anchor", "end");
    tx.textContent = rec.scaffold + ":" + Math.round(t[0]).toLocaleString();
    svg.appendChild(tx);
  });
  // orientation caption (left -> right reading direction)
  var ori = document.createElementNS(NS, "text");
  ori.setAttribute("x", W / 2); ori.setAttribute("y", yExon - 10);
  ori.setAttribute("text-anchor", "middle"); ori.setAttribute("class", "row-label");
  ori.textContent = flip
    ? "5'\u21923' protein N\u2192C  (minus-strand gene shown flipped)"
    : (rec.strand === "-" ? "genome 5'\u21923' (minus-strand: protein reads right\u2192left)" : "5'\u21923' protein N\u2192C");
  svg.appendChild(ori);
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
  STATE.gGeom = { svg: svg, X: X, vs: vs, ve: ve, innerW: innerW, PAD: PAD, W: W, flip: flip };
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
  if (g.flip) frac = 1 - frac;   // inverse of the flipped X() mapping
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
  if (g.flip) dg = -dg;   // flipped axis: dragging right moves toward lower coords
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
  // 3D hover emphasis: only when nothing is focus-locked, briefly dim-focus the
  // hovered element; restore to the locked focus (or none) on mouse-out.
  if (!STATE.selected) {
    STATE.dimElem = on ? elemId : null;
    recolor3D();
  }
}
function selectElement(elemId) {
  STATE.selected = (STATE.selected === elemId) ? null : elemId;
  document.querySelectorAll(".sel").forEach(function (n) { n.classList.remove("sel"); });
  if (STATE.selected) {
    document.querySelectorAll('[data-elem="' + elemId + '"]').forEach(function (n) { n.classList.add("sel"); });
    STATE.dimElem = elemId;      // focus in 3D: this element vivid, rest dimmed
    recolor3D();
    zoomGenomeToElement(elemId); // auto-zoom the genome track around this element
  } else {
    STATE.dimElem = null;
    recolor3D();
    resetGenomeView();
  }
  updateFocusBanner();
}
function updateFocusBanner() {
  var el = document.getElementById("focusInfo");
  if (!el) return;
  if (STATE.selected) {
    var e = STATE.rec.elements.filter(function (x) { return x.id === STATE.selected; })[0];
    el.textContent = e ? ("Focused " + e.id + " (" + e.kind + ", protein " + e.prot_start + "\u2013" + e.prot_end + ") \u2014 click again to clear") : "";
  } else { el.textContent = ""; }
}
// Re-run the color theme so the structure re-colors/dims WITHOUT re-adding the
// theme (the loci/selection APIs are not exported by the Viewer UMD; recoloring
// via the registered theme is the reliable path and also drives the dim effect).
function recolor3D() {
  // focus/dim only makes sense in element-color mode; leave pLDDT view untouched
  if (STATE.plddtMode) return;
  applyTheme(STATE.colorThemeName);
}
// Switch the structure between element colors and the custom pLDDT theme. The
// theme change reliably re-invokes the factory; if the in-place update path ever
// fails we rebuild the representation (the known-good initial-render path).
function setStructureColorMode() {
  if (!applyTheme(activeThemeName()) && STATE.rec) {
    loadStructure(STATE.rec.transcript_id);
  }
}
function applyTheme(themeName) {
  var plugin = STATE.plugin;
  if (!plugin) return false;
  try {
    var comps = plugin.managers.structure.hierarchy.current.structures[0];
    comps = comps ? comps.components : null;
    if (!comps || !comps.length) return false;
    // A changing colorParams (nonce) is REQUIRED: Mol*'s SW() builds the theme
    // descriptor as {name, params}, and the state reconciler skips the update
    // when {name, params} is unchanged — so recoloring with the same theme name
    // and empty params is a no-op and the factory never re-reads STATE (focus/
    // dim/range). Bumping a nonce makes params differ each call, forcing the
    // recompute. (Our theme getParams() ignores unknown keys, so nonce is inert
    // aside from busting the equality check.)
    STATE._recolorNonce = (STATE._recolorNonce || 0) + 1;
    plugin.managers.structure.component.updateRepresentationsTheme(comps, {
      color: themeName, colorParams: { nonce: STATE._recolorNonce },
    });
    return true;
  } catch (e) { if (window.console) console.warn("applyTheme failed", e); return false; }
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
  STATE.selected = null; STATE.dimElem = null; STATE.hoverElem = null;  // clear focus
  STATE.focusRange = null; STATE._seqDrag = null;
  updateFocusBanner();
  var selEl = document.getElementById("transcriptSelect");
  if (selEl && selEl.value !== tid) selEl.value = tid;
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
function transcriptSortKey(tid) {
  // numeric-aware key: split into runs of digits / non-digits so g113516.t1-00001
  // orders by gene number, then isoform, naturally (g9 before g100).
  var parts = String(tid).match(/\d+|\D+/g) || [tid];
  return parts.map(function (p) { return /^\d+$/.test(p) ? p.padStart(12, "0") : p; }).join("");
}
function optionLabel(t) {
  return t.transcript_id + "  [" + t.scaffold + ", " + t.n_helix + "H/" + t.n_strand + "S, " + t.protein_length + "aa]";
}
function populateSelector() {
  var sel = document.getElementById("transcriptSelect");
  var search = document.getElementById("transcriptSearch");
  var dlist = document.getElementById("transcriptOptions");
  // sort transcripts numerically and keep this order everywhere
  STATE.index.transcripts.sort(function (a, b) {
    var ka = transcriptSortKey(a.transcript_id), kb = transcriptSortKey(b.transcript_id);
    return ka < kb ? -1 : ka > kb ? 1 : 0;
  });
  sel.innerHTML = "";
  STATE.index.transcripts.forEach(function (t) {
    var o = document.createElement("option");
    o.value = t.transcript_id;
    o.textContent = optionLabel(t);
    sel.appendChild(o);
  });
  sel.addEventListener("change", function () {
    selectTranscript(sel.value);
    if (search) search.value = sel.value;
  });

  // ---- type-ahead search box (native datalist) ----
  if (!search || !dlist) return;
  function rebuildDatalist(q) {
    q = (q || "").toLowerCase();
    dlist.innerHTML = "";
    var shown = 0;
    for (var i = 0; i < STATE.index.transcripts.length && shown < 60; i++) {
      var t = STATE.index.transcripts[i];
      if (!q || t.transcript_id.toLowerCase().indexOf(q) !== -1) {
        var o = document.createElement("option");
        o.value = t.transcript_id;
        o.label = optionLabel(t);
        dlist.appendChild(o); shown++;
      }
    }
  }
  rebuildDatalist("");
  var idset = {};
  STATE.index.transcripts.forEach(function (t) { idset[t.transcript_id] = true; });
  function tryGo(v) {
    v = (v || "").trim();
    if (idset[v]) { sel.value = v; selectTranscript(v); return true; }
    // if a unique prefix/substring match exists, jump to it
    var matches = STATE.index.transcripts.filter(function (t) {
      return t.transcript_id.toLowerCase().indexOf(v.toLowerCase()) !== -1;
    });
    if (v && matches.length === 1) {
      sel.value = matches[0].transcript_id; search.value = matches[0].transcript_id;
      selectTranscript(matches[0].transcript_id); return true;
    }
    return false;
  }
  search.addEventListener("input", function () {
    rebuildDatalist(search.value);
    tryGo(search.value);       // picking a datalist entry fires 'input' with the full id
  });
  search.addEventListener("change", function () { tryGo(search.value); });
  search.addEventListener("keydown", function (ev) {
    if (ev.key === "Enter") { ev.preventDefault(); tryGo(search.value); }
  });
  // On focus, select all text and expose the FULL list so the user can browse
  // every protein without first deleting the current name.
  search.addEventListener("focus", function () {
    rebuildDatalist("");       // full list available immediately
    search.select();           // typing replaces the current name
  });
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
  var evb = document.getElementById("exportViewBtn");
  if (evb) evb.addEventListener("click", exportViewedRegionUI);
  // color-mode toggle: element colors <-> AF2 pLDDT confidence
  var ctog = document.getElementById("colorModeBtn");
  function labelColorBtn() {
    // buttons label the ACTION (what a press does), not the current state
    if (!ctog) return;
    ctog.textContent = STATE.plddtMode ? "Show element colors" : "Show pLDDT";
    ctog.title = STATE.plddtMode
      ? "Currently AF2 pLDDT confidence \u2014 switch to structural-element colors"
      : "Currently structural-element colors \u2014 switch to AF2 pLDDT confidence";
  }
  if (ctog) ctog.addEventListener("click", function () {
    STATE.plddtMode = !STATE.plddtMode;
    labelColorBtn();
    setStructureColorMode();
  });
  labelColorBtn();
  // genome orientation toggle
  var otog = document.getElementById("orientBtn");
  function labelOrientBtn() {
    if (!otog) return;
    otog.textContent = STATE.flipGenome ? "Show genome 5'\u21923'" : "Show protein N\u2192C";
    otog.title = STATE.flipGenome
      ? "Currently protein N\u2192C orientation \u2014 switch to native genome 5'\u21923'"
      : "Currently native genome 5'\u21923' \u2014 switch to protein N\u2192C orientation";
  }
  if (otog) otog.addEventListener("click", function () {
    STATE.flipGenome = !STATE.flipGenome;
    labelOrientBtn();
    renderGenome();
  });
  labelOrientBtn();
  window.addEventListener("resize", function () { if (STATE.rec) renderGenome(); });
  // once-installed genome pan handlers
  window.addEventListener("mousemove", onWindowMouseMove);
  window.addEventListener("mouseup", onWindowMouseUp);
  // once-installed sequence drag-select finish
  window.addEventListener("mouseup", finishSeqDrag);
  // double-click the sequence clears a range focus
  document.getElementById("seqTrack").addEventListener("dblclick", function () {
    if (STATE.focusRange) clearRangeFocus();
  });
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
  var searchBox = document.getElementById("transcriptSearch");
  if (searchBox) searchBox.value = first;
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
// Export the CURRENTLY-VIEWED genomic window as a GenBank feature table:
// includes only exon/CDS/element segments that overlap the visible view range.
function exportViewedRegionUI() {
  var rec = STATE.rec; if (!rec) return;
  var ext = genomeExtent();
  var view = STATE.gview || { start: ext[0], end: ext[1] };
  var vs = Math.round(Math.max(ext[0], view.start)), ve = Math.round(Math.min(ext[1], view.end));
  var gb = buildGenBankRegion(rec, vs, ve);
  var blob = new Blob([gb], { type: "text/plain" });
  var a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = rec.scaffold + "_" + vs + "-" + ve + ".gb";
  a.click();
}
function overlaps(a, b, s, e) { return a <= e && b >= s; }
function clipParts(ranges, s, e, origin) {
  // keep only the portions of each [a,b] range inside [s,e]; coords relative to origin (1-based)
  var out = [];
  ranges.forEach(function (r) {
    var a = Math.max(r[0], s), b = Math.min(r[1], e);
    if (a <= b) out.push((a - origin + 1) + ".." + (b - origin + 1));
  });
  return out;
}
function buildGenBankRegion(rec, vs, ve) {
  var L = ve - vs + 1;
  var lines = [];
  lines.push("LOCUS       " + (rec.scaffold.slice(0, 16)).padEnd(16) + " " + L + " bp    DNA     linear   UNK");
  lines.push("DEFINITION  " + rec.scaffold + ":" + vs + "-" + ve + " viewed region \u2014 " +
             rec.transcript_id + " structural-feature annotation.");
  lines.push("ACCESSION   " + rec.scaffold);
  lines.push("KEYWORDS    genostruct; AlphaFold; secondary structure; viewed region.");
  lines.push("FEATURES             Location/Qualifiers");
  lines.push("     source          1.." + L);
  lines.push("                     /note=\"" + rec.scaffold + ":" + vs + "-" + ve + " (viewed window)\"");
  // gene box clipped to view
  var exons = rec.exons.slice().sort(function (a, b) { return a[0] - b[0]; });
  var gmin = exons[0][0], gmax = exons[exons.length - 1][1];
  if (overlaps(gmin, gmax, vs, ve)) {
    var gs = Math.max(gmin, vs) - vs + 1, ge = Math.min(gmax, ve) - vs + 1;
    lines.push("     gene            " + gs + ".." + ge);
    lines.push("                     /gene=\"" + rec.gene_id + "\"");
    if (gmin < vs || gmax > ve) lines.push("                     /note=\"gene extends beyond viewed region\"");
  }
  // CDS parts within view
  var cds = rec.cds.slice().sort(function (a, b) { return a[0] - b[0]; });
  var cparts = clipParts(cds, vs, ve, vs);
  if (cparts.length) {
    var cj = cparts.length > 1 ? "join(" + cparts.join(",") + ")" : cparts[0];
    if (rec.strand === "-") cj = "complement(" + cj + ")";
    lines.push("     CDS             " + cj);
    lines.push("                     /transcript_id=\"" + rec.transcript_id + "\"");
    lines.push("                     /note=\"CDS segments within viewed region\"");
  }
  // structural elements overlapping the view
  var nEl = 0;
  rec.elements.forEach(function (el) {
    var parts = clipParts(el.genome_ranges, vs, ve, vs);
    if (!parts.length) return;
    nEl++;
    var ej = parts.length > 1 ? "join(" + parts.join(",") + ")" : parts[0];
    if (rec.strand === "-") ej = "complement(" + ej + ")";
    lines.push("     misc_feature    " + ej);
    lines.push("                     /note=\"" + el.kind + " " + el.id + " (protein " + el.prot_start + "-" + el.prot_end + ")\"");
    lines.push("                     /label=\"" + rec.transcript_id + ":" + el.id + "\"");
    lines.push("                     /color=\"" + el.color + "\"");
  });
  lines.push("ORIGIN");
  lines.push("//");
  lines.push("");
  lines.push("; Viewed region " + rec.scaffold + ":" + vs + "-" + ve + " (" + L + " bp), " + nEl + " structural feature(s).");
  lines.push("; Genomic sequence omitted to keep the viewer lightweight. For a complete");
  lines.push("; GenBank record WITH sequence, use genostruct_export.py export_region().");
  return lines.join("\n");
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
