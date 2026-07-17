"""
build_viewer.py — assemble the self-contained integrated viewer HTML.

Reads the pipeline output (out_dir/data/*.json + out_dir/structures/*.pdb) and
emits a single HTML file with Mol*, pako, all per-transcript data and PDB
coordinates embedded (gzip + base64). No network access required to open it.

Usage:
    python build_viewer.py --out viewer \
        --molstar molstar.js --molstar-css molstar.css --pako pako_inflate.min.js \
        --template viewer_src/template.html --app-js viewer_src/app.js \
        --app-css viewer_src/app.css --html viewer/genostruct_viewer.html
"""
import os, glob, json, gzip, base64, argparse, io

def gz_b64(text: str) -> str:
    if isinstance(text, str): text = text.encode("utf-8")
    return base64.b64encode(gzip.compress(text, 9)).decode("ascii")

# --- PDB -> mmCIF with model secondary structure -------------------------------
# Mol*'s "auto" secondary-structure only reads model SS from the mmCIF
# categories struct_conf (helices) / struct_sheet_range (strands); for a PDB it
# runs DSSP, whose WASM is fetched from a data: URL that the viewer's iframe CSP
# blocks. So we convert each AF2 PDB to mmCIF and add struct_conf /
# struct_sheet_range from the pipeline's own SSE elements. Then hasSecondaryStructure()
# returns true and Mol* uses the model SS -- no DSSP, no blocked fetch.
import numpy as _np
import biotite.structure.io.pdb as _pdbio
import biotite.structure.io.pdbx as _pdbx

def _residue_meta(arr):
    """auth res_id -> {comp, asym}; plus label_seq_id map (sequential per residue)."""
    meta = {}
    for i in range(arr.array_length()):
        rid = int(arr.res_id[i])
        if rid not in meta:
            meta[rid] = {"comp": str(arr.res_name[i]), "asym": (str(arr.chain_id[i]) or "A")}
    uniq = sorted(meta)
    label = {rid: str(k + 1) for k, rid in enumerate(uniq)}
    return meta, label

def _ss_spans(rec, meta):
    """Return (helices, strands) as lists of (beg_auth, end_auth) snapped to
    residue ids that actually exist in the model (AF2 models can have gaps)."""
    off = (rec.get("model") or {}).get("offset", 0) or 0
    if not meta:
        return [], []
    present = sorted(meta)
    lo, hi = present[0], present[-1]
    pset = set(present)

    def snap_up(x):    # smallest present id >= x
        while x <= hi and x not in pset:
            x += 1
        return x if x <= hi else None

    def snap_down(x):  # largest present id <= x
        while x >= lo and x not in pset:
            x -= 1
        return x if x >= lo else None

    helices, strands = [], []
    for el in rec.get("elements", []):
        a = el["prot_start"] - off
        b = el["prot_end"] - off
        if b < a:
            continue
        a = snap_up(max(lo, a))
        b = snap_down(min(hi, b))
        if a is None or b is None or b < a:
            continue
        (helices if el["kind"] == "helix" else strands).append((a, b))
    return helices, strands

def _add_ss_categories(block, rec, arr):
    meta, label = _residue_meta(arr)
    helices, strands = _ss_spans(rec, meta)
    if helices:
        cols = {k: [] for k in (
            "conf_type_id", "id", "pdbx_PDB_helix_class",
            "beg_label_comp_id", "beg_label_asym_id", "beg_label_seq_id",
            "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
            "end_label_comp_id", "end_label_asym_id", "end_label_seq_id",
            "end_auth_comp_id", "end_auth_asym_id", "end_auth_seq_id",
            "pdbx_PDB_helix_length")}
        for k, (a, b) in enumerate(helices, 1):
            cols["conf_type_id"].append("HELX_P"); cols["id"].append("HELX%d" % k); cols["pdbx_PDB_helix_class"].append("1")
            cols["beg_label_comp_id"].append(meta[a]["comp"]); cols["beg_label_asym_id"].append(meta[a]["asym"]); cols["beg_label_seq_id"].append(label[a])
            cols["beg_auth_comp_id"].append(meta[a]["comp"]); cols["beg_auth_asym_id"].append(meta[a]["asym"]); cols["beg_auth_seq_id"].append(str(a))
            cols["end_label_comp_id"].append(meta[b]["comp"]); cols["end_label_asym_id"].append(meta[b]["asym"]); cols["end_label_seq_id"].append(label[b])
            cols["end_auth_comp_id"].append(meta[b]["comp"]); cols["end_auth_asym_id"].append(meta[b]["asym"]); cols["end_auth_seq_id"].append(str(b))
            cols["pdbx_PDB_helix_length"].append(str(b - a + 1))
        cat = _pdbx.CIFCategory()
        for k, v in cols.items():
            cat[k] = _pdbx.CIFColumn(list(map(str, v)))
        block["struct_conf"] = cat
    if strands:
        cols = {k: [] for k in (
            "sheet_id", "id",
            "beg_label_comp_id", "beg_label_asym_id", "beg_label_seq_id",
            "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
            "end_label_comp_id", "end_label_asym_id", "end_label_seq_id",
            "end_auth_comp_id", "end_auth_asym_id", "end_auth_seq_id")}
        for k, (a, b) in enumerate(strands, 1):
            cols["sheet_id"].append("A"); cols["id"].append(str(k))
            cols["beg_label_comp_id"].append(meta[a]["comp"]); cols["beg_label_asym_id"].append(meta[a]["asym"]); cols["beg_label_seq_id"].append(label[a])
            cols["beg_auth_comp_id"].append(meta[a]["comp"]); cols["beg_auth_asym_id"].append(meta[a]["asym"]); cols["beg_auth_seq_id"].append(str(a))
            cols["end_label_comp_id"].append(meta[b]["comp"]); cols["end_label_asym_id"].append(meta[b]["asym"]); cols["end_label_seq_id"].append(label[b])
            cols["end_auth_comp_id"].append(meta[b]["comp"]); cols["end_auth_asym_id"].append(meta[b]["asym"]); cols["end_auth_seq_id"].append(str(b))
        cat = _pdbx.CIFCategory()
        for k, v in cols.items():
            cat[k] = _pdbx.CIFColumn(list(map(str, v)))
        block["struct_sheet_range"] = cat

def pdb_to_cif_with_ss(pdb_path, rec, block_name):
    """Convert an AF2 PDB to mmCIF text carrying model secondary structure."""
    arr = _pdbio.PDBFile.read(pdb_path).get_structure(model=1, extra_fields=["b_factor"])
    cif = _pdbx.CIFFile()
    _pdbx.set_structure(cif, arr, data_block=block_name)   # writes B_iso_or_equiv (pLDDT)
    _add_ss_categories(cif.block, rec, arr)
    buf = io.StringIO()
    cif.write(buf)
    return buf.getvalue()

def build_html(out_dir, molstar_js, molstar_css, pako_js,
               template, app_js, app_css, html_out, limit=None):
    data_dir = os.path.join(out_dir, "data")
    struct_dir = os.path.join(out_dir, "structures")
    index = json.load(open(os.path.join(data_dir, "index.json")))
    txs = index["transcripts"]
    if limit: txs = txs[:limit]

    embed = {"data": {}, "cif": {}}
    for t in txs:
        tid = t["transcript_id"]
        rec_path = os.path.join(data_dir, t["file"])
        rec_bytes = open(rec_path, "rb").read()
        embed["data"][tid] = gz_b64(rec_bytes)
        rec = json.loads(rec_bytes)
        # find pdb (safe-name)
        safe = "".join(c if (c.isalnum() or c in "._-") else "_" for c in tid)
        pdb_path = os.path.join(struct_dir, safe + ".pdb")
        if not os.path.exists(pdb_path):
            # fallback: any pdb whose stem contains tid
            cands = glob.glob(os.path.join(struct_dir, "*.pdb"))
            cands = [c for c in cands if safe in os.path.basename(c)]
            pdb_path = cands[0] if cands else None
        if pdb_path and os.path.exists(pdb_path):
            # convert to mmCIF WITH model secondary structure so Mol* never runs
            # the CSP-blocked DSSP WASM (see pdb_to_cif_with_ss).
            cif_text = pdb_to_cif_with_ss(pdb_path, rec, safe)
            embed["cif"][tid] = gz_b64(cif_text)

    # trim index to the transcripts we actually embedded
    index_slim = {"transcripts": [t for t in txs if t["transcript_id"] in embed["cif"]],
                  "warnings": index.get("warnings", [])[:50]}

    tmpl = open(template).read()
    out = (tmpl
        .replace("/*__MOLSTAR_CSS__*/", open(molstar_css).read())
        .replace("/*__APP_CSS__*/", open(app_css).read())
        .replace("__INDEX_B64__", gz_b64(json.dumps(index_slim, separators=(",", ":"))))
        .replace("/*__PAKO_JS__*/", open(pako_js).read())
        .replace("__EMBED_JSON__", json.dumps(embed, separators=(",", ":")))
        .replace("/*__MOLSTAR_JS__*/", open(molstar_js).read())
        .replace("/*__APP_JS__*/", open(app_js).read()))
    os.makedirs(os.path.dirname(os.path.abspath(html_out)), exist_ok=True)
    with open(html_out, "w") as fh:
        fh.write(out)
    size = os.path.getsize(html_out)
    print(f"Wrote {html_out}: {size/1e6:.1f} MB, {len(index_slim['transcripts'])} transcripts")
    return html_out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="pipeline output dir (has data/ and structures/)")
    ap.add_argument("--molstar", required=True)
    ap.add_argument("--molstar-css", required=True)
    ap.add_argument("--pako", required=True)
    ap.add_argument("--template", required=True)
    ap.add_argument("--app-js", required=True)
    ap.add_argument("--app-css", required=True)
    ap.add_argument("--html", required=True)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    build_html(a.out, a.molstar, a.molstar_css, a.pako, a.template,
               a.app_js, a.app_css, a.html, a.limit)
