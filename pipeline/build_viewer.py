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
import os, glob, json, gzip, base64, argparse

def gz_b64(text: str) -> str:
    if isinstance(text, str): text = text.encode("utf-8")
    return base64.b64encode(gzip.compress(text, 9)).decode("ascii")

# Standard 3-letter residue names for HELIX/SHEET records (fallback ALA).
_AA3 = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
        "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
        "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}

def _pdb_resnames(pdb_text):
    """auth_seq_id -> 3-letter residue name, from CA ATOM records (chain-agnostic)."""
    res, chain = {}, {}
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            try:
                n = int(line[22:26])
            except ValueError:
                continue
            res.setdefault(n, line[17:20].strip() or "ALA")
            chain.setdefault(n, (line[21:22].strip() or "A"))
    return res, chain

def _helix_record(serial, hid, init_resn, init_chain, init_seq,
                  end_resn, end_chain, end_seq, length):
    """Build a spec-column PDB HELIX record (offsets verified against Mol*'s parser).
    cols: serNum 7(3), helixID 11(3), initResName 15(3), initChain 19(1),
          initSeq 21(4), endResName 27(3), endChain 31(1), endSeq 33(4),
          helixClass 38(2), length 71(5)."""
    s = list(" " * 80)
    def put(text, start, width, right=False):
        text = str(text)[:width]
        text = text.rjust(width) if right else text.ljust(width)
        s[start:start+width] = list(text)
    put("HELIX", 0, 6)
    put(serial, 7, 3, right=True)
    put(hid, 11, 3, right=True)
    put(init_resn, 15, 3); put(init_chain, 19, 1)
    put(init_seq, 21, 4, right=True)
    put(end_resn, 27, 3); put(end_chain, 31, 1)
    put(end_seq, 33, 4, right=True)
    put("1", 38, 2, right=True)          # helixClass 1 = right-handed alpha
    put(length, 71, 5, right=True)
    return "".join(s).rstrip()

def _sheet_record(strand, sid, nstr, init_resn, init_chain, init_seq,
                  end_resn, end_chain, end_seq):
    """Build a spec-column PDB SHEET record (offsets verified against Mol*'s parser).
    cols: strand 7(3), sheetID 11(3), numStrands 14(2), initResName 17(3),
          initChain 21(1), initSeq 22(4), endResName 28(3), endChain 32(1),
          endSeq 33(4), sense 38(2)."""
    s = list(" " * 80)
    def put(text, start, width, right=False):
        text = str(text)[:width]
        text = text.rjust(width) if right else text.ljust(width)
        s[start:start+width] = list(text)
    put("SHEET", 0, 6)
    put(strand, 7, 3, right=True)
    put(sid, 11, 3)
    put(nstr, 14, 2, right=True)
    put(init_resn, 17, 3); put(init_chain, 21, 1)
    put(init_seq, 22, 4, right=True)
    put(end_resn, 28, 3); put(end_chain, 32, 1)
    put(end_seq, 33, 4, right=True)
    put("0", 38, 2, right=True)          # sense 0 for the first strand of a sheet
    return "".join(s).rstrip()

def pdb_with_ss(pdb_text, rec):
    """Prepend HELIX/SHEET records derived from the pipeline's SSE elements so
    Mol*'s 'auto' secondary-structure reads them from the model instead of
    running DSSP (whose WASM is fetched from a data: URL that the viewer's CSP
    blocks). Element prot_start/prot_end are full-protein 1-based; model auth =
    prot - offset (AF2 PDBs are numbered from 1)."""
    off = (rec.get("model") or {}).get("offset", 0) or 0
    res, chain = _pdb_resnames(pdb_text)
    if not res:
        return pdb_text
    lo_auth, hi_auth = min(res), max(res)
    hlines, slines, hser, sser = [], [], 0, 0
    for el in rec.get("elements", []):
        a = el["prot_start"] - off
        b = el["prot_end"] - off
        if b < a:
            continue
        a = max(lo_auth, a); b = min(hi_auth, b)   # clamp into the model
        if b < a:
            continue
        init_resn = res.get(a, "ALA"); end_resn = res.get(b, "ALA")
        init_chain = chain.get(a, "A"); end_chain = chain.get(b, "A")
        if el["kind"] == "helix":
            hser += 1
            hlines.append(_helix_record(hser, hser, init_resn, init_chain, a,
                                        end_resn, end_chain, b, b - a + 1))
        else:
            sser += 1
            slines.append(_sheet_record(sser, "A", 1, init_resn, init_chain, a,
                                        end_resn, end_chain, b))
    header = "\n".join(hlines + slines)
    return (header + "\n" + pdb_text) if header else pdb_text

def build_html(out_dir, molstar_js, molstar_css, pako_js,
               template, app_js, app_css, html_out, limit=None):
    data_dir = os.path.join(out_dir, "data")
    struct_dir = os.path.join(out_dir, "structures")
    index = json.load(open(os.path.join(data_dir, "index.json")))
    txs = index["transcripts"]
    if limit: txs = txs[:limit]

    embed = {"data": {}, "pdb": {}}
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
            pdb_text = open(pdb_path, "r").read()
            # inject HELIX/SHEET so Mol* uses model SS (avoids the DSSP WASM fetch)
            pdb_text = pdb_with_ss(pdb_text, rec)
            embed["pdb"][tid] = gz_b64(pdb_text)

    # trim index to the transcripts we actually embedded
    index_slim = {"transcripts": [t for t in txs if t["transcript_id"] in embed["pdb"]],
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
