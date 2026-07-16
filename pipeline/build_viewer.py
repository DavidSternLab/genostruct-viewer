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
        embed["data"][tid] = gz_b64(open(rec_path, "rb").read())
        # find pdb (safe-name)
        safe = "".join(c if (c.isalnum() or c in "._-") else "_" for c in tid)
        pdb_path = os.path.join(struct_dir, safe + ".pdb")
        if not os.path.exists(pdb_path):
            # fallback: any pdb whose stem contains tid
            cands = glob.glob(os.path.join(struct_dir, "*.pdb"))
            cands = [c for c in cands if safe in os.path.basename(c)]
            pdb_path = cands[0] if cands else None
        if pdb_path and os.path.exists(pdb_path):
            embed["pdb"][tid] = gz_b64(open(pdb_path, "rb").read())

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
