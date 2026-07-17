# Genostruct — integrated genome · gene · protein · structure viewer

A pipeline and self-contained web viewer that link **protein 3D structure**,
**protein sequence**, and **the underlying genome** through a shared, per-element
color code. For each protein that has a structure model, every secondary-structure
element (α-helix / β-strand) is detected, assigned a unique color, and that same
color is used to paint:

1. the residues in the **3D structure** (Mol*),
2. the residues in the **protein sequence track**, and
3. the exact **genomic bases** (CDS codons, split across exons) that encode the
   element, in a purpose-built **genome track**.

The result: you can see a helix in 3D, find the same-colored stretch in the
sequence, and read off the exon(s) and genomic coordinates that encode it — with
introns and strand correctly represented.

---

## Contents

| File | What it is |
|---|---|
| `pipeline/genostruct.py` | The data pipeline (generic; no organism/file names hardcoded). |
| `pipeline/build_viewer.py` | Assembles the self-contained viewer HTML from pipeline output. |
| `pipeline/genostruct_export.py` | GenBank export of any transcript locus or genomic region, annotated with colored structural features. |
| `viewer_src/` | Front-end source: `app.js`, `app.css`, `template.html`. |
| `genostruct_viewer.html` | **The deliverable** — one self-contained file, open it in a browser. |
| `viewer/data/*.json`, `viewer/data/index.json` | Precomputed per-transcript records (SSE, residue→genome maps, colors). |
| `exports/*.gb` | Example GenBank exports. |
| `validation_colorlink.png` | Static proof that sequence colors map to the correct genomic segments. |

---

## Quick start

**Just view the results:** open `genostruct_viewer.html` in any modern browser.
No server, no network — Mol*, the decompressor, all structure coordinates and
per-transcript data are embedded (gzip + base64, ~10.7 MB for 187 transcripts).

Pick a transcript from the dropdown. The four stacked panels are color-linked:
hovering or clicking an element in any panel highlights it everywhere.

### Viewer controls

- **3D structure (Mol\*)** — rotate (drag), zoom (scroll), pan (right-drag). Native Mol* interaction.
  - **Show pLDDT / Show element colors** — toggle the cartoon between the per-element color scheme (matched to the sequence and genome tracks) and the AF2 pLDDT confidence ramp (from the model B-factors).
- **Protein sequence** — zoomable: use the **− / reset / +** buttons, or **Ctrl/⌘ + scroll** over the sequence.
  - **Drag-select a residue range** — dims the rest of the structure, recenters the 3D camera on the selection, and zooms the genome track to the encoding bases. Double-click the sequence to clear.
- **Genome track** — **scroll to zoom** (centered on the cursor), **drag to pan**, **double-click to reset** to the full locus. A readout shows current zoom factor and window size in bp.
  - **Show genome 5'→3' / Show protein N→C** — toggle the genome-track orientation. By default the track is drawn in the protein's N→C direction (reversed for minus-strand genes); the toggle switches to native 5'→3'.
- **Structural elements (legend or sequence)** — click an element (a legend chip, a colored sequence block, or a genome feature) to **focus it in 3D** (that element stays vivid, the rest dims), **recenter the camera** on it, and **zoom the genome track** to its bases. Click again or click empty space to clear. Hovering an element highlights it across all three panels.
- **Search box** — type a gene number to filter the transcript list; the full numerically-ordered dropdown stays available.
- **Export GenBank (transcript / viewed region)** — downloads an annotated GenBank record; the viewed-region export includes the full sequence of the current genome window (see note below).

---

## Running the pipeline on your own data

The pipeline is generic — nothing about the organism, file names, or gene IDs is
hardcoded. It needs four inputs:

- a **genome FASTA** (`*.fa`),
- a **GFF3** annotation (genes / mRNA / exon / CDS with phase),
- a **protein FASTA** (`pep.fa`),
- a folder of **PDB** structure models (e.g. AlphaFold predictions), for any
  subset of the proteins.

```python
import sys; sys.path.insert(0, "pipeline")
import genostruct

genostruct.build(
    genome_fa   = "path/to/genome.fa",
    gff3        = "path/to/annotation.gff3",
    pep_fa      = "path/to/pep.fa",
    pdb_dir     = "path/to/pdbs",
    out_dir     = "viewer",
    # id_regex extracts the transcript id from each PDB filename; the default
    # tries several strategies and falls back to matching known transcript ids.
    id_regex    = r"(.+)_ranked_0",   # example for *_ranked_0.pdb
)
```

This writes `viewer/data/<transcript>.json`, `viewer/data/index.json`, and copies
matched PDBs into `viewer/structures/`.

Then build the standalone HTML:

```bash
python pipeline/build_viewer.py \
  --out viewer \
  --molstar molstar.js --molstar-css molstar.css --pako pako_inflate.min.js \
  --template viewer_src/template.html \
  --app-js viewer_src/app.js --app-css viewer_src/app.css \
  --html viewer/genostruct_viewer.html
```

(`molstar.js`, `molstar.css`, `pako_inflate.min.js` are bundled here; they are the
only third-party assets and are embedded into the output.)

---

## GenBank export

`genostruct_export.py` produces a GenBank record — **with genomic sequence** — for
either a single transcript's locus or an arbitrary region. Structural elements are
written as `misc_feature`s whose location is the exact codon ranges (as
`join(...)`, `complement(...)` for minus strand), carrying:

- `/label` = `transcript:element` (e.g. `g113516.t1-00001:H2`),
- `/note` = element kind + protein residue range,
- `/color` = the viewer hex color (e.g. `#e6194B`),
- `/colour` = Artemis/Geneious-style `R G B` (for tools that read that convention),
- `/structure_type` = `helix` or `strand`.

```bash
# a single transcript's locus, with 200 bp flanks
python pipeline/genostruct_export.py \
  --genome genome.fa --data-dir viewer/data \
  --transcript g113516.t1-00001 --flank 200 \
  --species "Acyrthosiphon pisum" --assembly JIC1_v1.0 \
  --out exports/g113516_locus.gb

# an arbitrary region (all overlapping transcripts included)
python pipeline/genostruct_export.py \
  --genome genome.fa --data-dir viewer/data \
  --region scaffold_2:89950000-90000000 \
  --out exports/region.gb
```

The `Export GenBank` button in the viewer produces a lightweight **feature-table-only**
version (no sequence — kept out of the browser to hold file size down). For a
complete record with sequence, use `genostruct_export.py`.

---

## Method notes

- **Secondary structure** is computed from the model's Cα coordinates with
  biotite's `annotate_sse` (the P-SEA algorithm), then segmented into contiguous
  helix/strand elements. Each element gets a unique color from a warm palette
  (helices) or cool palette (strands), cycled if there are many.
- **Residue → genome mapping** walks each CDS block with its GFF phase, in
  translation order, on the correct strand, so every protein residue maps to its
  three genomic codon positions — correctly handling introns and codons that
  straddle exon boundaries.
- **Model ↔ protein alignment.** AlphaFold models here were predicted from
  signal-peptide-trimmed sequences and are numbered from 1, so each model is
  aligned to its full protein to recover an N-terminal offset (exact substring →
  gapless scan → local alignment fallback). Colors are indexed by full-protein
  residue, and the Mol* color theme applies the offset so 3D coloring lines up
  with the sequence and genome.

---

## Validation

- **Coordinate correctness:** for all 187 transcripts, the exported CDS translates
  back to the exact deposited protein (**187/187**), and extracting + translating
  every structural-element `misc_feature` reproduces its protein residues exactly
  (**1113/1113** elements, **187/187** transcripts fully correct — no CDS
  mismatches, no out-of-range elements, no empty ranges).
- Element residue ranges are clamped to the coding protein length in the pipeline,
  so models that are longer than their protein (from a gapped local alignment) can
  no longer map an element past the C-terminus into the stop codon or out-of-range
  codon positions.
- `validation_colorlink.png` visually confirms sequence-color → genome-segment
  correspondence, including multi-exon elements on the minus strand.

### Known limitations

- **3 PDB models** in the source set correspond to stale isoform IDs absent from the
  protein FASTA and are excluded (documented in the pipeline output warnings).
- **10 models** differ from their deposited protein by internal indels and are placed
  by local alignment (identity ≥ 0.6); their per-residue mapping is still validated
  correct by the CDS/element round-trip above.
- The viewer runs entirely offline in a single self-contained HTML file. If it is
  opened inside a **sandboxed iframe with a restrictive Content-Security-Policy**
  (some hosted environments do this), the file-download of a GenBank export may be
  blocked by the sandbox; a copy-to-clipboard / open-in-new-tab fallback is provided
  for that case. Opened directly in a browser (double-clicking the file, or served
  from a normal web server), downloads work normally.
