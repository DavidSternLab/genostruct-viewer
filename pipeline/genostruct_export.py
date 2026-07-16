"""
genostruct_export.py — GenBank export of genomic regions annotated with
colored protein structural features.

Given the pipeline output (per-transcript JSON records + the genome FASTA), emit
a GenBank file for either:
  * a single transcript's locus (its span +/- flank), or
  * an arbitrary genomic region (scaffold:start-end), including every transcript
    record that overlaps it.

The GenBank feature table includes, per transcript:
  gene, mRNA (exon join), CDS (join, phase-correct), and one misc_feature per
  structural element (alpha-helix / beta-strand) whose location is the exact
  genomic base ranges encoding that element (split across exons via join()),
  carrying /label, /note and a /color qualifier (hex) matching the viewer.

Colors are written both as a plain `/color="#rrggbb"` qualifier (read by many
genome/annotation tools) and as the Artemis/Geneious-style `/colour="R G B"`
qualifier for maximum interoperability.

Requires: biopython.
"""
from __future__ import annotations
import os, json, glob, argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio import SeqIO

# reuse FASTA loader semantics (standalone to avoid hard import dependency)
def load_fasta(path):
    d = {}; cur = None; seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None: d[cur] = "".join(seq)
                cur = line[1:].split()[0].strip(); seq = []
            else: seq.append(line.strip())
        if cur is not None: d[cur] = "".join(seq)
    return d

def _hex_to_rgb_str(hx):
    hx = hx.lstrip("#")
    return f"{int(hx[0:2],16)} {int(hx[2:4],16)} {int(hx[4:6],16)}"

def load_records(data_dir):
    recs = {}
    for f in glob.glob(os.path.join(data_dir, "*.json")):
        if os.path.basename(f) == "index.json": continue
        r = json.load(open(f))
        recs[r["transcript_id"]] = r
    return recs

def _feat_location(ranges, region_start, strand):
    """ranges: list of [gstart,gend] 1-based inclusive genomic. region_start: 1-based
    genomic coordinate of region base 1. Returns a Bio FeatureLocation/CompoundLocation
    in 0-based region coordinates with the given strand (+1/-1), or None if empty."""
    parts = []
    for a, b in sorted(ranges):
        s = a - region_start          # 0-based start within region
        e = b - region_start + 1      # exclusive end
        parts.append(FeatureLocation(s, e, strand=strand))
    if not parts:
        return None
    if len(parts) == 1:
        return parts[0]
    # For minus strand, CompoundLocation order should be 5'->3' (descending genomic);
    # BioPython writes join() and handles complement on output when strand=-1.
    if strand == -1:
        parts = parts[::-1]
    return CompoundLocation(parts)

def build_region_record(scaffold, region_start, region_end, genome, recs,
                        species="", assembly=""):
    """
    region_start/region_end: 1-based inclusive genomic coordinates.
    Returns a Bio.SeqRecord for the region with all overlapping annotations.
    """
    seq = genome[scaffold][region_start - 1:region_end]
    rec = SeqRecord(Seq(seq), id=scaffold,
                    name=(scaffold[:16] or "region"),
                    description=f"{species} {scaffold}:{region_start}-{region_end} "
                                f"structural-feature annotation".strip())
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "linear"
    if assembly: rec.annotations["comment"] = f"Assembly: {assembly}"
    rec.annotations["source"] = species or "genostruct"

    # source feature
    src = SeqFeature(FeatureLocation(0, len(seq), strand=1), type="source")
    src.qualifiers["organism"] = [species or "unknown"]
    src.qualifiers["note"] = [f"{scaffold}:{region_start}-{region_end}"]
    rec.features.append(src)

    for tid, r in recs.items():
        if r["scaffold"] != scaffold: continue
        if not r["span"]: continue
        gs, ge = r["span"]
        if ge < region_start or gs > region_end: continue   # no overlap
        strand = 1 if r["strand"] == "+" else -1
        fully_contained = (gs >= region_start and ge <= region_end)

        # gene (clipped to region if it runs off an edge)
        gf = SeqFeature(FeatureLocation(max(0, gs - region_start),
                                        min(len(seq), ge - region_start + 1),
                                        strand=strand), type="gene")
        gf.qualifiers["gene"] = [r["gene_id"]]
        gf.qualifiers["note"] = [f"transcript {tid}"]
        if not fully_contained:
            # Detailed feature coordinates would fall outside the exported window,
            # so we emit only the clipped gene box and flag it partial rather than
            # writing sequence-incorrect CDS/element features.
            gf.qualifiers["note"].append("partial: transcript extends beyond exported region")
            rec.features.append(gf)
            continue
        rec.features.append(gf)

        # mRNA (exon join)
        if r["exons"]:
            mloc = _feat_location(r["exons"], region_start, strand)
            mf = SeqFeature(mloc, type="mRNA")
            mf.qualifiers["gene"] = [r["gene_id"]]
            mf.qualifiers["note"] = [tid]
            rec.features.append(mf)

        # CDS (join, phase-correct)
        if r["cds"]:
            cloc = _feat_location([[c[0], c[1]] for c in r["cds"]], region_start, strand)
            cf = SeqFeature(cloc, type="CDS")
            cf.qualifiers["gene"] = [r["gene_id"]]
            cf.qualifiers["note"] = [f"{tid}; {r['n_helix']} helices, {r['n_strand']} strands"]
            cf.qualifiers["translation"] = [r["protein_sequence"]]
            cf.qualifiers["transl_table"] = [1]
            rec.features.append(cf)

        # structural elements
        for el in r["elements"]:
            eloc = _feat_location(el["genome_ranges"], region_start, strand)
            if eloc is None:
                continue   # element has no CDS-mappable genomic range; skip
            ef = SeqFeature(eloc, type="misc_feature")
            ef.qualifiers["label"] = [f"{tid}:{el['id']}"]
            ef.qualifiers["note"] = [
                f"{el['kind']} {el['id']} of {tid} "
                f"(protein residues {el['prot_start']}-{el['prot_end']})"]
            ef.qualifiers["color"] = [el["color"]]                    # hex, viewer-matched
            ef.qualifiers["colour"] = [_hex_to_rgb_str(el["color"])]  # Artemis/Geneious R G B
            ef.qualifiers["structure_type"] = [el["kind"]]
            rec.features.append(ef)

    return rec

def export_transcript(tid, genome, recs, out_path, flank=0, species="", assembly=""):
    r = recs[tid]
    gs, ge = r["span"]
    region_start = max(1, gs - flank)
    region_end = min(len(genome[r["scaffold"]]), ge + flank)
    rec = build_region_record(r["scaffold"], region_start, region_end, genome, recs,
                              species=species, assembly=assembly)
    rec.description = f"{species} {tid} locus".strip()
    SeqIO.write(rec, out_path, "genbank")
    return out_path, region_start, region_end

def export_region(scaffold, start, end, genome, recs, out_path, species="", assembly=""):
    rec = build_region_record(scaffold, start, end, genome, recs,
                              species=species, assembly=assembly)
    SeqIO.write(rec, out_path, "genbank")
    return out_path

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="GenBank export of genome regions with colored structural features")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-dir", required=True, help="viewer/data dir with per-transcript JSON")
    ap.add_argument("--out", required=True, help="output .gb path")
    ap.add_argument("--species", default="")
    ap.add_argument("--assembly", default="")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--transcript", help="transcript id to export (its locus)")
    g.add_argument("--region", help="scaffold:start-end (1-based inclusive)")
    ap.add_argument("--flank", type=int, default=0, help="flank bp for --transcript")
    a = ap.parse_args()

    genome = load_fasta(a.genome)
    recs = load_records(a.data_dir)
    if a.transcript:
        p, s, e = export_transcript(a.transcript, genome, recs, a.out,
                                    flank=a.flank, species=a.species, assembly=a.assembly)
        print(f"Wrote {p}  ({a.transcript}, region {s}-{e})")
    else:
        sc, rng = a.region.split(":"); s, e = rng.split("-")
        export_region(sc, int(s), int(e), genome, recs, a.out,
                      species=a.species, assembly=a.assembly)
        print(f"Wrote {a.out}  ({a.region})")
