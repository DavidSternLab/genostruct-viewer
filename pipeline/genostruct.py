"""
genostruct — Integrated genome / gene / protein / structure pipeline.

Generic, input-agnostic pipeline that, given
    * a genome FASTA
    * a GFF3 annotation (with mRNA / exon / CDS features)
    * a peptide FASTA (protein sequences keyed by transcript ID)
    * a directory of PDB structures (one model per transcript, e.g. AlphaFold2)

produces, for every transcript that has a structure model:
    * per-residue secondary-structure assignment (P-SEA, CA-based)
    * segmentation into discrete structural elements (helix-1, strand-1, ...)
    * a unique color per element (cyclic palette), consistent across panels
    * a residue -> genomic-coordinate map (via CDS + phase + strand)
    * per-element genomic base ranges (split across exons where needed)

Nothing about the organism, file names, or gene IDs is hardcoded.
The PDB<->transcript mapping is auto-detected (see map_pdbs_to_transcripts).

Outputs one compact JSON per transcript plus an index.json, suitable for the
self-contained HTML viewer, and supports GenBank export of any genomic region.
"""
from __future__ import annotations
import os, re, glob, json, argparse
from collections import defaultdict
from dataclasses import dataclass, field, asdict

import numpy as np

# ----------------------------------------------------------------------------
# Genetic code (standard table 1)
# ----------------------------------------------------------------------------
CODON_TABLE = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
 'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
 'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','SEC':'U','PYL':'O','MSE':'M'}

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]

def translate(nt: str) -> str:
    return "".join(CODON_TABLE.get(nt[i:i+3].upper(), 'X') for i in range(0, len(nt) - 2, 3))

# ----------------------------------------------------------------------------
# FASTA
# ----------------------------------------------------------------------------
def load_fasta(path: str) -> dict:
    """Load a FASTA file; keys are the first whitespace-delimited token of each header."""
    d = {}; cur = None; seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None: d[cur] = "".join(seq)
                cur = line[1:].split()[0].strip(); seq = []
            else:
                seq.append(line.strip())
        if cur is not None: d[cur] = "".join(seq)
    return d

# ----------------------------------------------------------------------------
# GFF3
# ----------------------------------------------------------------------------
@dataclass
class Transcript:
    tid: str
    scaffold: str
    strand: str
    gene: str = ""
    cds: list = field(default_factory=list)    # list of (start, end, phase)   1-based inclusive
    exons: list = field(default_factory=list)  # list of (start, end)          1-based inclusive
    span: tuple = None                          # (start, end)

def parse_gff3(path: str):
    """Return {tid: Transcript}. Recognises mRNA/transcript + exon + CDS by Parent linkage."""
    tx = {}
    def get(attr, key):
        m = re.search(rf"{key}=([^;]+)", attr); return m.group(1) if m else None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            chrom, src, ftype, start, end, score, strand, phase, attr = f
            start, end = int(start), int(end)
            if ftype in ("mRNA", "transcript"):
                tid = get(attr, "ID")
                if tid is None: continue
                t = tx.setdefault(tid, Transcript(tid, chrom, strand))
                t.scaffold, t.strand = chrom, strand
                t.gene = get(attr, "Parent") or ""
                t.span = (start, end)
            elif ftype == "CDS":
                pid = get(attr, "Parent")
                if pid is None: continue
                ph = int(phase) if phase.isdigit() else 0
                for p in pid.split(","):
                    t = tx.setdefault(p, Transcript(p, chrom, strand))
                    t.scaffold, t.strand = chrom, strand
                    t.cds.append((start, end, ph))
            elif ftype == "exon":
                pid = get(attr, "Parent")
                if pid is None: continue
                for p in pid.split(","):
                    t = tx.setdefault(p, Transcript(p, chrom, strand))
                    t.scaffold, t.strand = chrom, strand
                    t.exons.append((start, end))
    # sort features
    for t in tx.values():
        t.cds.sort(); t.exons.sort()
    return tx

# ----------------------------------------------------------------------------
# Residue -> genome coordinate mapping
# ----------------------------------------------------------------------------
def cds_codon_map(t: Transcript):
    """
    Build the ordered list of coding genomic positions (5'->3' along the mRNA),
    honoring strand. Returns list of 1-based genomic coordinates, one per coding
    nucleotide, in translation order. Also returns the assembled coding nt string.
    """
    if t.strand == '+':
        blocks = sorted(t.cds, key=lambda x: x[0])
        positions = []
        for a, b, _ in blocks:
            positions.extend(range(a, b + 1))
    else:
        blocks = sorted(t.cds, key=lambda x: x[0], reverse=True)
        positions = []
        for a, b, _ in blocks:
            positions.extend(range(b, a - 1, -1))   # descending on '-' strand
    return positions

def coding_nt(t: Transcript, genome: dict):
    sc = genome[t.scaffold]
    if t.strand == '+':
        blocks = sorted(t.cds, key=lambda x: x[0])
        return "".join(sc[a - 1:b] for a, b, _ in blocks)
    else:
        blocks = sorted(t.cds, key=lambda x: x[0], reverse=True)
        return "".join(revcomp(sc[a - 1:b]) for a, b, _ in blocks)

def residue_genome_map(t: Transcript, genome: dict):
    """
    Map protein residue index (1-based, over the FULL translated CDS) to the
    genomic coordinates of its codon.

    Returns:
      res2codon: list where res2codon[i] = [g1,g2,g3] genomic (1-based) positions
                 of codon for residue i+1, in translation (5'->3') order.
    Uses CDS phase of the first block to skip any leading partial codon.
    """
    positions = cds_codon_map(t)
    # leading phase: number of bases to skip at the 5' end so frame starts at a codon
    first_phase = t.cds[0][2] if t.strand == '+' else sorted(t.cds, key=lambda x: x[0], reverse=True)[0][2]
    start = first_phase
    res2codon = []
    i = start
    while i + 3 <= len(positions):
        res2codon.append(positions[i:i + 3])
        i += 3
    return res2codon

# ----------------------------------------------------------------------------
# PDB parsing (sequence + residue ids) — no external deps
# ----------------------------------------------------------------------------
def pdb_sequence(path: str):
    """Return (seq, resids, bfactors) from the first model / chain of a PDB.
    seq: one-letter string; resids: list of residue numbers; bfactors: mean CA plddt per residue."""
    seq = []; resids = []; bfac = []
    seen = set(); last = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"): break
            if not line.startswith("ATOM"): continue
            resnum = int(line[22:26]); resname = line[17:20].strip()
            key = (resnum, line[21])           # (resseq, chain)
            if key in seen:
                if line[12:16].strip() == "CA":
                    bfac[-1] = float(line[60:66])
                continue
            seen.add(key)
            seq.append(THREE2ONE.get(resname, 'X'))
            resids.append(resnum)
            bfac.append(float(line[60:66]) if line[12:16].strip() == "CA" else 0.0)
    return "".join(seq), resids, bfac

# ----------------------------------------------------------------------------
# PDB <-> transcript auto-mapping
# ----------------------------------------------------------------------------
def map_pdbs_to_transcripts(pdb_dir: str, transcript_ids, pep: dict,
                            id_regex: str | None = None):
    """
    Auto-detect which transcript each PDB corresponds to, generically.

    Strategy (in order):
      1. If id_regex given, use its first capture group as the transcript id.
      2. Try to find any known transcript id as a substring of the file stem
         (longest match wins) — robust to prefixes/suffixes like
         'Species_<tid>_ranked_0.pdb'.
      3. Fall back to the file stem itself.

    Returns {tid: pdb_path} for those that resolve to a known transcript.
    """
    pdbs = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))) + \
           sorted(glob.glob(os.path.join(pdb_dir, "*.cif")))
    tid_set = set(transcript_ids)
    # sort ids by length desc for greedy longest-substring matching
    ids_by_len = sorted(tid_set, key=len, reverse=True)
    mapping = {}
    unresolved = []
    for p in pdbs:
        stem = os.path.splitext(os.path.basename(p))[0]
        tid = None
        if id_regex:
            m = re.search(id_regex, stem)
            if m: tid = m.group(1)
        if tid is None or tid not in tid_set:
            # greedy substring search
            for cand in ids_by_len:
                if cand in stem:
                    tid = cand; break
        if tid is None or tid not in tid_set:
            unresolved.append(os.path.basename(p)); continue
        # keep first (or prefer one whose pep matches better later)
        mapping.setdefault(tid, p)
    return mapping, unresolved

# ----------------------------------------------------------------------------
# Align the (possibly signal-peptide-trimmed) model sequence to the full protein
# ----------------------------------------------------------------------------
def align_offset(model_seq: str, full_seq: str):
    """
    Find the offset such that model residue k (0-based) corresponds to
    full protein residue (offset + k), 0-based.

    Returns (offset, method, identity):
      * exact substring  -> ('substring', 1.0)
      * else local (Smith-Waterman-ish) alignment; offset = full index aligned
        to model position 0, identity = matched/len(model).
    Returns offset=None if no reasonable alignment (identity < 0.6).
    """
    full = full_seq.rstrip('*')
    idx = full.find(model_seq)
    if idx >= 0:
        return idx, "substring", 1.0
    # try trimming trailing X / mismatched ends via a simple gapless best-offset scan
    best_off, best_score = None, -1
    m = len(model_seq)
    for off in range(0, len(full) - m + 1):
        window = full[off:off + m]
        score = sum(1 for a, b in zip(model_seq, window) if a == b)
        if score > best_score:
            best_score, best_off = score, off
    if best_off is not None and best_score / max(m, 1) >= 0.6:
        return best_off, "gapless", best_score / m
    # last resort: gapped local alignment (needleman on model vs full, keep offset)
    off, ident = _local_align_offset(model_seq, full)
    if off is not None and ident >= 0.6:
        return off, "local", ident
    return None, "none", best_score / max(m, 1)

def _local_align_offset(model_seq: str, full: str):
    """Lightweight Smith-Waterman to recover the full-protein index aligned to
    model position 0. Returns (offset, identity)."""
    try:
        from Bio import pairwise2
    except Exception:
        return None, 0.0
    aln = pairwise2.align.localms(full, model_seq, 2, -1, -5, -0.5,
                                  one_alignment_only=True)
    if not aln: return None, 0.0
    a = aln[0]
    # walk to first aligned (non-gap in model) column
    fi = 0; off = None; matches = 0; mlen = 0
    for ca, cb in zip(a.seqA, a.seqB):
        if cb != '-':
            if off is None and ca != '-':
                off = fi
            if ca != '-':
                mlen += 1
                if ca == cb: matches += 1
        if ca != '-': fi += 1
    ident = matches / max(len(model_seq), 1)
    return off, ident

# ----------------------------------------------------------------------------
# Secondary structure -> discrete elements + colors
# ----------------------------------------------------------------------------
# A perceptually-spread cyclic palette (helices warm, strands cool families).
HELIX_PALETTE = ["#e6194B","#f58231","#ffe119","#bfef45","#fabed4","#9A6324",
                 "#800000","#aaffc3","#ff4d6d","#ffa600"]
STRAND_PALETTE = ["#4363d8","#42d4f4","#911eb4","#469990","#000075","#3cb44b",
                  "#008080","#6a5acd","#00bfff","#7f00ff"]

def segment_sse(sse_labels):
    """
    Given a per-residue SSE array ('a'/'b'/'c'/''), produce ordered discrete
    elements. Returns list of dicts:
      {kind:'helix'|'strand', index:int(1-based within kind),
       start:int, end:int}  (start/end are 0-based residue indices into model)
    """
    elems = []
    n = len(sse_labels)
    i = 0
    h = s = 0
    while i < n:
        lab = sse_labels[i]
        if lab in ('a', 'b'):
            j = i
            while j < n and sse_labels[j] == lab:
                j += 1
            if lab == 'a':
                h += 1; kind = 'helix'; idx = h
            else:
                s += 1; kind = 'strand'; idx = s
            elems.append({"kind": kind, "index": idx, "start": i, "end": j - 1})
            i = j
        else:
            i += 1
    return elems

def color_elements(elems):
    """Assign a unique color per element from the cyclic palettes."""
    hc = sc = 0
    for e in elems:
        if e["kind"] == "helix":
            e["color"] = HELIX_PALETTE[hc % len(HELIX_PALETTE)]; hc += 1
        else:
            e["color"] = STRAND_PALETTE[sc % len(STRAND_PALETTE)]; sc += 1
    return elems

# ----------------------------------------------------------------------------
# Element -> genomic ranges
# ----------------------------------------------------------------------------
def _merge_positions_to_ranges(positions):
    """Collapse a sorted list of ints into contiguous [start,end] inclusive ranges."""
    if not positions: return []
    positions = sorted(positions)
    ranges = []
    a = b = positions[0]
    for p in positions[1:]:
        if p == b + 1:
            b = p
        else:
            ranges.append([a, b]); a = b = p
    ranges.append([a, b])
    return ranges

def element_genome_ranges(elem, model_offset, res2codon, prot_len=None):
    """
    For one structural element, return the genomic base ranges it is encoded by.
    model residue e['start']..e['end'] (0-based) -> full protein residue
    (model_offset + k) -> codon genomic positions.

    `prot_len` (number of coding residues, excluding the stop codon) bounds the
    mapping: a model that is longer than its protein (e.g. from a gapped local
    alignment) can produce elements whose residue indices run past the protein
    C-terminus; without this bound they would incorrectly pull in the stop-codon
    or out-of-range codon positions. Residues at/after prot_len are dropped.
    Returns list of [gstart, gend] inclusive (merged, ascending).
    """
    upper = len(res2codon) if prot_len is None else min(prot_len, len(res2codon))
    gpos = []
    for k in range(elem["start"], elem["end"] + 1):
        full_res = model_offset + k          # 0-based index into full protein
        if 0 <= full_res < upper:
            gpos.extend(res2codon[full_res])
    return _merge_positions_to_ranges(gpos)

def annotate_sse_from_pdb(pdb_path):
    """Compute P-SEA SSE labels aligned to the model's residue order.
    Returns (labels list, resids list, per-residue plddt list, model one-letter seq)."""
    import biotite.structure as struc
    from biotite.structure.io.pdb import PDBFile
    pf = PDBFile.read(pdb_path)
    arr = pf.get_structure(model=1, extra_fields=["b_factor"])
    # first chain only
    chains = np.unique(arr.chain_id)
    arr = arr[arr.chain_id == chains[0]]
    sse = struc.annotate_sse(arr)            # one per residue, in residue order
    # residue ids & sequence in same order
    res_ids, res_names = struc.get_residues(arr)
    seq = "".join(THREE2ONE.get(r, 'X') for r in res_names)
    # per-residue mean CA b-factor (plddt)
    plddt = []
    ca = arr[arr.atom_name == "CA"]
    ca_map = {int(rid): float(b) for rid, b in zip(ca.res_id, ca.b_factor)}
    for rid in res_ids:
        plddt.append(ca_map.get(int(rid), None))
    return list(sse), [int(x) for x in res_ids], plddt, seq

# ----------------------------------------------------------------------------
# Main build
# ----------------------------------------------------------------------------
def build_transcript_record(tid, t, pdb_path, full_prot, genome, copy_pdb_to=None):
    """
    Assemble the full per-transcript record. Returns (record_dict, warnings).
    """
    warns = []
    sse, resids, plddt, model_seq = annotate_sse_from_pdb(pdb_path)

    # align model to full protein (recover signal-peptide/offset)
    offset, method, ident = align_offset(model_seq, full_prot)
    if offset is None:
        warns.append(f"{tid}: model seq could not be aligned to protein (ident={ident:.2f})")
        offset = 0

    # residue -> genome codon map (full protein coordinates)
    res2codon = residue_genome_map(t, genome)
    # sanity: translated CDS should match full protein (trim stop)
    trans = translate(coding_nt(t, genome)).rstrip('*')
    cds_matches = (trans == full_prot.rstrip('*'))
    if not cds_matches:
        warns.append(f"{tid}: CDS translation != pep ({len(trans)} vs {len(full_prot)})")

    # SSE elements
    elems = color_elements(segment_sse(sse))

    # number of coding residues (excludes any trailing stop codon)
    prot_len = len(full_prot.rstrip('*'))

    # attach full-protein residue ranges + genomic ranges to each element.
    # prot_start/prot_end are clamped to [1, prot_len]: a model longer than its
    # protein (gapped local alignment) can otherwise emit residue indices past the
    # C-terminus, which would map into the stop codon / out-of-range positions.
    keep = []
    for e in elems:
        e["prot_start"] = max(1, offset + e["start"] + 1)          # 1-based, full protein
        e["prot_end"]   = min(prot_len, offset + e["end"] + 1)
        if e["prot_end"] < e["prot_start"]:
            continue   # element lies entirely outside the protein; drop it
        e["genome_ranges"] = element_genome_ranges(e, offset, res2codon, prot_len)
        e["id"] = f"{e['kind'][0].upper()}{e['index']}"   # H1, S1, ...
        keep.append(e)
    elems = keep

    # per-residue color track over the FULL protein (default coil = None)
    res_color = [None] * prot_len
    res_sse   = ['c'] * prot_len
    res_elem  = [None] * prot_len
    for e in elems:
        for pr in range(e["prot_start"], e["prot_end"] + 1):
            if 1 <= pr <= prot_len:
                res_color[pr - 1] = e["color"]
                res_sse[pr - 1] = 'a' if e["kind"] == 'helix' else 'b'
                res_elem[pr - 1] = e["id"]

    # per-residue genomic codon (for the full protein), for the genome track
    res_genome = []
    for i in range(prot_len):
        res_genome.append(res2codon[i] if i < len(res2codon) else None)

    pdb_ref = pdb_path
    if copy_pdb_to:
        import shutil
        os.makedirs(copy_pdb_to, exist_ok=True)
        dst = os.path.join(copy_pdb_to, f"{_safe(tid)}.pdb")
        shutil.copy(pdb_path, dst)
        pdb_ref = os.path.relpath(dst, os.path.dirname(copy_pdb_to.rstrip('/')))

    rec = {
        "transcript_id": tid,
        "gene_id": t.gene,
        "scaffold": t.scaffold,
        "strand": t.strand,
        "span": list(t.span) if t.span else None,
        # forward-strand genomic sequence over the locus span (1-based inclusive).
        # Lets the viewer emit GenBank with the actual sequence for any sub-window
        # without shipping the whole genome. Loci are small (~kb), gzips well.
        "locus_sequence": (str(genome[t.scaffold][t.span[0]-1:t.span[1]]).upper()
                            if t.span and t.scaffold in genome else None),
        "protein_length": prot_len,
        "protein_sequence": full_prot.rstrip('*'),
        "model": {
            "path": os.path.basename(pdb_ref),
            "offset": offset,               # model res 0 -> full protein res (0-based) offset
            "align_method": method,
            "align_identity": round(ident, 3),
            "model_length": len(model_seq),
            "plddt": [round(p, 1) if p is not None else None for p in plddt],
        },
        "cds_translation_ok": cds_matches,
        "exons": [list(e) for e in t.exons],
        "cds": [list(c) for c in t.cds],
        "elements": elems,                  # each: id,kind,index,color,prot_start,prot_end,genome_ranges
        "residue_color": res_color,         # len prot_len, hex or None
        "residue_sse": res_sse,             # len prot_len, a/b/c
        "residue_element": res_elem,        # len prot_len, element id or None
        "residue_genome": res_genome,       # len prot_len, [g1,g2,g3] or None
        "n_helix": sum(1 for e in elems if e["kind"] == "helix"),
        "n_strand": sum(1 for e in elems if e["kind"] == "strand"),
    }
    return rec, warns

def _safe(s):
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)

def build(genome_fa, gff3, pep_fa, pdb_dir, out_dir,
          id_regex=None, limit=None, verbose=True):
    """
    Run the full pipeline. Writes out_dir/data/<tid>.json and out_dir/data/index.json,
    and copies the matched PDBs into out_dir/structures/.
    Returns the index dict.
    """
    os.makedirs(out_dir, exist_ok=True)
    data_dir = os.path.join(out_dir, "data"); os.makedirs(data_dir, exist_ok=True)
    struct_dir = os.path.join(out_dir, "structures")

    if verbose: print("Loading genome ..."); 
    genome = load_fasta(genome_fa)
    if verbose: print(f"  {len(genome)} sequences")
    if verbose: print("Loading peptides ...")
    pep = load_fasta(pep_fa)
    if verbose: print(f"  {len(pep)} proteins")
    if verbose: print("Parsing GFF3 ...")
    tx = parse_gff3(gff3)
    if verbose: print(f"  {len(tx)} transcripts")

    if verbose: print("Mapping PDBs to transcripts ...")
    # only transcripts that have both CDS and a peptide are usable
    usable_ids = [tid for tid, t in tx.items() if t.cds and tid in pep]
    mapping, unresolved = map_pdbs_to_transcripts(pdb_dir, usable_ids, pep, id_regex)
    if verbose:
        print(f"  {len(mapping)} PDBs mapped, {len(unresolved)} unresolved")

    index = {"transcripts": [], "unresolved_pdbs": unresolved, "warnings": []}
    items = list(mapping.items())
    if limit: items = items[:limit]
    for n, (tid, pdb_path) in enumerate(sorted(items), 1):
        t = tx[tid]; full_prot = pep[tid]
        try:
            rec, warns = build_transcript_record(tid, t, pdb_path, full_prot, genome,
                                                 copy_pdb_to=struct_dir)
        except Exception as ex:
            index["warnings"].append(f"{tid}: FAILED {type(ex).__name__}: {ex}")
            continue
        index["warnings"].extend(warns)
        with open(os.path.join(data_dir, f"{_safe(tid)}.json"), "w") as fh:
            json.dump(rec, fh, separators=(",", ":"))
        index["transcripts"].append({
            "transcript_id": tid, "gene_id": t.gene, "scaffold": t.scaffold,
            "strand": t.strand, "span": list(t.span) if t.span else None,
            "protein_length": rec["protein_length"],
            "n_helix": rec["n_helix"], "n_strand": rec["n_strand"],
            "file": f"{_safe(tid)}.json",
            "cds_translation_ok": rec["cds_translation_ok"],
        })
        if verbose and n % 25 == 0: print(f"  processed {n}/{len(items)}")
    index["transcripts"].sort(key=lambda r: (r["scaffold"], r["span"][0] if r["span"] else 0))
    with open(os.path.join(data_dir, "index.json"), "w") as fh:
        json.dump(index, fh, separators=(",", ":"))
    if verbose:
        print(f"Done: {len(index['transcripts'])} transcript records -> {data_dir}")
        print(f"  {len(index['warnings'])} warnings")
    return index

# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Integrated genome/gene/protein/structure pipeline")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--pep", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--id-regex", default=None,
                    help="regex whose group(1) extracts the transcript id from a PDB filename stem")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    build(args.genome, args.gff, args.pep, args.pdb_dir, args.out,
          id_regex=args.id_regex, limit=args.limit)



