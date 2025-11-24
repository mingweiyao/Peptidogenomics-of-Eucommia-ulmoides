#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
start_codon_scanner.py

A lightweight pipeline to:
1) Filter peptide rows (proteins == 1 and type != "exon_diff") from an Excel file.
2) Extract original sequences and compute coverage per accession.
3) For each peptide, scan for candidate start codons (ATG, CTG, GTG, TTG, ACG) in-frame
   upstream (strand '+') or downstream (strand '-') from the peptide genomic position.
4) Compute a Kozak score around the hit (in transcript orientation).
5) Write found.tsv (hits) and not_found.tsv (misses), and coverage.tsv per accession.

Usage:
  python start_codon_scanner.py \
      --peptides peptides.xlsx \
      --genome genome.fasta \
      --originals originals.fasta \
      --outdir results \
      [--sheet SHEETNAME] \
      [--max-scan-nt 3000]

Assumptions / Notes:
- The Excel file includes at least: accessions, start, end, chrom, strand, proteins, type.
- Coordinates start/end are assumed 1-based inclusive genome coordinates for scanning.
- 'strand' is '+' or '-'.
- If a column 'peptide_seq' is available and originals.fasta contains protein sequences,
  coverage is computed by mapping peptide substrings onto the full protein sequence.
- Otherwise, if per-accession rows have start/end presumed as positions within the originals FASTA
  sequence, coverage will be computed as the union of those intervals (1-based, clamped to length).
- If neither applies, coverage is left as NaN for that accession.
"""

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}

def is_protein_sequence(seq: str) -> bool:
    """Heuristic: consider protein if >=80% of chars are not A/C/G/T/U/N."""
    s = seq.upper()
    if not s:
        return False
    nts = set("ACGTUN")
    non_nt = sum(1 for c in s if c not in nts and c.isalpha())
    alpha = sum(1 for c in s if c.isalpha())
    if alpha == 0:
        return False
    return (non_nt / alpha) >= 0.8

def read_fasta_dict(path, key_attr="id"):
    """Read a FASTA into a dict by chosen key (id or name); value is uppercase string sequence."""
    d = {}
    with open(path, "r") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            key = getattr(rec, key_attr)
            d[key] = str(rec.seq).upper()
    return d

def reverse_complement(s: str) -> str:
    return str(Seq(s).reverse_complement())

def kozak_score(context_seq: str, codon_pos: int) -> dict:
    """
    Compute Kozak score given a transcript-oriented sequence and the index codon_pos
    as the first base of the candidate start codon within context_seq.

    We expect that context_seq includes at least 6 nt upstream and 6 nt downstream,
    but we handle edges by treating missing as 'N'. Positions are relative to the AUG 'A' at index codon_pos.

    Scoring:
    - Core: -3 in {A,G} => +2 ; +4 == G => +2
    - Extended (optional simplified): matches against selected positions of GCCRCCAUGG (vertebrate-like)
      We'll check: -6:'G', -5:'C', -4:'C', -2:'C', +4:'G', +5:'G' ; each match +0.5

    Returns dict with components and total.
    """
    # Helper to fetch base with safety
    def base(rel):
        idx = codon_pos + rel
        if 0 <= idx < len(context_seq):
            return context_seq[idx]
        return 'N'

    core = 0.0
    if base(-3) in ('A','G'):
        core += 2.0
    # +4 relative to first base of codon means rel=+3 (0,1,2 are codon)
    if base(+3) == 'G':
        core += 2.0

    extended = 0.0
    consensus = { -6:'G', -5:'C', -4:'C', -2:'C', +4:'G', +5:'G' }
    for rel, ref in consensus.items():
        if base(rel) == ref:
            extended += 0.5

    return {
        "kozak_core": core,
        "kozak_ext": extended,
        "kozak_total": core + extended
    }

CODON_PRIOR = {
    "ATG": 0.0,
    "CTG": -1.0,
    "GTG": -1.2,
    "TTG": -1.3,
    "ACG": -1.7,
}

def total_score(kozak_total: float, codon: str) -> float:
    return kozak_total + CODON_PRIOR.get(codon, -2.0)

def compute_coverage(peps_df: pd.DataFrame, originals: dict) -> pd.DataFrame:
    """
    Compute per-accession coverage. Three modes:
    A) If column 'peptide_seq' exists and originals[accessions] look like proteins,
       map peptide substrings onto full sequence.
    B) Else if rows provide start/end assumed to be positions within originals sequence,
       union intervals (1-based inclusive) clamped to sequence length.
    C) Else, coverage = NaN.

    Returns dataframe columns: accession, length, covered_len, coverage
    """
    rows = []
    # Guess if originals are proteins by sampling a few sequences
    sample_seqs = [seq for _, seq in list(originals.items())[:5]]
    originals_are_protein = any(is_protein_sequence(s) for s in sample_seqs)

    by_acc = defaultdict(list)
    for _, r in peps_df.iterrows():
        acc = str(r.get("accessions"))
        by_acc[acc].append(r)

    for acc, items in by_acc.items():
        full = originals.get(acc)
        if not full:
            rows.append({"accessions": acc, "length": np.nan, "covered_len": np.nan, "coverage": np.nan})
            continue
        L = len(full)
        covered = np.zeros(L, dtype=bool)

        have_pepseq = ("peptide_seq" in peps_df.columns) and any(isinstance(x, str) and x for x in [it.get("peptide_seq") for it in items])
        if originals_are_protein and have_pepseq:
            # Mode A: peptide strings
            for it in items:
                pep = it.get("peptide_seq")
                if isinstance(pep, str) and pep:
                    pep = pep.strip().upper()
                    # Find all occurrences (allow multiple)
                    start_idx = 0
                    while True:
                        pos = full.find(pep, start_idx)
                        if pos == -1:
                            break
                        covered[pos:pos+len(pep)] = True
                        start_idx = pos + 1
        else:
            # Mode B: intervals as positions within full sequence
            # Assume 'start','end' are 1-based inclusive positions within the originals sequence.
            # If they are genomic instead, this may not be correct; we still offer a best-effort.
            for it in items:
                s = it.get("start")
                e = it.get("end")
                if pd.notnull(s) and pd.notnull(e):
                    try:
                        s = int(s)
                        e = int(e)
                    except Exception:
                        continue
                    if s > e:
                        s, e = e, s
                    s0 = max(1, s)
                    e0 = min(L, e)
                    if s0 <= e0:
                        covered[s0-1:e0] = True

        cov_len = int(covered.sum())
        coverage = float(cov_len)/float(L) if L > 0 else np.nan
        rows.append({"accessions": acc, "length": L, "covered_len": cov_len, "coverage": coverage})

    return pd.DataFrame(rows)

def get_transcript_oriented_context(genome_seq: str, codon_first_idx0: int, strand: str, flank: int = 6) -> (str, int):
    """
    Extract a context window [-flank .. +flank] around the start codon, in transcript 5'->3' orientation.
    Returns (context_seq, codon_pos_in_context).
    The input codon_first_idx0 is 0-based index on the genome_seq for the '+' strand case.
    For '-' strand, we still provide transcript-oriented context (reverse-complement).
    """
    if strand == '+':
        left = max(0, codon_first_idx0 - flank)
        right = min(len(genome_seq), codon_first_idx0 + 3 + flank)  # include the whole codon
        context = genome_seq[left:right]
        codon_pos = codon_first_idx0 - left  # index of first base of codon in context
        return context, codon_pos
    else:
        # For '-' we need to take reverse complement centered around codon (which itself on '-' is the reverse comp)
        left = max(0, codon_first_idx0 - flank)
        right = min(len(genome_seq), codon_first_idx0 + 3 + flank)
        genomic_context = genome_seq[left:right]
        rc = reverse_complement(genomic_context)
        # In reverse complement, the codon starts at position: length - (codon_end - left) - 3
        # Original codon spans [codon_first_idx0, codon_first_idx0+3) on genome '+' indexing;
        # In the genomic_context window, codon start is at (codon_first_idx0 - left),
        # codon end is at (codon_first_idx0 - left + 3)
        rel_start = codon_first_idx0 - left
        rel_end = rel_start + 3
        # In RC string, positions are reversed: index i in original maps to (len-1-i) in RC.
        # The first base of codon in transcript orientation is the RC of the third base in genomic window,
        # so position in RC = len(rc) - rel_end
        codon_pos_rc = len(rc) - rel_end
        return rc, codon_pos_rc

def scan_for_start_codon(genome_seq: str, pep_start_genomic_1based: int, pep_end_genomic_1based: int,
                         strand: str, max_scan_nt: int = 3000) -> dict:
    """
    Scan in-frame for the nearest start codon (from the five-codon set).
    For '+' strand: scan upstream (towards lower genomic coords) from peptide start.
    For '-' strand: scan downstream (towards higher genomic coords) from peptide end.
    Respect reading frame and stop when encountering an in-frame stop codon before any start codon.

    Returns a dict with keys if found:
      codon, genomic_start_1b, genomic_end_1b, distance_nt, frame_ok
    Or {} if not found.
    """
    L = len(genome_seq)
    stops = {"TAA","TAG","TGA"}

    if strand == '+':
        # Reference is peptide start (5' side)
        ref = pep_start_genomic_1based - 1  # 0-based index of first nt of peptide codon
        # Determine frame relative to ref: ref is assumed to be codon boundary
        # We'll walk upstream in steps of 3: positions [ref-3, ref), [ref-6, ref-3), ...
        steps = 0
        pos = ref
        while pos - 3 >= 0 and steps*3 <= max_scan_nt:
            codon_start = pos - 3
            triplet = genome_seq[codon_start: pos]
            if len(triplet) != 3:
                break
            if triplet in stops:
                # stop codon before any start codon in-frame
                return {}
            if triplet in START_CODONS:
                dist = (ref - codon_start)
                return {
                    "codon": triplet,
                    "genomic_start_1b": codon_start + 1,
                    "genomic_end_1b": codon_start + 3,
                    "distance_nt": dist,
                    "frame_ok": True
                }
            pos -= 3
            steps += 1
        return {}
    else:
        # '-' strand, reference is peptide end (closer to 5' in transcript)
        ref_end = pep_end_genomic_1based  # 1-based inclusive
        ref = ref_end  # we'll consider the next codon downstream
        steps = 0
        pos = ref - 1  # 0-based boundary
        while pos + 3 <= L and steps*3 <= max_scan_nt:
            codon_start = pos
            triplet = genome_seq[codon_start: codon_start+3]
            if len(triplet) != 3:
                break
            # On '-' strand, the transcript codon is reverse complement
            rc_triplet = reverse_complement(triplet)
            if rc_triplet in {"TAA","TAG","TGA"}:
                return {}
            if rc_triplet in START_CODONS:
                dist = (codon_start + 1 - (ref_end))  # positive downstream distance
                return {
                    "codon": rc_triplet,
                    "genomic_start_1b": codon_start + 1,
                    "genomic_end_1b": codon_start + 3,
                    "distance_nt": dist,
                    "frame_ok": True
                }
            pos += 3
            steps += 1
        return {}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", required=True, help="Excel file with peptide data")
    ap.add_argument("--genome", required=True, help="Genome FASTA (by chromosome names matching 'chrom' column)")
    ap.add_argument("--originals", required=True, help="Original sequences FASTA (by accessions)")
    ap.add_argument("--sheet", default=None, help="Excel sheet name (optional)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--max-scan-nt", type=int, default=3000, help="Max nt to scan upstream/downstream")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load inputs
    df = pd.read_excel(args.peptides, sheet_name=args.sheet)
    # Basic column normalization
    for col in ["accessions","chrom","strand","type"]:
        if col in df.columns:
            df[col] = df[col].astype(str)

    # Filter rows: proteins == 1 and type != "exon_diff"
    # 'proteins' might be int or string; normalize
    if "proteins" in df.columns:
        df["proteins_num"] = pd.to_numeric(df["proteins"], errors="coerce")
        mask_prot = df["proteins_num"] == 1
    else:
        print("WARNING: 'proteins' column not found; keeping all rows.", file=sys.stderr)
        mask_prot = pd.Series([True]*len(df))

    if "type" in df.columns:
        mask_type = df["type"].str.lower() != "exon_diff"
    else:
        print("WARNING: 'type' column not found; keeping all rows.", file=sys.stderr)
        mask_type = pd.Series([True]*len(df))

    df_filt = df[mask_prot & mask_type].copy()
    if df_filt.empty:
        print("No rows after filtering; exiting.", file=sys.stderr)
        # still write empty outputs
        df_filt.to_csv(os.path.join(args.outdir, "filtered_empty.tsv"), sep="\t", index=False)
        return

    # Load FASTAs
    genome = read_fasta_dict(args.genome, key_attr="id")
    originals = read_fasta_dict(args.originals, key_attr="id")

    # Compute coverage per accession (best-effort modes)
    cov_df = compute_coverage(df_filt, originals)
    cov_df.to_csv(os.path.join(args.outdir, "coverage.tsv"), sep="\t", index=False)

    # Per-row scanning for start codons
    found_rows = []
    not_found_rows = []

    required_cols = ["accessions","chrom","strand","start","end"]
    for rc in required_cols:
        if rc not in df_filt.columns:
            print(f"ERROR: required column '{rc}' missing.", file=sys.stderr)
            sys.exit(1)

    # Join coverage for quick lookup
    cov_map = {r["accessions"]: r for _, r in cov_df.iterrows()}

    for idx, row in df_filt.iterrows():
        acc = str(row["accessions"])
        chrom = str(row["chrom"])
        strand = str(row["strand"]).strip()
        start = row["start"]
        end = row["end"]

        # Prepare baseline record
        base = {
            "accessions": acc,
            "chrom": chrom,
            "strand": strand,
            "peptide_start": start,
            "peptide_end": end,
        }
        cov = cov_map.get(acc, {})
        base["coverage"] = cov.get("coverage", np.nan)

        # Fetch chromosome sequence
        chrom_seq = genome.get(chrom)
        if not chrom_seq:
            base_nf = dict(base)
            base_nf["reason"] = "missing_chrom"
            not_found_rows.append(base_nf)
            continue

        # Validate coordinates
        try:
            s = int(start)
            e = int(end)
        except Exception:
            base_nf = dict(base)
            base_nf["reason"] = "invalid_coords"
            not_found_rows.append(base_nf)
            continue

        # Clamp to chromosome length
        L = len(chrom_seq)
        if not (1 <= s <= L and 1 <= e <= L):
            # best-effort clamp
            s = max(1, min(L, s))
            e = max(1, min(L, e))

        # Decide reference depending on strand
        if strand == '+':
            hit = scan_for_start_codon(chrom_seq, s, e, strand='+', max_scan_nt=args.max_scan_nt)
        elif strand == '-':
            hit = scan_for_start_codon(chrom_seq, s, e, strand='-', max_scan_nt=args.max_scan_nt)
        else:
            base_nf = dict(base)
            base_nf["reason"] = "invalid_strand"
            not_found_rows.append(base_nf)
            continue

        if hit:
            # Transcript-oriented context and Kozak
            codon_start0 = hit["genomic_start_1b"] - 1
            ctx_seq, codon_pos = get_transcript_oriented_context(chrom_seq, codon_start0, strand, flank=6)
            kz = kozak_score(ctx_seq, codon_pos)
            sc_total = total_score(kz["kozak_total"], hit["codon"])

            out = dict(base)
            out.update({
                "codon": hit["codon"],
                "codon_genomic_start": hit["genomic_start_1b"],
                "codon_genomic_end": hit["genomic_end_1b"],
                "distance_nt": hit["distance_nt"],
                "frame_ok": hit["frame_ok"],
                "context_seq": ctx_seq,
                "kozak_core": kz["kozak_core"],
                "kozak_ext": kz["kozak_ext"],
                "kozak_total": kz["kozak_total"],
                "codon_prior": CODON_PRIOR.get(hit["codon"], -2.0),
                "score_total": sc_total,
                "is_top": 1  # one per peptide in this pipeline
            })
            found_rows.append(out)
        else:
            base_nf = dict(base)
            base_nf["reason"] = "no_inframe_start_or_stop_before"
            not_found_rows.append(base_nf)

    found_df = pd.DataFrame(found_rows)
    not_found_df = pd.DataFrame(not_found_rows)

    found_path = os.path.join(args.outdir, "found.tsv")
    not_found_path = os.path.join(args.outdir, "not_found.tsv")

    found_df.to_csv(found_path, sep="\t", index=False)
    not_found_df.to_csv(not_found_path, sep="\t", index=False)

    # Also save filtered input for traceability
    df_filt.to_csv(os.path.join(args.outdir, "filtered_input.tsv"), sep="\t", index=False)

    print(f"Wrote: {found_path}")
    print(f"Wrote: {not_found_path}")
    print(f"Wrote: {os.path.join(args.outdir, 'coverage.tsv')}")
    print(f"Wrote: {os.path.join(args.outdir, 'filtered_input.tsv')}")

if __name__ == "__main__":
    main()
