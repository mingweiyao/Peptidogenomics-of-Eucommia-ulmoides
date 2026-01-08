import os, re, math, json
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np

# ===========
# CONFIG
# ===========
CDS_FA = "/media/wanglab/caca/work_mechanism/mer/Eu_CDS.fasta"
GENOME_FA = "/media/wanglab/caca/work_mechanism/mer/Eu_genome.fasta"
GFF3_FA = "/media/wanglab/caca/work_mechanism/mer/GWHBISF00000000.gff"
OUT_DIR = "/media/wanglab/caca/work_mechanism/mer/out"

# CDS_FA = r"F:\work_mechanism\mer\Eu_CDS.fasta"
# GENOME_FA = r"F:\work_mechanism\mer\Eu_genome.fasta"
# GFF3_FA = r"F:\work_mechanism\mer\GWHBISF00000000.gff"
# OUT_DIR = r"F:\work_mechanism\mer\out"

# upstream window definition: [tis_pos-UP_START, tis_pos-UP_END)
# Here: upstream 100 nt right before the codon (exclude the codon itself)
UP_LEN_EXPECT = 100
UP_START = 100
UP_END = 0

# kmer settings
KMER_KS = [1, 2, 3, 4, 5]
TOPN_3 = 30
TOPN_5 = 30
PSEUDOCOUNT = 1.0

# TN sampling
TN_PER_TX = 5
MIN_INTERNAL_NT = 30
SEED = 13

# output
KMER_STATS_XLSX = os.path.join(OUT_DIR, "tp_tn_kmer_stats.xlsx")
TOP_MOTIFS_JSON = os.path.join(OUT_DIR, "top_motifs_3mer_5mer.json")

# =========================================================
# 1) Parse genome + gff3 + build merged seq (upstream 100nt + CDS)
# =========================================================
def parse_genome_file(genome_file):
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(genome_file, "fasta")}

def parse_gff_file(gff3_file):
    """Return set of transcript IDs that have annotated five_prime_UTR."""
    has_5utr = set()
    with open(gff3_file) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "five_prime_UTR":
                continue
            attr = parts[8]
            # faster than regex for big gff3
            key = "Parent_Accession="
            i = attr.find(key)
            if i == -1:
                continue
            val = attr[i + len(key):].split(";", 1)[0]
            if val:
                has_5utr.add(val)
    return has_5utr

def parse_tid_details(description, cds_seq, genome_dict):
    """
    Read Position=chrom:exons:strand from fasta description.
    Build merged sequence: upstream100nt (genomic) + cds_seq.
    """
    chrom, exon_part, strand = description.split("\t")[3].split("=")[1].split(": ")
    try:
        exon_coords = [(int(s), int(e)) for s, e in (b.split("-") for b in exon_part.split(","))]
    except Exception:
        return None

    g = genome_dict[chrom]
    cds_seq = cds_seq.upper()

    if strand == "+":
        start = min(s for s, _ in exon_coords)
        # need genome[start-101 : start-1] (100 nt), start is 1-based in typical annotations
        if start - 101 < 0:
            return None
        up = g[start - 101: start - 1]
        if len(up) != 100:
            return None
        return (up + cds_seq).upper()

    elif strand == "-":
        end = max(e for _, e in exon_coords)
        # need genome[end : end+100] (0-based slicing), then reverse-complement
        if end + 100 > len(g):
            return None
        up = str(Seq(g[end: end + 100]).reverse_complement())
        if len(up) != 100:
            return None
        return (up + cds_seq).upper()

    return None

def load_cds_fasta(cds_file, gff3_file, genome_dict):
    has_5utr = parse_gff_file(gff3_file)
    cds_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        tid = record.id
        if tid not in has_5utr:
            continue
        cds_seq = str(record.seq).upper()
        merged = parse_tid_details(record.description, cds_seq, genome_dict)
        if merged:
            cds_dict[tid] = merged
    return cds_dict

# =========================================================
# 2) TP / TN sampling
# TP: start ATG at position 100 (0-based) in merged seq
# TN: in-frame internal ATG in CDS part
# =========================================================
def sample_tp_tn(cds_dict, tn_per_tx=TN_PER_TX, min_internal_nt=MIN_INTERNAL_NT, seed=SEED):
    import random
    random.seed(seed)
    rows = []
    for tid, seq in cds_dict.items():
        # TP at merged index 100
        if seq[100:103] != "ATG":
            continue
        rows.append((tid, 100, seq, 1))

        cand = []
        for i in range(103, len(seq) - 2, 3):
            if i - 100 < min_internal_nt:
                continue
            if seq[i:i+3] == "ATG":
                cand.append(i)

        if cand:
            picks = random.sample(cand, k=min(tn_per_tx, len(cand)))
            for i in picks:
                rows.append((tid, i, seq, 0))

    return pd.DataFrame(rows, columns=["tid", "tis_pos", "seq", "label"])

def upstream_window(seq, tis_pos0):
    """
    Return upstream sequence [tis_pos-UP_START, tis_pos-UP_END), expected length 100.
    """
    s = tis_pos0 - UP_START
    e = tis_pos0 - UP_END
    if s < 0 or e > len(seq) or e <= s:
        return None
    up = seq[s:e]
    if len(up) != UP_LEN_EXPECT:
        return None
    return up

# =========================================================
# 3) FAST k-mer counting (numpy base-4 encoding + bincount)
# =========================================================
_BASE = np.full(256, -1, dtype=np.int16)
_BASE[ord('A')] = 0
_BASE[ord('C')] = 1
_BASE[ord('G')] = 2
_BASE[ord('T')] = 3

def _encode_seq_to_ints(seq: str):
    """Map A/C/G/T -> 0/1/2/3; return None if invalid char exists."""
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    v = _BASE[arr]
    if (v < 0).any():
        return None
    return v.astype(np.int16, copy=False)

def _kmer_codes(v: np.ndarray, k: int):
    """Compute base-4 codes for all k-mers in v (values 0..3)."""
    L = v.shape[0]
    if L < k:
        return np.empty(0, dtype=np.int32)
    n = L - k + 1
    codes = np.zeros(n, dtype=np.int32)
    for i in range(k):
        codes = (codes << 2) | v[i:i+n]
    return codes

def aggregate_kmer_counts_fast(df_tp_tn, ks=(1,2,3,4,5)):
    """
    Return:
      tp_counts[k] = np.array size 4**k with counts
      tn_counts[k] = np.array size 4**k with counts
      tp_total[k] / tn_total[k] = total number of k-mer windows summed
    """
    ks = list(ks)
    tp_counts = {k: np.zeros(4**k, dtype=np.int64) for k in ks}
    tn_counts = {k: np.zeros(4**k, dtype=np.int64) for k in ks}
    tp_total  = {k: 0 for k in ks}
    tn_total  = {k: 0 for k in ks}

    tis_arr = df_tp_tn["tis_pos"].values
    seq_arr = df_tp_tn["seq"].values
    lab_arr = df_tp_tn["label"].values

    for tis_pos, seq, lab in zip(tis_arr, seq_arr, lab_arr):
        up = upstream_window(seq, int(tis_pos))
        if up is None:
            continue
        v = _encode_seq_to_ints(up)
        if v is None:
            continue

        is_tp = int(lab) == 1
        for k in ks:
            codes = _kmer_codes(v, k)
            if codes.size == 0:
                continue
            bc = np.bincount(codes, minlength=4**k)
            if is_tp:
                tp_counts[k] += bc
                tp_total[k]  += codes.size
            else:
                tn_counts[k] += bc
                tn_total[k]  += codes.size

    return tp_counts, tn_counts, tp_total, tn_total

def codes_to_kmers(k: int):
    """Generate all kmers in ACGT order mapping to code index."""
    bases = np.array(list("ACGT"))
    out = []
    for code in range(4**k):
        x = code
        s = []
        for _ in range(k):
            s.append(bases[x & 3])
            x >>= 2
        out.append("".join(reversed(s)))
    return out

def build_kmer_stats_table_fast(tp_arr, tn_arr, tp_total, tn_total, pseudocount=PSEUDOCOUNT):
    """Build full stats table for all 4**k kmers (including zeros)."""
    V = tp_arr.size
    k = int(round(np.log(V) / np.log(4)))
    kmers = codes_to_kmers(k)

    tp_freq = (tp_arr + pseudocount) / (tp_total + pseudocount * V)
    tn_freq = (tn_arr + pseudocount) / (tn_total + pseudocount * V)
    log2_enr = np.log2(tp_freq / tn_freq)

    df = pd.DataFrame({
        "kmer": kmers,
        "tp_count": tp_arr,
        "tn_count": tn_arr,
        "tp_freq": tp_freq,
        "tn_freq": tn_freq,
        "log2_enrichment_tp_over_tn": log2_enr
    }).sort_values("log2_enrichment_tp_over_tn", ascending=False)

    return df

def pick_topN(df_stats, topN, mode="positive"):
    """
    mode:
      - "positive": 选 TP 富集最大的 topN（log2_enrichment 最大）
      - "abs":      选 |log2_enrichment| 最大的 topN（双向差异）
    """
    if mode == "abs":
        df2 = df_stats.copy()
        df2["abs_log2"] = df2["log2_enrichment_tp_over_tn"].abs()
        return df2.sort_values("abs_log2", ascending=False).head(topN)["kmer"].tolist()
    return df_stats.head(topN)["kmer"].tolist()

# =========================================================
# 4) Main: write k-mer stats (1–5) + save top 3mer/5mer
# =========================================================
def run_tp_tn_kmer_stats_and_select_topN():
    os.makedirs(OUT_DIR, exist_ok=True)

    genome_dict = parse_genome_file(GENOME_FA)
    cds_dict = load_cds_fasta(CDS_FA, GFF3_FA, genome_dict)
    df = sample_tp_tn(cds_dict)

    # Fast aggregate for all ks in one pass
    tp_counts, tn_counts, tp_total, tn_total = aggregate_kmer_counts_fast(df, ks=KMER_KS)

    top_motifs = {}
    with pd.ExcelWriter(KMER_STATS_XLSX) as w:
        for k in KMER_KS:
            stats = build_kmer_stats_table_fast(
                tp_counts[k], tn_counts[k],
                tp_total[k], tn_total[k],
                pseudocount=PSEUDOCOUNT
            )
            stats.to_excel(w, sheet_name=f"{k}mer", index=False)

            if k == 3:
                top_motifs["top_3mer"] = pick_topN(stats, TOPN_3, mode="positive")
            if k == 5:
                top_motifs["top_5mer"] = pick_topN(stats, TOPN_5, mode="positive")

        summary = pd.DataFrame({
            "metric": ["TP_samples", "TN_samples", "TOPN_3", "TOPN_5", "UP_len", "UP_START", "UP_END"],
            "value": [
                int((df["label"] == 1).sum()),
                int((df["label"] == 0).sum()),
                TOPN_3, TOPN_5, UP_LEN_EXPECT, UP_START, UP_END
            ]
        })
        summary.to_excel(w, sheet_name="summary", index=False)

    with open(TOP_MOTIFS_JSON, "w", encoding="utf-8") as f:
        json.dump({
            "top_3mer": top_motifs.get("top_3mer", []),
            "top_5mer": top_motifs.get("top_5mer", []),
            "params": {
                "UP_START": UP_START,
                "UP_END": UP_END,
                "UP_LEN_EXPECT": UP_LEN_EXPECT,
                "TN_PER_TX": TN_PER_TX,
                "MIN_INTERNAL_NT": MIN_INTERNAL_NT,
                "SEED": SEED,
                "PSEUDOCOUNT": PSEUDOCOUNT,
                "TOPN_3": TOPN_3,
                "TOPN_5": TOPN_5
            }
        }, f, indent=2)

    print(f"[OK] kmer stats written: {KMER_STATS_XLSX}")
    print(f"[OK] top motifs written: {TOP_MOTIFS_JSON}")
    print(f"[INFO] TP={int((df['label']==1).sum())}, TN={int((df['label']==0).sum())}, transcripts_kept={len(set(df['tid']))}")

if __name__ == "__main__":
    run_tp_tn_kmer_stats_and_select_topN()
