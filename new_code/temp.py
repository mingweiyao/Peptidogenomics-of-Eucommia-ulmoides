import os, re, random, math
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from collections import Counter

# =========================================================
# CONFIG
# =========================================================
INPUT_DIR = "/media/wanglab/caca/work_mechanism/mer"
CDS_FA = os.path.join(INPUT_DIR, "Eu_CDS.fasta")
GENOME_FA = os.path.join(INPUT_DIR, "Eu_genome.fasta")
CANDIDATES_XLSX = os.path.join(INPUT_DIR, "output_candidates.xlsx")
GFF3_FA = os.path.join(INPUT_DIR, "GWHBISF00000000.gff")
OUT_DIR = os.path.join(INPUT_DIR, "codon_prediction_v2")
OUT_XLSX = os.path.join(OUT_DIR, "candidates_scored.xlsx")
SHEET_NAME = "output_candidates"

# upstream window: [-100, 0) relative to codon first base (tis_pos0)
UP_LEN = 100
UP_START = 100
UP_END = 0

# TN sampling (只用于 LR 学权重；你仍可用 multi-seed 来弱化抽样偶然性)
TN_PER_TX = 5
MIN_INTERNAL_NT = 30
TN_SEEDS = range(1, 11)
DEDUP_TN = True
KMER = (3, 5)

# Motifs of interest (RNA -> DNA)
# UUC UCU UCC UCUUC UCUCU -> TTC TCT TCC TCTTC TCTCT
MOTIFS_DNA = ["TTC", "TCT", "TCC", "TCTTC", "TCTCT"]
MOTIF_EPS = 1e-9
CODON_BONUS_SCHEMES = {
    "weak":   {"ATG": 0.0, "CTG": -0.1, "GTG": -0.2, "TTG": -0.2, "ACG": -0.3},
    "medium": {"ATG": 0.0, "CTG": -0.5, "GTG": -1.0, "TTG": -1.0, "ACG": -2.0},
    "strong": {"ATG": 0.0, "CTG": -1.0, "GTG": -2.0, "TTG": -2.0, "ACG": -3.0},
}
DEFAULT_OTHER_CODON_BONUS = -3.0

# =========================================================
# Genome background (A/C/G/T)
# =========================================================
def compute_genome_bg(genome_dict):
    c = Counter()
    for chrom, seq in genome_dict.items():
        seq = seq.upper()
        for ch in seq:
            if ch in "ACGT":
                c[ch] += 1
    tot = sum(c[b] for b in "ACGT")
    return {b: (c[b] / tot) for b in "ACGT"}
# =========================================================
# CDS parsing
# =========================================================
def parse_gff_file(gff3_file):
    has_5utr = set()
    with open(gff3_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            if parts[2] != "five_prime_UTR": continue
            m = re.search(r"Parent_Accession=([^;]+)", parts[8])
            if m: has_5utr.add(m.group(1))
    return has_5utr
def parse_tid_details(description, cds_seq, genome_dict):
    chrom, exon_part, strand = description.split("\t")[3].split("=")[1].split(": ")
    exon_coords = [(int(s), int(e)) for s, e in (b.split("-") for b in exon_part.split(","))]
    if chrom not in genome_dict: return None
    g = genome_dict[chrom]
    cds_seq = str(cds_seq).upper()
    if strand == "+":
        start = min(s for s, _ in exon_coords)
        if start - (UP_START + 1) < 0: return None
        up_seq = g[start - (UP_START + 1): start - 1] + cds_seq
    elif strand == "-":
        end = max(e for _, e in exon_coords)
        if end + UP_START > len(g): return None
        up_seq = str(Seq(g[end: end + UP_START]).reverse_complement()) + cds_seq
    else: return None
    return up_seq.upper()
def load_cds_fasta(cds_file, gff3_file, genome_dict, require_5utr=False):
    has_5utr = parse_gff_file(gff3_file) if require_5utr else None
    cds_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        tid = record.id
        if require_5utr and (tid not in has_5utr): continue
        merged = parse_tid_details(record.description, str(record.seq), genome_dict)
        if merged is None:
            continue
        cds_dict[tid] = merged
    return cds_dict
# =========================================================
# k-mer LM
# =========================================================
def count_kmers(seq, k):
    c = Counter()
    L = len(seq)
    for i in range(L - k + 1):
        km = seq[i:i+k]
        if all(ch in "ACGT" for ch in km): c[km] += 1
    return c
def upstream_window(seq, tis_pos0):
    s = tis_pos0 - UP_START
    e = tis_pos0 - UP_END
    if s < 0 or e > len(seq) or e <= s: return None
    up = seq[s:e].upper()
    return up if len(up) == UP_LEN else None
def build_motif(tp_df, kmer = KMER, alpha=1.0):
    out = {}
    for k in kmer:
        obs = Counter()
        total = 0
        for _, r in tp_df.iterrows():
            s = r["seq"]
            up = upstream_window(s, int(r["tis_pos"]))
            if up is None: continue
            cnt = count_kmers(up, k)
            obs.update(cnt)
            total += max(0, len(up) - k + 1)
        V = 4 ** k
        denom = total + alpha * V
        default_p = alpha / denom
        default_logp = math.log(default_p)
        logp = {km: math.log((c + alpha) / denom) for km, c in obs.items()}
        out[k] = {
            "params": {"k": k, "alpha": alpha, "total_windows": int(tp_df.shape[0]), "total_kmer_positions": int(total), "V": V},
            "logp": logp,
            "default_logp": float(default_logp),
        }
    return out
def train_tp_tn_lm(df):
    lm_tp = build_motif(df[df["label"]==1], KMER)
    tp_lm3, tp_lm5 = lm_tp[3], lm_tp[5]
    lm_tn = build_motif(df[df["label"]==0], KMER)
    tn_lm3, tn_lm5 = lm_tn[3], lm_tn[5]
    return tp_lm3, tp_lm5, tn_lm3, tn_lm5
def sample_tp_tn_multi_seed(cds_dict, tn_per_tx=TN_PER_TX, min_internal_nt=MIN_INTERNAL_NT, seeds=TN_SEEDS, dedup_tn=DEDUP_TN):
    rows_tp = []
    for tid, seq in cds_dict.items():
        seq = seq.upper()
        if seq[UP_START:UP_START+3] != "ATG": continue
        rows_tp.append((tid, UP_START, seq, 1))
    tn_set = set()
    rows_tn = []
    for seed in seeds:
        rng = random.Random(seed)
        for tid, seq in cds_dict.items():
            seq = seq.upper()
            if seq[UP_START:UP_START+3] != "ATG": continue
            cand = []
            for i in range(UP_START + 3, len(seq) - 2, 3):
                if (i - UP_START) < min_internal_nt: continue
                if seq[i:i+3] == "ATG": cand.append(i)
            if not cand: continue
            picks = rng.sample(cand, k=min(tn_per_tx, len(cand)))
            for i in picks:
                if dedup_tn:
                    key = (tid, i)
                    if key in tn_set: continue
                    tn_set.add(key)
                rows_tn.append((tid, i, seq, 0))
    return pd.DataFrame(rows_tp + rows_tn, columns=["tid", "tis_pos", "seq", "label"])
# =========================================================
# Main scoring pipeline
# =========================================================
def score_candidates_excel(require_5utr=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    # ---- 0) Load genome & CDS ----
    genome_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(GENOME_FA, "fasta")}
    genome_bg = compute_genome_bg(genome_dict)
    cds_dict = load_cds_fasta(CDS_FA, GFF3_FA, genome_dict, require_5utr=require_5utr)
    # ---- 1) Train TP LM + BG LM (3mer/5mer) ----
    df = sample_tp_tn_multi_seed(cds_dict)
    lm_tp3, lm_tp5, lm_tn3, lm_tn5 = train_tp_tn_lm(df)


if __name__ == "__main__":
    score_candidates_excel(require_5utr=False)
