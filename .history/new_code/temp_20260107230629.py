import os, re, math
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression

# ===========
# CONFIG
# ===========
CDS_FA = "data/CDS.fa"
GENOME_FA = "data/genome.fa"
CANDIDATES_XLSX = "data/candidates.xlsx"
GFF3_FA = "data/genome.gff3"
OUT_XLSX = "out/candidates_scored.xlsx"

MOTIFS = ["TCTTC", "TCTCT"]  # UCUUC / UCUCU
MOTIF_MODE = "binary" 

def parse_genome_file(genome_file):
    seq_dict = {}
    for rec in SeqIO.parse(genome_file, "fasta"):
        seq_dict[rec.id] = str(rec.seq)
    return seq_dict
def parse_gff3_file(gff3_path):
    has_5utr = set()
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attributes = parts[8]
            if feature_type == "five_prime_UTR":
                match = re.search(r'Parent_Accession=([^;]+)', attributes)
                tid = match.group(1)
                has_5utr.add(tid)
    return has_5utr
def parse_tid_details(description, genome_dict):
    m = re.search(r"Position=([^ ]+)", description)
    pos_str = m.group(1)
    chrom, exon_part, strand = pos_str.split(":", 2)
    exon_coords = []
    for block in exon_part.split(","):
        s, e = block.split("-")
        exon_coords.append((int(s), (int(e))))
    if strand == "+":
        A_genome = min(s for s, _ in exon_coords)
        start, end = A_genome, A_genome + 2
        kozak_seq = genome_dict[chrom][start-7:end+4]
        up = genome_dict[chrom][start-121:start-20]
    else:
        A_genome = max(e for _, e in exon_coords)
        start, end = A_genome - 2, A_genome
        kozak_seq = Seq(genome_dict[chrom][start-5:end+6]).reverse_complement()
        up = Seq(genome_dict[chrom][end+19:end+120]).reverse_complement()
    return kozak_seq, up
def load_cds_fasta(cds_file, gff3_file, genome):
    has_5utr = parse_gff3_file(gff3_file)
    d = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        tid = record.id
        if tid not in has_5utr:
            continue
        kozak_seq, up_seq = parse_tid_details(record.description, genome)
        d[tid] = {"kozak_seq": kozak_seq, "up_seq": up_seq, "cds_seq": str(record.seq)}
    return d

def compute_bg(cds_dict):
    c = {b:0 for b in "ACGT"}; tot=0
    for rec in cds_dict.values():
        for ch in rec["cds_seq"]:
            if ch in c:
                c[ch] += 1
                tot += 1
    return {b:(c[b]/tot if tot else 0.25) for b in "ACGT"}  
def sample_tp_tn(cds_dict, tn_per_tx=5, min_internal_nt=30, seed=13):
    import random
    random.seed(seed)
    rows = []
    for tid, values in cds_dict.items():
        s = values['cds_seq']
        if s[:3] == "ATG":
            rows.append((tid, 0 ,s, 1))
        else: continue
        cand = []
        for i in range(3, len(s)-2, 3):
            if i <min_internal_nt: continue
            if s[i:i+3] == "ATG":
                cand.append(i)
        if cand:
            picks = random.sample(cand, k=min(tn_per_tx, len(cand)))
            for i in picks:
                rows.append((tid, i, s, 0))
    return pd.DataFrame(rows, columns=["tid", "tis_pos", "seq", "label"])
def build_pwm(tp_windows, pseudocount=1.0):
    L = len(tp_windows[0])
    counts = {b: np.full(L, pseudocount) for b in "ACGT"}
    for w in tp_windows:
        for i, ch in enumerate(w):
            if ch in counts: counts[ch][i] += 1
    totals = sum(counts[b] for b in "AGCT")
    return {b: (counts[b]/totals).tolist() for b in "AGCT"}
def kozak_logodds(w, pwm, bg, eps=1e-12):
    s = 0.0
    for i, ch in enumerate(w):
        if ch in "ACGT":
            ratio = pwm[ch][i] / max(bg[ch], eps)
            s += math.log2(max(ratio, eps))
    return s
def cu_fraction(up):
    return (sum(1 for ch in up if ch in "CT")/len(up)) if up else float("nan")
def cu_motif_score(up):
    if MOTIF_MODE=="binary":
        return float(sum(1 for m in MOTIFS if m in up))
    tot=0
    for m in MOTIFS:
        tot += len(re.findall(f"(?={m})", up))
    return float(tot)
def train_pwm_and_weights(cds_dict, genome_dict):
    bg = compute_bg(cds_dict)
    df = sample_tp_tn(cds_dict)
    tp_ws = []
    for _, r in df.iterrows():
        seq = r["seq"]; tis=int(r["tis_pos"])
        w = cds_dict[r['tid']]["kozak_seq"]
        tp_ws.append(w)
    pwm = build_pwm(tp_ws)
    # features for LR
    X = []; y = []
    for _, r in df.iterrows():
        seq = r['seq']; tis = int(r['tis_pos'])
        w = cds_dict[r['tid']]["kozak_seq"]
        up = cds_dict[r['tid']]["up_seq"]
        X.append([kozak_logodds(w,pwm,bg), cu_fraction(up), cu_motif_score(up)])
        y.append(int(r["label"]))
    X=np.asarray(X,float); y=np.asarray(y,int)
    mu=X.mean(axis=0); sd=X.std(axis=0); sd[sd==0]=1.0
    Z=(X-mu)/sd
    lr=LogisticRegression(penalty="l2", solver="liblinear", max_iter=2000)
    lr.fit(Z,y)
    weights={
        "w1": float(lr.coef_[0][0]),
        "w2": float(lr.coef_[0][1]),
        "w3": float(lr.coef_[0][2]),
        "beta0": float(lr.intercept_[0]),
    }
    std={"mu": mu.tolist(), "sd": sd.tolist(), "feat":["kozak","cu_fraction","cu_motif"]}
    return pwm, bg, weights, std

def score_candidates_excel():
    os.makedirs(os.path.dirname(OUT_XLSX), exist_ok=True)
    genome = parse_genome_file(GENOME_FA)
    cds = load_cds_fasta(CDS_FA, GFF3_FA, genome)
    pwm, bg, weights, std = train_pwm_and_weights(cds, genome)
if __name__ == "__main__":
    score_candidates_excel()