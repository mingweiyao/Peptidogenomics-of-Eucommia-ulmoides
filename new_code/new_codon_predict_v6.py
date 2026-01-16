import os, re, math, json
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from collections import Counter
from itertools import product
import numpy as np
from sklearn.linear_model import LogisticRegression

# =========================================================
# CONFIG
# =========================================================
INPUT_DIR = "F:/work_mechanism/new_file/02figure/Eu_genome_modified"
CDS_FA = os.path.join(INPUT_DIR, "Eu_CDS.fasta")
GENOME_FA = os.path.join(INPUT_DIR, "Eu_genome.fasta")
TRANSCRIPT_FA = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\hit_transcripts_extra.fa"
CANDIDATES_XLSX = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\hit_transcript_predict_orf.xlsx"
CODON_EFFICIENCY = r"F:\work_mechanism\new_file\02figure\figure4\codon\codon_efficiency.xlsx"
GFF3_FA = os.path.join(INPUT_DIR, "GWHBISF00000000.gff")
OUT_DIR = r"F:\work_mechanism\new_file\02figure\figure4\codon\codon_prediction\codon_prediction_v6"
OUT_XLSX = os.path.join(OUT_DIR, "candidates_scored.xlsx")
SHEET_NAME = "hit_transcript_predict_orf"
# upstream window: [-100, 0) relative to codon first base (tis_pos0)
UP_LEN = 100
UP_START = 100
UP_END = 0

CODON = ['ATG', 'CTG', 'ACG', 'GTG', 'TTG']
MIN_INTERNAL_NT = 30
DNA_RE_100 = re.compile(r"^[ACGT]{100}$")
DNA_RE = re.compile(r"^[ACGT]+$")
KMER = (3, 5)

MOTIFS_DNA = ["TTC", "TCT", "CTC", "TCTTC", "TCTCT"]
MOTIF_EPS = 1e-9

# =========================================================
# CODON efficiency parsing
# =========================================================
def parse_codon_file(CODON_EFFICIENCY):
    df = pd.read_excel(CODON_EFFICIENCY)
    codon_efficiency = {}
    for _, r in df.iterrows():
        codon_efficiency[str(r['TIS Sequence'])] = r['TIS Efficiency']
    return codon_efficiency
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
    if cds_seq[:3] != "ATG": return None
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
        if merged is None: continue
        cds_dict[tid] = merged
    return cds_dict
# =========================================================
# Building TP/TN
# =========================================================
def sample_tp_tn_all_internal_atg(cds_dict, min_internal_nt=MIN_INTERNAL_NT):
    rows_tp = []
    rows_tn = []
    for tid, seq in cds_dict.items():
        seq = seq.upper()
        if seq[UP_START:UP_START+3] != "ATG": continue
        rows_tp.append(('ATG', tid, UP_START, seq, 1))
        for i in range(UP_START + 3, len(seq) - 2, 3):
            if (i - UP_START) < min_internal_nt: continue
            codon = seq[i:i+3]
            if codon in CODON: rows_tn.append((codon, tid, i, seq, 0))
    return pd.DataFrame(rows_tp + rows_tn, columns=["codon", "tid", "tis_pos", "seq", "label"])
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
def build_motif(df, kmer=KMER, alpha=0.01):
    out = {}
    for k in kmer:
        obs = Counter()
        total = 0
        used_windows = 0
        for _, r in df.iterrows():
            s = r["seq"]
            up = upstream_window(s, int(r["tis_pos"]))
            if up is None: continue
            up = up.upper()
            if not DNA_RE_100.fullmatch(up): continue
            used_windows += 1
            cnt = count_kmers(up, k)
            obs.update(cnt)
            total += max(0, len(up) - k + 1)
        V = 4 ** k
        denom = total + alpha * V
        default_p = alpha / denom
        default_logp = math.log(default_p)
        logp = {km: math.log((c + alpha) / denom) for km, c in obs.items()}
        out[k] = {
            "params": {"k": k, "alpha": alpha, "total_windows": int(used_windows), "total_kmer_positions": int(total), "V": V },
            "logp": logp,
            "default_logp": float(default_logp),
        }
    return out
def train_tp_tn_lm(df):
    lm_tp = build_motif(df[df["label"]==1], KMER)
    tp_lm3, tp_lm5 = lm_tp[3], lm_tp[5]
    tn_lm3 = {}; tn_lm5 = {}
    df_tn = df[df["label"] == 0].copy()
    for codon, sub in df_tn.groupby("codon"):
        lm_bg = build_motif(sub, KMER)
        tn_lm3[codon] = lm_bg[3]
        tn_lm5[codon] = lm_bg[5]
    return tp_lm3, tp_lm5, tn_lm3, tn_lm5
# =========================================================
# Motif weights (TP vs CDS background)
# =========================================================
def lm_prob(lm_obj, kmer):
    lp = lm_obj["logp"].get(kmer, lm_obj["default_logp"])
    return math.exp(lp)
def compute_motif_weights(tp_lm3, tp_lm5, tn_lm3_dict, tn_lm5_dict, motifs=MOTIFS_DNA, eps=MOTIF_EPS):
    motif_weights_by_codon = {}
    codons = sorted(set(tn_lm3_dict.keys()) | set(tn_lm5_dict.keys()))
    for codon in codons:
        motif_weights_by_codon[codon] = {}
        for m in motifs:
            k = len(m)
            if k == 3:
                p_tp = lm_prob(tp_lm3, m)
                if codon not in tn_lm3_dict:
                    p_bg = float("nan")
                    log_ratio = float("nan")
                else:
                    p_bg = lm_prob(tn_lm3_dict[codon], m)
                    log_ratio = math.log((p_tp + eps) / (p_bg + eps))
            elif k == 5:
                p_tp = lm_prob(tp_lm5, m)
                if codon not in tn_lm5_dict:
                    p_bg = float("nan")
                    log_ratio = float("nan")
                else:
                    p_bg = lm_prob(tn_lm5_dict[codon], m)
                    log_ratio = math.log((p_tp + eps) / (p_bg + eps))
            else:
                raise ValueError(f"Unsupported motif length: {m}")
            motif_weights_by_codon[codon][m] = {
                "p_tp": float(p_tp),
                "p_bg": float(p_bg),
                "log_ratio": float(log_ratio),
            }
    return motif_weights_by_codon
# =========================================================
# Training (PWM + LR)
# =========================================================
def cu_fraction(up):
    return (sum(1 for ch in up if ch in "CT") / len(up)) if up else float("nan")
def score_seq_ll(seq, lm_obj, k):
    if seq is None or len(seq) < k: return float("nan")
    seq = seq.upper()
    L = len(seq)
    denom = (L - k + 1)
    if denom <= 0: return float("nan")
    logp = lm_obj["logp"]
    dlog = lm_obj["default_logp"]
    s = 0.0
    for i in range(denom):
        km = seq[i:i+k]
        if all(ch in "ACGT" for ch in km): s += logp.get(km, dlog)
        else: return float("nan")
    return float(s / denom)
def motif_score(up, motif_weights):
    if up is None: return float("nan")
    up = up.upper()
    if not DNA_RE.fullmatch(up): return float("nan")
    s = 0.0
    for m, d in motif_weights.items():
        w = float(d["log_ratio"])
        s += motif_rate(up, m) * w
    return float(s)
def motif_rate(seq, motif):
    L = len(seq)
    k = len(motif)
    if L < k: return 0.0
    denom = L - k + 1
    cnt = 0
    for i in range(denom):
        if seq[i:i+k] == motif:
            cnt += 1
    return cnt / denom if denom > 0 else 0.0
def train_pwm_and_lr(df, motif_weights, codon):
    X = []; y = []
    for _, r in df[(df['label']==1) | ((df['label']==0 )&(df['codon']==codon))].iterrows():
        seq = r["seq"].upper()
        tis = int(r["tis_pos"])
        up = upstream_window(seq, tis)
        if not DNA_RE_100.fullmatch(up): continue
        cu = cu_fraction(up)
        ms = motif_score(up, motif_weights[codon])
        if any(np.isnan(v) for v in [cu, ms]): continue
        X.append([cu, ms])
        y.append(int(r["label"]))
    X = np.asarray(X, float)
    y = np.asarray(y, int)
    if X.shape[0] < 50: raise ValueError(f"训练样本过少：{X.shape[0]}")
    # standardize
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    lr = LogisticRegression(penalty="l2", solver="liblinear", max_iter=2000, class_weight="balanced")
    lr.fit(Z, y)
    weights = {
        "beta0": float(lr.intercept_[0]),
        "w1": float(lr.coef_[0][0]),  # cu_fraction
        "w2": float(lr.coef_[0][1]),  # motif_score
    }
    std = {
        "mu": mu.tolist(),
        "sd": sd.tolist(),
        "feat": ["cu_fraction", "motif_score"],
    }
    return weights, std
# =========================================================
# Candidate scoring helpers
# =========================================================
def make_candidate_key(row):
    return f'{row["chrom"]}|{row["strand"]}|{int(row["phy_start"])}|{int(row["phy_end"])}|{str(row["codon"]).upper()}'
def write_lm_tsv(lm_obj, out_tsv):
    rows = []
    for kmer, logp in lm_obj["logp"].items():
        rows.append({
            "kmer": kmer,
            "logp": logp
        })
    df = pd.DataFrame(rows)
    df = df.sort_values("kmer")
    df.to_csv(out_tsv, sep="\t", index=False)

def score_candidates_excel(require_5utr=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    # ---- 0) Load genome, codon_efficiency & CDS ----
    genome_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(GENOME_FA, "fasta")}
    transcript_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(TRANSCRIPT_FA, "fasta")}
    codon_efficiency = parse_codon_file(CODON_EFFICIENCY)
    cds_dict = load_cds_fasta(CDS_FA, GFF3_FA, genome_dict, require_5utr=require_5utr)
    # ---- 1) Train TP LM + BG LM (3mer/5mer) ----
    df = sample_tp_tn_all_internal_atg(cds_dict)
    tp_lm3, tp_lm5, tn_lm3_dict, tn_lm5_dict = train_tp_tn_lm(df)
    # ---- 2) Compute motif weights (TP vs CDS background) ----
    motif_weights = compute_motif_weights(tp_lm3, tp_lm5, tn_lm3_dict, tn_lm5_dict, motifs=MOTIFS_DNA)
    # ---- 3) Train PWM + LR weights ----
    model = {}
    for codon in CODON:
        weights, std = train_pwm_and_lr(df, motif_weights, codon)
        model[codon] = {"weights": weights, "std": std}
    # ---- 4) Read candidates ----
    df_sp = pd.read_excel(CANDIDATES_XLSX, sheet_name=SHEET_NAME)
    candidate_keys = []
    cu_fracs = []
    motif_scores = []
    upstream_score = []
    codon_score = []
    for _, r in df_sp.iterrows():
        tran_id = r['trans_id']
        strand = str(r.get("strand", ""))
        start = int(r.get("phy_start", -1))
        end = int(r.get("phy_end", -1))
        trans_start = start + 100
        trans_end = end + 100
        codon = r['codon']
        candidate_keys.append(make_candidate_key(r))
        mu = np.array(model[codon]['std']["mu"], float)
        sd = np.array(model[codon]['std']["sd"], float)
        g = transcript_dict[tran_id]
        if strand == "+":
            left = trans_start - (UP_START + 1)
            right = trans_start - 1
            up = g[left:right] if (left >= 0 and right <= len(g)) else None
        else:
            left = trans_end
            right = trans_end + UP_START
            up = str(Seq(g[left:right]).reverse_complement()) if (left >= 0 and right <= len(g)) else None
        if up is None or len(up) != UP_LEN or (not DNA_RE_100.fullmatch(up)):
            cu_fracs.append(np.nan); motif_scores.append(np.nan)
            continue
        cu = cu_fraction(up)
        ms = motif_score(up, motif_weights[codon])
        z = (np.array([cu, ms], float) - mu) / sd
        score_base = (
            model[codon]['weights']["beta0"]
            + model[codon]['weights']["w1"] * z[0]
            + model[codon]['weights']["w2"] * z[1]
        )
        cu_fracs.append(cu)
        motif_scores.append(ms)
        upstream_score.append(score_base) 
        # codon and background
        w0 = str(r.get("kozak_seq", "")).upper()[2:10]
        if DNA_RE.fullmatch(w0):
            eff = codon_efficiency.get(w0)
            if eff is not None and not pd.isna(eff):
                codon_score.append(codon_efficiency[w0])
            else:
                codon_score.append(np.nan)
        else:
            codon_score.append(np.nan)
    df_sp["cu_fraction"] = cu_fracs
    df_sp["motif_score"] = motif_scores
    df_sp['upstream_score'] = upstream_score
    df_sp['codon_efficiency'] = codon_score
    df_sp["codon_rank"] = (
        df_sp
        .groupby("trans_id")["codon_efficiency"]
        .rank(method="min", ascending=False)
        .astype("Int64")
    )
    df_sp["upstream_rank"] = (
        df_sp
        .groupby("trans_id")["upstream_score"]
        .rank(method="min", ascending=False)
        .astype("Int64")
    )
    df_sp.to_excel(OUT_XLSX, index=False)
    # ---- 5) save model artifacts ----
    # 5.1 保存 LR 模型（每种 codon 的 weights + standardizer）
    with open(os.path.join(OUT_DIR, "lr_model_by_codon.json"), "w", encoding="utf-8") as f:
        json.dump(model, f, ensure_ascii=False, indent=2)
    # 5.2 保存 TP 的 3mer/5mer LM
    with open(os.path.join(OUT_DIR, "lm_tp_3.json"), "w", encoding="utf-8") as f:
        json.dump(tp_lm3, f, ensure_ascii=False, indent=2)
    with open(os.path.join(OUT_DIR, "lm_tp_5.json"), "w", encoding="utf-8") as f:
        json.dump(tp_lm5, f, ensure_ascii=False, indent=2)
    # 5.3 保存 TN 背景 LM（按 codon 分开）
    with open(os.path.join(OUT_DIR, "lm_bg_3_by_codon.json"), "w", encoding="utf-8") as f:
        json.dump(tn_lm3_dict, f, ensure_ascii=False, indent=2)
    with open(os.path.join(OUT_DIR, "lm_bg_5_by_codon.json"), "w", encoding="utf-8") as f:
        json.dump(tn_lm5_dict, f, ensure_ascii=False, indent=2)
    # 5.4 保存 motif 权重（按 codon）
    with open(os.path.join(OUT_DIR, "motif_weights_by_codon.json"), "w", encoding="utf-8") as f:
        json.dump(motif_weights, f, ensure_ascii=False, indent=2)
    # ---- 6) save tp/tn.tsv ----
    # 6.1 TP LM tsv
    write_lm_tsv(tp_lm3, os.path.join(OUT_DIR, "tp_lm3.tsv"))
    write_lm_tsv(tp_lm5, os.path.join(OUT_DIR, "tp_lm5.tsv"))
    # 6.2 TN LM tsv（按 codon 分文件）
    for codon, lm3 in tn_lm3_dict.items():
        write_lm_tsv(lm3, os.path.join(OUT_DIR, f"bg_{codon}_lm3.tsv"))
    for codon, lm5 in tn_lm5_dict.items():
        write_lm_tsv(lm5, os.path.join(OUT_DIR, f"bg_{codon}_lm5.tsv"))
    print("[OK] wrote:", OUT_XLSX)
    print("[OK] artifacts in:", OUT_DIR)
if __name__ == "__main__":
    score_candidates_excel(require_5utr=False)