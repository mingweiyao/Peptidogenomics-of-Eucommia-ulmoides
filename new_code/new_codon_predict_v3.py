import os, re, random, math,json
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from collections import Counter
import numpy as np
from sklearn.linear_model import LogisticRegression

# =========================================================
# CONFIG
# =========================================================
INPUT_DIR = "/media/wanglab/caca/work_mechanism/mer"
CDS_FA = os.path.join(INPUT_DIR, "Eu_CDS.fasta")
GENOME_FA = os.path.join(INPUT_DIR, "Eu_genome.fasta")
CANDIDATES_XLSX = os.path.join(INPUT_DIR, "output_candidates.xlsx")
GFF3_FA = os.path.join(INPUT_DIR, "GWHBISF00000000.gff")
OUT_DIR = os.path.join(INPUT_DIR, "codon_prediction_v3")
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

DNA_RE_100 = re.compile(r"^[ACGT]{100}$")
DNA_RE = re.compile(r"^[ACGT]+$")

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
    lm_tn = build_motif(df[df["label"]==0], KMER)
    tn_lm3, tn_lm5 = lm_tn[3], lm_tn[5]
    return tp_lm3, tp_lm5, tn_lm3, tn_lm5
# =========================================================
# Building TP/TN
# =========================================================
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
# Motif weights (TP vs CDS background)
# =========================================================
def lm_prob(lm_obj, kmer):
    lp = lm_obj["logp"].get(kmer, lm_obj["default_logp"])
    return math.exp(lp)
def compute_motif_weights(lm_tp3, lm_tp5, lm_tn3, lm_tn5, motifs=MOTIFS_DNA, eps=MOTIF_EPS):
    weights = {}
    for m in motifs:
        k = len(m)
        if k==3:
            p_tp = lm_prob(lm_tp3, m)
            p_tn = lm_prob(lm_tn3, m)
        elif k == 5:
            p_tp = lm_prob(lm_tp5, m)
            p_tn = lm_prob(lm_tn5, m)
        else:
            raise ValueError(f"Unsupported motif length: {m}")
        w = math.log((p_tp + eps) / (p_tn + eps))
        weights[m] = {
            "p_tp": float(p_tp),
            "p_bg": float(p_tn),
            "log_ratio": float(w),
        }
    return weights  
# =========================================================
# Training (PWM + LR)
# =========================================================
def kozak_window(seq, tis_pos):
    s = tis_pos - 6
    e = tis_pos + 7
    if s < 0 or e > len(seq): return None
    w = seq[s:e].upper()
    return w if len(w) == 13 else None
def build_pwm(tp_windows, pseudocount=1.0):
    if not tp_windows: raise ValueError("TP kozak windows 为空，无法建 PWM")
    counts = {b: np.full(13, pseudocount, dtype=float) for b in "ACGT"}
    for w in tp_windows:
        if len(w) != 13: continue
        for i, ch in enumerate(w):
            if ch in counts: counts[ch][i] += 1.0
    totals = sum(counts[b] for b in "ACGT")
    return {b: (counts[b] / totals).tolist() for b in "ACGT"}
def kozak_logodds(w, pwm, bg, eps=1e-12):
    s = 0.0
    for i, ch in enumerate(w):
        if ch in "ACGT":
            p = max(pwm[ch][i], eps)
            q = max(bg.get(ch, 0.25), eps)
            s += math.log(p / q)
    return float(s)
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
def train_pwm_and_lr(df, genome_bg, motif_weights):
    # 1) PWM from TP windows
    tp_ws = []
    for _, r in df[df["label"] == 1].iterrows():
        w = kozak_window(r["seq"], int(r["tis_pos"]))
        if w is not None and DNA_RE.fullmatch(w): tp_ws.append(w)
    pwm = build_pwm(tp_ws)
    # 2) Features
    X = []; y = []
    for _, r in df.iterrows():
        seq = r["seq"].upper()
        tis = int(r["tis_pos"])
        w = kozak_window(seq, tis)
        up = upstream_window(seq, tis)
        if w is None or up is None: continue
        if not DNA_RE.fullmatch(w): continue
        if not DNA_RE_100.fullmatch(up): continue
        kz = kozak_logodds(w, pwm, genome_bg)
        cu = cu_fraction(up)
        # ll_tp3 = score_seq_ll(up, lm_tp3, 3)
        # ll_tp5 = score_seq_ll(up, lm_tp5, 5)
        # ll_bg3 = score_seq_ll(up, lm_tn3, 3)
        # ll_bg5 = score_seq_ll(up, lm_tn5, 5)
        # llr3 = ll_tp3 - ll_bg3
        # llr5 = ll_tp5 - ll_bg5
        ms = motif_score(up, motif_weights)
        if any(np.isnan(v) for v in [kz, cu, ms]): continue
        X.append([kz, cu, ms])
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
        "w1": float(lr.coef_[0][0]),  # kozak_logodds
        "w2": float(lr.coef_[0][1]),  # cu_fraction
        "w3": float(lr.coef_[0][2]),  # motif_score
    }
    std = {
        "mu": mu.tolist(),
        "sd": sd.tolist(),
        "feat": ["kozak_logodds", "cu_fraction", "motif_score"],
    }
    return pwm, weights, std
# =========================================================
# Candidate scoring helpers
# =========================================================
def codon_bonus_by_scheme(codon, scheme_dict):
    c = str(codon).upper()
    return float(scheme_dict.get(c, DEFAULT_OTHER_CODON_BONUS))
def make_candidate_key(row):
    return f'{row["chrom"]}|{row["strand"]}|{int(row["phy_start"])}|{int(row["phy_end"])}|{str(row["codon"]).upper()}'
# =========================================================
# Main scoring pipeline
# =========================================================
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
    # ---- 0) Load genome & CDS ----
    genome_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(GENOME_FA, "fasta")}
    genome_bg = compute_genome_bg(genome_dict)
    cds_dict = load_cds_fasta(CDS_FA, GFF3_FA, genome_dict, require_5utr=require_5utr)
    # ---- 1) Train TP LM + BG LM (3mer/5mer) ----
    df = sample_tp_tn_multi_seed(cds_dict)
    lm_tp3, lm_tp5, lm_tn3, lm_tn5 = train_tp_tn_lm(df)
    # ---- 2) Compute motif weights (TP vs CDS background) ----
    motif_weights = compute_motif_weights(lm_tp3, lm_tp5, lm_tn3, lm_tn5, motifs=MOTIFS_DNA)
    # ---- 3) Train PWM + LR weights ----
    pwm, weights, std = train_pwm_and_lr(df, genome_bg, motif_weights)
    # ---- 4) Read candidates ----
    df_sp = pd.read_excel(CANDIDATES_XLSX, sheet_name=SHEET_NAME)
    mu = np.array(std["mu"], float)
    sd = np.array(std["sd"], float)
    candidate_keys = []
    kozak_scores = []
    cu_fracs = []
    motif_scores = []
    tis_scores = {name: [] for name in CODON_BONUS_SCHEMES}  
    for _, r in df_sp.iterrows():
        chrom = str(r.get("chrom", ""))
        strand = str(r.get("strand", ""))
        start = int(r.get("phy_start", -1))
        end = int(r.get("phy_end", -1))
        codon = str(r.get("codon", "")).upper()
        candidate_keys.append(make_candidate_key(r))
        g = genome_dict[chrom]
        w = str(r.get("kozak", "")).upper()[:13]
        if len(w) != 13 or (not DNA_RE.fullmatch(w)):
            kozak_scores.append(np.nan); cu_fracs.append(np.nan)
            motif_scores.append(np.nan)
            for name in CODON_BONUS_SCHEMES:
                tis_scores[name].append(np.nan)
            continue
        if strand == "+":
            left = start - (UP_START + 1)
            right = start - 1
            up = g[left:right] if (left >= 0 and right <= len(g)) else None
        else:
            left = end
            right = end + UP_START
            up = str(Seq(g[left:right]).reverse_complement()) if (left >= 0 and right <= len(g)) else None
        if up is None or len(up) != UP_LEN or (not DNA_RE_100.fullmatch(up)):
            kozak_scores.append(np.nan); cu_fracs.append(np.nan)
            motif_scores.append(np.nan)
            for name in CODON_BONUS_SCHEMES:
                tis_scores[name].append(np.nan)
            continue
        kz = kozak_logodds(w, pwm, genome_bg)
        cu = cu_fraction(up)
        # ll_tp3 = score_seq_ll(up, lm_tp3, 3)
        # ll_tp5 = score_seq_ll(up, lm_tp5, 5)
        # ll_bg3 = score_seq_ll(up, lm_tn3, 3)
        # ll_bg5 = score_seq_ll(up, lm_tn5, 5)
        # llr3 = ll_tp3 - ll_bg3
        # llr5 = ll_tp5 - ll_bg5
        ms = motif_score(up, motif_weights)
        z = (np.array([kz, cu, ms], float) - mu) / sd
        score_base = (
            weights["beta0"]
            + weights["w1"] * z[0]
            + weights["w2"] * z[1]
            + weights["w3"] * z[2]
        )
        kozak_scores.append(kz)
        cu_fracs.append(cu)
        motif_scores.append(ms)
        for scheme_name, scheme_dict in CODON_BONUS_SCHEMES.items():
            tis_scores[scheme_name].append(score_base + codon_bonus_by_scheme(codon, scheme_dict))
    # ---- 5) write outputs ----
    df_sp["candidate_key"] = candidate_keys
    df_sp["kozak_score"] = kozak_scores
    df_sp["cu_fraction"] = cu_fracs
    df_sp["motif_score"] = motif_scores
    for scheme_name in CODON_BONUS_SCHEMES:
        score_col = f"tis_score_{scheme_name}"
        rank_col = f"rank_{scheme_name}"
        df_sp[score_col] = tis_scores[scheme_name]
        df_sp[rank_col] = df_sp.groupby("accession")[score_col].rank(ascending=False, method="first")
        top3 = df_sp[df_sp[rank_col] <= 3].sort_values(["accession", rank_col])
        top3.to_excel(OUT_XLSX.replace(".xlsx", f"_top3_{scheme_name}.xlsx"), index=False)
    df_sp.to_excel(OUT_XLSX, index=False)
    # sensitivity summary
    top1 = {}
    for scheme_name in CODON_BONUS_SCHEMES:
        score_col = f"tis_score_{scheme_name}"
        idx = df_sp.groupby("accession")[score_col].idxmax()
        idx = idx.dropna().astype(int)
        tmp = df_sp.loc[idx, ["accession", "candidate_key"]].set_index("accession")["candidate_key"]
        top1[scheme_name] = tmp
    summary = pd.DataFrame({
        "top1_weak": top1.get("weak"),
        "top1_medium": top1.get("medium"),
        "top1_strong": top1.get("strong"),
    })
    summary["stable_top1"] = (summary["top1_weak"] == summary["top1_medium"]) & (summary["top1_medium"] == summary["top1_strong"])
    summary = summary.reset_index()
    summary.to_excel(OUT_XLSX.replace(".xlsx", "_sensitivity_summary.xlsx"), index=False)
    # ---- 8) save model artifacts ----
    with open(os.path.join(OUT_DIR, "genome_bg.json"), "w", encoding="utf-8") as f:
        json.dump(genome_bg, f, indent=2)
    with open(os.path.join(OUT_DIR, "pwm.json"), "w", encoding="utf-8") as f:
        json.dump(pwm, f, indent=2)
    with open(os.path.join(OUT_DIR, "weights.json"), "w", encoding="utf-8") as f:
        json.dump(weights, f, indent=2)
    with open(os.path.join(OUT_DIR, "standardizer.json"), "w", encoding="utf-8") as f:
        json.dump(std, f, indent=2)
    with open(os.path.join(OUT_DIR, "lm_tp_3.json"), "w", encoding="utf-8") as f:
        json.dump(lm_tp3, f, indent=2)
    with open(os.path.join(OUT_DIR, "lm_tp_5.json"), "w", encoding="utf-8") as f:
        json.dump(lm_tp5, f, indent=2)
    with open(os.path.join(OUT_DIR, "lm_bg_3.json"), "w", encoding="utf-8") as f:
        json.dump(lm_tn3, f, indent=2)
    with open(os.path.join(OUT_DIR, "lm_bg_5.json"), "w", encoding="utf-8") as f:
        json.dump(lm_tn5, f, indent=2)
    with open(os.path.join(OUT_DIR, "motif_weights.json"), "w", encoding="utf-8") as f:
        json.dump(motif_weights, f, indent=2)
    # ---- 9) save tp/tn.tsv ----
    write_lm_tsv(lm_tp3, os.path.join(OUT_DIR, "tp_lm3.tsv"))
    write_lm_tsv(lm_tp5, os.path.join(OUT_DIR, "tp_lm5.tsv"))
    write_lm_tsv(lm_tn3, os.path.join(OUT_DIR, "tn_lm3.tsv"))
    write_lm_tsv(lm_tn5, os.path.join(OUT_DIR, "tn_lm5.tsv"))    
    print("[OK] wrote:", OUT_XLSX)
    print("[OK] artifacts in:", OUT_DIR)
if __name__ == "__main__":
    score_candidates_excel(require_5utr=False)
