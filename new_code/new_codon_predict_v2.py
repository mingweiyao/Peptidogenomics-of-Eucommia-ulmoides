import os, re, json, random, math
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

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

# 输出 upstream kmer 模型（kmer -> logP）
UP_LM_3 = os.path.join(OUT_DIR, "upstream_3mer_logp.tsv")
UP_LM_5 = os.path.join(OUT_DIR, "upstream_5mer_logp.tsv")
UP_LM_JSON = os.path.join(OUT_DIR, "upstream_kmer_lm.json")

# upstream window: [-100, 0) relative to codon first base (tis_pos0)
UP_LEN = 100
UP_START = 100
UP_END = 0

# k-mer LM settings
LM_KS = (3, 5)
LM_ALPHA = 1.0 
LM_LOG_BASE = math.e

# TN sampling (只用于 LR 学权重；你仍可用 multi-seed 来弱化抽样偶然性)
TN_PER_TX = 5
MIN_INTERNAL_NT = 30
TN_SEEDS = range(1, 11)
DEDUP_TN = True

CODON_BONUS_SCHEMES = {
    "weak":   {"ATG": 0.0, "CTG": -0.1, "GTG": -0.2, "TTG": -0.2, "ACG": -0.3},
    "medium": {"ATG": 0.0, "CTG": -0.5, "GTG": -1.0, "TTG": -1.0, "ACG": -2.0},
    "strong": {"ATG": 0.0, "CTG": -1.0, "GTG": -2.0, "TTG": -2.0, "ACG": -3.0},
}
DEFAULT_OTHER_CODON_BONUS = -3.0

DNA_RE_100 = re.compile(r"^[ACGT]{100}$")

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

    if chrom not in genome_dict:
        return None
    g = genome_dict[chrom]
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
def load_cds_fasta(cds_file, gff3_file, genome_dict, require_5utr=True):
    has_5utr = parse_gff_file(gff3_file) if require_5utr else None
    cds_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        tid = record.id
        if tid not in has_5utr: continue
        merged = parse_tid_details(record.description, str(record.seq), genome_dict)
        if merged is None: continue
        cds_dict[tid] = merged
    return cds_dict
# =========================================================
# Upstream k-mer Language Model (NEW)
# =========================================================
def count_kmers(seq, k):
    c = Counter()
    seq = seq.upper()
    L = len(seq)
    for i in range(L - k + 1):
        km = seq[i:i+k]
        if all(ch in "ACGT" for ch in km): c[km] += 1
    return c
def build_upstream_kmer_lm(tp_ups, k, alpha=1.0):
    obs = Counter()
    total = 0
    for up in tp_ups:
        cnt = count_kmers(up, k)
        obs.update(cnt)
        total += max(0, len(up) - k + 1)
    V = 4 ** k
    denom = total + alpha * V
    default_p = alpha / denom
    default_logp = math.log(default_p)
    logp = {km: math.log((c + alpha) / denom) for km, c in obs.items()}
    return {
        "params": {"k": k, "alpha": alpha, "total_windows": len(tp_ups), "total_kmer_positions": total, "V": V},
        "logp": logp,
        "default_logp": float(default_logp),
    }
def collect_tp_upstreams(cds_dict):
    tp_ups = []
    for _, seq in cds_dict.items():
        seq = seq.upper()
        if seq[UP_START:UP_START+3] != "ATG": continue
        up = upstream_window(seq, UP_START)
        if up is not None and DNA_RE_100.fullmatch(up): tp_ups.append(up)
    if len(tp_ups) < 50: raise ValueError(f"TP upstream windows 太少（{len(tp_ups)}），无法建 upstream k-mer 模型")
    return tp_ups
# =========================================================
# LR training (features now use upstream 3/5-mer LL)
# =========================================================
def sample_tp_tn_multi_seed(cds_dict, tn_per_tx=TN_PER_TX, min_internal_nt=MIN_INTERNAL_NT,
                            seeds=TN_SEEDS, dedup_tn=DEDUP_TN):
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
def build_pwm(tp_windows, pseudocount=1.0):
    if not tp_windows: raise ValueError("TP kozak windows 为空，无法建 PWM")
    counts = {b: np.full(13, pseudocount, dtype=float) for b in "ACGT"}
    for w in tp_windows:
        if len(w) != 13: continue
        for i, ch in enumerate(w):
            if ch in counts: counts[ch][i] += 1.0
    totals = sum(counts[b] for b in "ACGT")
    return {b: (counts[b] / totals).tolist() for b in "ACGT"}
def kozak_window(seq, tis_pos):
    s = tis_pos - 6
    e = tis_pos + 7
    if s < 0 or e > len(seq): return None
    w = seq[s:e].upper()
    return w if len(w) == 13 else None
def upstream_window(seq, tis_pos0):
    s = tis_pos0 - UP_START
    e = tis_pos0 - UP_END
    if s < 0 or e > len(seq) or e <= s: return None
    up = seq[s:e].upper()
    return up if len(up) == UP_LEN else None
def kozak_logprob(w, pwm, eps=1e-12):
    s = 0.0
    for i, ch in enumerate(w):
        if ch in "ACGT":
            s += math.log(max(pwm[ch][i], eps))
    return float(s)
def cu_fraction(up):
    return (sum(1 for ch in up if ch in "CT") / len(up)) if up else float("nan")
def score_upstream_ll(up, lm_obj, k):
    if up is None or len(up) < k: return float("nan")
    up = up.upper()
    L = len(up)
    denom = (L - k + 1)
    if denom <= 0: return float("nan")
    logp = lm_obj["logp"]
    dlog = lm_obj["default_logp"]
    s = 0.0
    for i in range(denom):
        km = up[i:i+k]
        if all(ch in "ACGT" for ch in km): s += logp.get(km, dlog)
        else: return float("nan")
    return float(s / denom)
def train_pwm_and_lr(cds_dict, lm3, lm5):
    df = sample_tp_tn_multi_seed(cds_dict)
    tp_ws = []
    for _, r in df[df["label"] == 1].iterrows():
        w = kozak_window(r["seq"], int(r["tis_pos"]))
        if w is not None: tp_ws.append(w)
    pwm = build_pwm(tp_ws)
    X = []
    y = []
    for _, r in df.iterrows():
        seq = r["seq"].upper()
        tis = int(r["tis_pos"])
        w = kozak_window(seq, tis)
        up = upstream_window(seq, tis)
        if w is None or up is None: continue
        if not DNA_RE_100.fullmatch(up): continue
        kz = kozak_logprob(w, pwm)
        cu = cu_fraction(up)
        ll3 = score_upstream_ll(up, lm3, 3)
        ll5 = score_upstream_ll(up, lm5, 5)
        if any(np.isnan(v) for v in [kz, cu, ll3, ll5]): continue
        X.append([kz, cu, ll3, ll5])
        y.append(int(r["label"]))
    X = np.asarray(X, float)
    y = np.asarray(y, int)
    if X.shape[0] < 50: raise ValueError(f"训练样本过少：{X.shape[0]}")
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    lr = LogisticRegression(
        penalty="l2",
        solver="liblinear",
        max_iter=2000,
        class_weight="balanced",
    )
    lr.fit(Z, y)
    weights = {
        "beta0": float(lr.intercept_[0]),
        "w1": float(lr.coef_[0][0]),
        "w2": float(lr.coef_[0][1]),
        "w3": float(lr.coef_[0][2]),
        "w4": float(lr.coef_[0][3]),
    }
    std = {
        "mu": mu.tolist(),
        "sd": sd.tolist(),
        "feat": ["kozak_logprob", "cu_fraction", "upstream_3mer_ll", "upstream_5mer_ll"],
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
def score_candidates_excel(require_5utr=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    genome_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(GENOME_FA, "fasta")}
    cds_dict = load_cds_fasta(CDS_FA, GFF3_FA, genome_dict, require_5utr=require_5utr)
    # ---- 1) build upstream LM from TP-only upstream windows ----
    tp_ups = collect_tp_upstreams(cds_dict)
    lm3 = build_upstream_kmer_lm(tp_ups, k=3, alpha=LM_ALPHA)
    lm5 = build_upstream_kmer_lm(tp_ups, k=5, alpha=LM_ALPHA)
    # save LM artifacts (json + tsv)
    with open(UP_LM_JSON, "w", encoding="utf-8") as f: json.dump({"3mer": lm3, "5mer": lm5}, f, indent=2)
    df3 = pd.DataFrame(lm3["logp"].items(), columns=["3mer", "logp"]).sort_values("logp", ascending=False)
    df5 = pd.DataFrame(lm5["logp"].items(), columns=["5mer", "logp"]).sort_values("logp", ascending=False)
    df3.to_csv(UP_LM_3, sep="\t", index=False)
    df5.to_csv(UP_LM_5, sep="\t", index=False)
    # ---- 2) train PWM + LR weights ----
    pwm, weights, std = train_pwm_and_lr(cds_dict, lm3, lm5)
    # ---- 3) read candidates ----
    df_sp = pd.read_excel(CANDIDATES_XLSX, sheet_name=SHEET_NAME)
    mu = np.array(std["mu"], float)
    sd = np.array(std["sd"], float)
    candidate_keys = []
    kozak_scores = []
    cu_fracs = []
    ll3_scores = []
    ll5_scores = []
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
        if strand == "+":
            left = start - (UP_START + 1); right = start - 1
            up = g[left:right] if (left >= 0 and right <= len(g)) else None
        else:
            left = end; right = end + UP_START
            up = str(Seq(g[left:right]).reverse_complement()) if (left >= 0 and right <= len(g)) else None
        if up is None or len(up) != UP_LEN or (not DNA_RE_100.fullmatch(up)):
            kozak_scores.append(np.nan); cu_fracs.append(np.nan); ll3_scores.append(np.nan); ll5_scores.append(np.nan)
            for name in CODON_BONUS_SCHEMES:
                tis_scores[name].append(np.nan)
            continue
        kz = kozak_logprob(w, pwm)
        cu = cu_fraction(up)
        ll3 = score_upstream_ll(up, lm3, 3)
        ll5 = score_upstream_ll(up, lm5, 5)
        z = (np.array([kz, cu, ll3, ll5], float) - mu) / sd
        score_base = (
            weights["beta0"]
            + weights["w1"] * z[0]
            + weights["w2"] * z[1]
            + weights["w3"] * z[2]
            + weights["w4"] * z[3]
        )
        kozak_scores.append(kz)
        cu_fracs.append(cu)
        ll3_scores.append(ll3)
        ll5_scores.append(ll5)
        for scheme_name, scheme_dict in CODON_BONUS_SCHEMES.items():
            tis_scores[scheme_name].append(score_base + codon_bonus_by_scheme(codon, scheme_dict))
    # ---- 4) write outputs ----
    df_sp["candidate_key"] = candidate_keys
    df_sp["kozak_score"] = kozak_scores
    df_sp["cu_fraction"] = cu_fracs
    df_sp["upstream_3mer_ll"] = ll3_scores
    df_sp["upstream_5mer_ll"] = ll5_scores
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
    # save model artifacts
    with open(os.path.join(OUT_DIR, "pwm.json"), "w", encoding="utf-8") as f:
        json.dump(pwm, f, indent=2)
    with open(os.path.join(OUT_DIR, "weights.json"), "w", encoding="utf-8") as f:
        json.dump(weights, f, indent=2)
    with open(os.path.join(OUT_DIR, "standardizer.json"), "w", encoding="utf-8") as f:
        json.dump(std, f, indent=2)
    print("[OK] wrote:", OUT_XLSX)
    print("[OK] wrote:", UP_LM_3)
    print("[OK] wrote:", UP_LM_5)
    print("[OK] wrote:", UP_LM_JSON)
if __name__ == "__main__":
    # require_5utr=True 更严格；如果你确认 gff 的 five_prime_UTR 很可靠，建议 True
    score_candidates_excel(require_5utr=True)
