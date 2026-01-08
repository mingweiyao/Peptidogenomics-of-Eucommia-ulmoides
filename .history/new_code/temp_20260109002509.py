import os, re, json, random, math
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import multiprocessing as mp

# =========================================================
# CONFIG
# =========================================================
INPUT_DIR = "/data/Eu/Eu_genome"
CDS_FA = os.path.join(INPUT_DIR, "GWHBISF00000000.CDS.fasta")
GENOME_FA = os.path.join(INPUT_DIR, "GWHBISF00000000.genome.fasta")
CANDIDATES_XLSX = os.path.join(INPUT_DIR, "output_candidates.xlsx")
GFF3_FA = os.path.join(INPUT_DIR, "GWHBISF00000000.gff")
OUT_DIR = os.path.join(INPUT_DIR, "codon_prediction")
OUT_XLSX = os.path.join(OUT_DIR, "candidates_scored.xlsx")
SHEET_NAME = "output_candidates"

THREADS = 50
TP_ONLY_KMER_SHUFFLES = 200
TP_ONLY_KMER_MIN_Z = 2.5
TP_ONLY_KMER_ARTIFACT = os.path.join(OUT_DIR, "tp_only_kmer_shuffle_weights.json")

CODON_BONUS_SCHEMES = {
    "weak":   {"ATG": 0.0, "CTG": -0.1, "GTG": -0.2, "TTG": -0.2, "ACG": -0.3},
    "medium": {"ATG": 0.0, "CTG": -0.5, "GTG": -1.0, "TTG": -1.0, "ACG": -2.0},
    "strong": {"ATG": 0.0, "CTG": -1.0, "GTG": -2.0, "TTG": -2.0, "ACG": -3.0},
}
DEFAULT_OTHER_CODON_BONUS = -3.0

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
    g = genome_dict[chrom]
    if strand == "+":
        start = min(s for s, _ in exon_coords)
        if start - 101 < 0: return None
        up_seq = g[start - 101: start - 1] + cds_seq
    elif strand == "-":
        end = max(e for _, e in exon_coords)
        if end + 100 > len(g): return None
        up_seq = str(Seq(g[end: end + 100]).reverse_complement()) + cds_seq
    return up_seq
def load_cds_fasta(cds_file, gff3_file, genome_dict):
    # has_5utr = parse_gff_file(gff3_file)
    cds_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        tid = record.id
        # if tid not in has_5utr: continue
        merged = parse_tid_details(record.description, str(record.seq), genome_dict)
        if merged is None:
            continue
        cds_dict[tid] = merged
    return cds_dict    
# =========================================================
# Training data: TP/TN sampling (only used for LR training on context, not for motif discovery)
# =========================================================
def upstream_window(seq, tis_pos0):
    s = tis_pos0 - 100; e = tis_pos0 - 0
    if s < 0 or e > len(seq) or e <= s: return None
    up = seq[s:e]
    return up if len(up) == 100 else None
def sample_tp_tn_multi_seed(cds_dict, tn_per_tx=5, min_internal_nt=30, seeds=range(1,11), dedup_tn = True):
    rows_tp = []
    for tid, seq in cds_dict.items():
        if seq[100:103] != "ATG": continue
        rows_tp.append((tid, 100, seq, 1))
    tn_set = set()
    rows_tn = []        
    for seed in seeds:
        rng = random.Random(seed)
        for tid, seq in cds_dict.items():
            if seq[100:103] != "ATG": continue
            cand = []
            for i in range(103, len(seq) - 2, 3):
                if (i - 100) < min_internal_nt: continue
                if seq[i:i+3] == "ATG": cand.append(i)
            if not cand: continue
            picks = rng.sample(cand, k=min(tn_per_tx, len(cand)))
            for i in picks:
                if dedup_tn:
                    key = (tid, i)
                    if key in tn_set: continue
                    tn_set.add(key)
                rows_tn.append((tid, i, seq, 0))
    df = pd.DataFrame(rows_tp + rows_tn, columns=["tid", "tis_pos", "seq", "label"])
    return df
def compute_bg(cds_dict):
    c = {b: 0 for b in "ACGT"}
    for seq in cds_dict.values():
        for ch in seq[100:]:
            if ch in c: c[ch] += 1
    tot = sum(c.values())
    return {b: (c[b] / tot if tot else 0.25) for b in c}
def kozak_window(seq, tis_pos):
    s = tis_pos - 6
    e = tis_pos + 7
    if s < 0 or e > len(seq):
        return None
    w = seq[s:e]
    return w if len(w) == 13 else None
def build_pwm(tp_windows, pseudocount=1.0):
    if not tp_windows:
        raise ValueError("TP kozak windows 为空，无法建 PWM")
    L = len(tp_windows[0])
    counts = {b: np.full(L, pseudocount, dtype=float) for b in "ACGT"}
    for w in tp_windows:
        if len(w) != L: continue
        for i, ch in enumerate(w):
            if ch in counts:
                counts[ch][i] += 1.0
    totals = sum(counts[b] for b in "ACGT")
    return {b: (counts[b] / totals).tolist() for b in "ACGT"}
def kozak_logodds(w, pwm, bg, eps=1e-12):
    s = 0.0
    for i, ch in enumerate(w):
        if ch in "ACGT":
            ratio = pwm[ch][i] / max(bg[ch], eps)
            s += math.log2(max(ratio, eps))
    return s
def cu_fraction(up):
    return (sum(1 for ch in up if ch in "CT") / len(up)) if up else float("nan")
def tp_only_kmer_score(up, weights_z, k):
    if up is None or len(up) < k:
        return float("nan")
    up = up.upper()
    s = 0.0
    L = len(up)
    for i in range(L - k + 1):
        km = up[i:i+k]
        w = weights_z.get(km)
        if w is not None:
            s += w
    return float(s)
def train_pwm_and_weights(cds_dict, w3, w5):
    bg = compute_bg(cds_dict)
    df = sample_tp_tn_multi_seed(cds_dict)
    tp_ws = []
    for _, r in df[df["label"] == 1].iterrows():
        w = kozak_window(r["seq"], int(r["tis_pos"]))
        if w is not None:
            tp_ws.append(w)
    pwm = build_pwm(tp_ws)
    # LR features
    X = []; y = []
    for _, r in df.iterrows():
        seq = r["seq"]; tis = int(r["tis_pos"])
        w = kozak_window(seq, tis)
        up = upstream_window(seq, tis)
        if w is None or up is None: continue
        kz = kozak_logodds(w, pwm, bg)
        cu = cu_fraction(up)
        s3 = tp_only_kmer_score(up, w3, 3)
        s5 = tp_only_kmer_score(up, w5, 5)
        X.append([kz, cu, s3, s5])
        y.append(int(r["label"]))
    X = np.asarray(X, float)
    y = np.asarray(y, int)
    if X.shape[0] < 50: raise ValueError(f"训练样本过少：{X.shape[0]}")
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    lr = LogisticRegression(penalty="l2", solver="liblinear", max_iter=2000, class_weight="balanced")
    lr.fit(Z, y)
    weights = {
        "w1": float(lr.coef_[0][0]),  # kozak
        "w2": float(lr.coef_[0][1]),  # cu_fraction
        "w3": float(lr.coef_[0][2]),  # tp_only_3mer_score
        "w4": float(lr.coef_[0][3]),  # tp_only_5mer_score
        "beta0": float(lr.intercept_[0]),
    }
    std = {
        "mu": mu.tolist(),
        "sd": sd.tolist(),
        "feat": ["kozak", "cu_fraction", "tp_only_3mer_score", "tp_only_5mer_score"],
    }
    return pwm, bg, weights, std
# =========================================================
# Normalized background influence of TP (multiprocessing)
# =========================================================
def _eulerian_trail_from_adj(adj, start):
    local = {u: adj[u][:] for u in adj}
    for u in local: random.shuffle(local[u])
    stack = [start]
    path = []
    while stack:
        v = stack[-1]
        if local[v]: stack.append(local[v].pop())
        else: path.append(stack.pop())
    path.reverse()
    return path
def _build_dinucl_graph(seq):
    seq = seq.upper()
    adj = {b: [] for b in "ACGT"}
    for a, b in zip(seq[:-1], seq[1:]):
        if a in adj and b in "ACGT": adj[a].append(b)
        else: return None
    return adj
def dinuc_shuffle(seq, max_tries=20):
    seq = seq.upper()
    if len(seq) < 2 or any(ch not in "ACGT" for ch in seq): return None
    adj = _build_dinucl_graph(seq)
    if adj is None: return None
    start = seq[0]
    for _ in range(max_tries):
        path = _eulerian_trail_from_adj(adj, start)
        if len(path) == len(seq):
            return "".join(path)
    return None    
def count_kmers(seq, k):
    c = Counter()
    seq = seq.upper()
    L = len(seq)
    pat = re.compile(rf"[ACGT]{{{k}}}")
    for i in range(L - k + 1):
        s = seq[i:i+k]
        if pat.fullmatch(s): c[s] += 1
    return c
_G_TP_UPS = None
_G_K = None
_G_K2I = None
_G_KMERS_LEN = None
_G_DINUC_SHUFFLE = None
def _init_worker(tp_ups, k, k2i, kmers_len):
    global _G_TP_UPS, _G_K, _G_K2I, _G_KMERS_LEN, _G_DINUC_SHUFFLE
    _G_TP_UPS = tp_ups
    _G_K = k
    _G_K2I = k2i
    _G_KMERS_LEN = kmers_len
    _G_DINUC_SHUFFLE = dinuc_shuffle
def _one_shuffle_rep(seed):
    rng = random.Random(seed)
    c = np.zeros(_G_KMERS_LEN, dtype=np.float64)
    random.seed(rng.randint(0, 2**31 - 1))
    k = _G_K
    k2i = _G_K2I
    for up in _G_TP_UPS:
        sh = _G_DINUC_SHUFFLE(up)
        if sh is None: continue
        L = len(sh)
        for i in range(L - k + 1):
            km = sh[i:i+k]
            j = k2i.get(km)
            if j is not None:
                c[j] += 1.0
    return c
def learn_tp_only_kmer_weights(tp_ups, k, n_shuffle=200, min_z=2.5, seed=7, top_keep=None, n_jobs=THREADS):
    obs = Counter()
    for up in tp_ups:
        obs.update(count_kmers(up, k))
    kmers = [''.join(p) for p in __import__("itertools").product("ACGT", repeat=k)]
    k2i = {km: i for i, km in enumerate(kmers)}
    rep_seeds = [seed + 1000003 * r for r in range(n_shuffle)]
    n_jobs = max(1, min(n_jobs, n_shuffle))
    chunksize = 1 if n_shuffle <= 400 else max(1, n_shuffle // (n_jobs * 4))
    with mp.get_context("fork").Pool(
        processes=n_jobs, initializer=_init_worker,
        initargs=(tp_ups, k, k2i, len(kmers))) as pool:
        reps = pool.map(_one_shuffle_rep, rep_seeds, chunksize=chunksize)
    rep_sums = np.vstack(reps)
    exp = rep_sums.mean(axis=0)
    sd = rep_sums.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    obs_arr = np.zeros(len(kmers), dtype=np.float64)
    for km, i in k2i.items():
        obs_arr[i] = float(obs.get(km, 0))
    z = (obs_arr - exp) / sd
    weights_z = {km: float(z[k2i[km]]) for km in kmers if z[k2i[km]] >= min_z}
    if top_keep is not None and len(weights_z) > top_keep:
        weights_z = dict(sorted(weights_z.items(), key=lambda x: x[1], reverse=True)[:top_keep])
    return {
        "params": {"k": k, "n_shuffle": n_shuffle, "min_z": min_z, "seed": seed,
                   "window": f"[-{100}, -{0})", "up_len": 100, "top_keep": top_keep,
                   "n_jobs": n_jobs},
        "weights_z": weights_z,
    }
def ensure_tp_only_artifact(cds_dict):
    os.makedirs(os.path.dirname(TP_ONLY_KMER_ARTIFACT), exist_ok=True)
    if os.path.exists(TP_ONLY_KMER_ARTIFACT):
        with open(TP_ONLY_KMER_ARTIFACT, "r", encoding="utf-8") as f:
            return json.load(f)
    tp_ups = []
    for _, seq in cds_dict.items():
        if seq[100:103] != "ATG": continue
        up = upstream_window(seq, 100)
        if up is not None and re.fullmatch(r"[ACGT]{100}", up): tp_ups.append(up)
    if len(tp_ups) < 50: raise ValueError(f"TP upstream windows 太少（{len(tp_ups)}），无法建 TP-only kmer 背景模型")
    art3 = learn_tp_only_kmer_weights(tp_ups, k=3, n_shuffle=TP_ONLY_KMER_SHUFFLES, min_z=TP_ONLY_KMER_MIN_Z, seed=7, top_keep=None)
    art5 = learn_tp_only_kmer_weights(tp_ups, k=5, n_shuffle=TP_ONLY_KMER_SHUFFLES, min_z=TP_ONLY_KMER_MIN_Z, seed=11, top_keep=500)
    obj = {
        "params": {"window": f"[-{100}, -{0})", "up_len": 100},
        "3mer": art3,
        "5mer": art5,
    }
    with open(TP_ONLY_KMER_ARTIFACT, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
    return obj                  
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
def score_candidates_excel():
    os.makedirs(os.path.dirname(OUT_XLSX), exist_ok=True)
    genome_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(GENOME_FA, "fasta")}
    cds_dict = load_cds_fasta(CDS_FA, GFF3_FA, genome_dict)
    art = ensure_tp_only_artifact(cds_dict)
    w3 = art["3mer"]["weights_z"]
    w5 = art["5mer"]["weights_z"]
    pwm, bg, weights, std = train_pwm_and_weights(cds_dict, w3, w5)
    df_sp = pd.read_excel(CANDIDATES_XLSX, sheet_name=SHEET_NAME)
    mu = np.array(std["mu"], float)
    sd = np.array(std["sd"], float)
    candidate_keys = []
    kozak_scores = []
    cu_fracs = []
    tp3_scores = []
    tp5_scores = []
    tis_scores = {name: [] for name in CODON_BONUS_SCHEMES}   
    for _, r in df_sp.iterrows():
        chrom = str(r["chrom"])
        strand = str(r["strand"])
        start = int(r["phy_start"])
        end = int(r["phy_end"])
        codon = str(r.get("codon", "")).upper()
        candidate_keys.append(make_candidate_key(r))
        w = str(r.get("kozak", "")).upper()[:13]
        g = genome_dict[chrom]
        if strand == "+":
            left = start - 101; right = start - 1
            up = g[left:right] if (left >= 0 and right <= len(g)) else None
        else:
            left = end; right = end + 100
            up = str(Seq(g[left:right]).reverse_complement()) if (left >= 0 and right <= len(g)) else None
        kz = kozak_logodds(w, pwm, bg)
        cu = cu_fraction(up)
        s3 = tp_only_kmer_score(up, w3, 3)
        s5 = tp_only_kmer_score(up, w5, 5)
        z = (np.array([kz, cu, s3, s5], float) - mu) / sd
        # LR logit score
        score_base = (
            weights["beta0"]
            + weights["w1"] * z[0]
            + weights["w2"] * z[1]
            + weights["w3"] * z[2]
            + weights["w4"] * z[3]
        )
        kozak_scores.append(kz)
        cu_fracs.append(cu)
        tp3_scores.append(s3)
        tp5_scores.append(s5)
        for scheme_name, scheme_dict in CODON_BONUS_SCHEMES.items():
            tis_scores[scheme_name].append(score_base + codon_bonus_by_scheme(codon, scheme_dict))
    df_sp["candidate_key"] = candidate_keys
    df_sp["kozak_score"] = kozak_scores
    df_sp["cu_fraction"] = cu_fracs
    df_sp["tp_only_3mer_score"] = tp3_scores
    df_sp["tp_only_5mer_score"] = tp5_scores
    for scheme_name in CODON_BONUS_SCHEMES:
        score_col = f"tis_score_{scheme_name}"
        rank_col = f"rank_{scheme_name}"
        df_sp[score_col] = tis_scores[scheme_name]
        df_sp[rank_col] = df_sp.groupby("accession")[score_col].rank(ascending=False, method="first")
        top3 = df_sp[df_sp[rank_col] <= 3].sort_values(["accession", rank_col])
        top3.to_excel(OUT_XLSX.replace(".xlsx", f"_top3_{scheme_name}.xlsx"), index=False)
    df_sp.to_excel(OUT_XLSX, index=False)
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
    out_dir = os.path.dirname(OUT_XLSX)
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "pwm.json"), "w", encoding="utf-8") as f:
        json.dump(pwm, f, indent=2)
    with open(os.path.join(out_dir, "weights.json"), "w", encoding="utf-8") as f:
        json.dump(weights, f, indent=2)
    with open(os.path.join(out_dir, "standardizer.json"), "w", encoding="utf-8") as f:
        json.dump(std, f, indent=2)                    
if __name__ == "__main__":
    score_candidates_excel()