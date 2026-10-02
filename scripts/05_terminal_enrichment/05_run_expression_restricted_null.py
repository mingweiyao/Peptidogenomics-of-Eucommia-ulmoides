# -*- coding: utf-8 -*-
"""R2-Major10 follow-up step 4: expression-restricted terminal-enrichment permutation.

Null: intervals (length preserved) may only land in non-coding space with transcriptional evidence - i.e.
merged intervals of features with CPM>=1 in any sample of total_sp_matrix_all.csv (223,870 non-overlapping
NCP features x 122 transcriptomes). All other conventions match step 1.
"""
import numpy as np
import pandas as pd
import re
import os as _os
import argparse
import os

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal", help="upstream archive root")
_ap.add_argument("--gen13", default=r"./data/figure_middleware", help="figure middleware root")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

ROOT = os.path.abspath(os.path.expanduser(args.rej))
REJ = os.path.abspath(os.path.expanduser(args.rej))
RNASEQ = REJ + r"/01location/rnaseq/03ouput"
FIG3 = REJ + r"/02figure/figure3"
DATA = OUTROOT + r"/_base_data_chain"
GLM = OUTROOT + r"/SupFig5_terminal_enrichment"
NPERM = 2000
RNG = np.random.default_rng(20260906)

L = {}
for line in open(FIG3 + "/figure3c_seq_len.txt", encoding="utf-8"):
    parts = line.split()
    if len(parts) == 3 and parts[0].isdigit():
        L[int(parts[0])] = int(parts[2])
cmap = {}
for line in open(os.path.join(args.gen13, '01figure/figure3/chrom_count.csv'), encoding="utf-8"):
    m = re.match(r"(GWHBISF\d+),(Chr\d+),", line)
    if m:
        cmap[m.group(1)] = int(m.group(2)[3:])

# NCP intervals (same as step 1)
s4 = pd.read_csv(DATA + "/supplementary_table4_high_confidence.csv")
s4 = s4[[cmap.get(c) in L for c in s4["Chr"]]].copy()
s4["chrN"] = [cmap[c] for c in s4["Chr"]]
s4["s"] = s4["Start"].astype(int) - 1
s4["e"] = s4["End"].astype(int)
ncp = s4.drop_duplicates(["chrN", "s", "e"])[["chrN", "s", "e"]]

# Expressed features: total_sp_matrix_all.csv, CPM>=1 in any sample
print("Reading expression matrix...")
m = pd.read_csv(RNASEQ + "/total_sp_matrix_all.csv", index_col=0)
cs = m.sum(axis=0).values.astype(float)
cpm = m.values / cs[None, :] * 1e6
expressed = (cpm >= 1).any(axis=1)
print(f"{m.shape[0]} features; with expression evidence (CPM>=1 in any sample): {int(expressed.sum())}")

iv = {"c": [], "s": [], "e": []}
for gid, isexp in zip(m.index[expressed], expressed[expressed]):
    parts = str(gid).split("_")
    if len(parts) < 3:
        continue
    c = cmap.get("GWHBISF" + parts[0].zfill(8))
    if c is None:
        continue
    iv["c"].append(c); iv["s"].append(int(parts[1]) - 1); iv["e"].append(int(parts[2]))
E = pd.DataFrame(iv)
print("Expressed-feature intervals mapped to Chr1-17:", len(E))

# Merge expressed intervals per chromosome -> restricted placement space
expr_runs = {}
for c in sorted(L):
    sub = E[E["c"] == c].sort_values("s")
    merged = []
    for a, b in zip(sub["s"], sub["e"]):
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    expr_runs[c] = [(int(a), int(b) - int(a)) for a, b in merged]

def windows_for(clen):
    w = clen * 10 // 100
    return [(0, w), (clen - w, clen)]

def overlap_stats(runs, ilen, windows):
    t = o = 0
    for rs, rl in runs:
        if rl < ilen:
            continue
        t += rl - ilen + 1
        lo_valid, hi_valid = rs, rs + rl - ilen
        for w1, w2 in windows:
            lo, hi = max(lo_valid, w1 - ilen + 1), min(hi_valid, w2 - 1)
            if hi >= lo:
                o += hi - lo + 1
    return t, o

rows = []
pooled_u = np.zeros(NPERM, dtype=np.int64)
pooled_e = np.zeros(NPERM, dtype=np.int64)
pooled_obs = 0
for c in sorted(L):
    clen = L[c]
    wins = windows_for(clen)
    sub = ncp[ncp["chrN"] == c]
    n = len(sub)
    if n == 0:
        continue
    obs = 0
    for s, e in zip(sub["s"], sub["e"]):
        for w1, w2 in wins:
            if s < w2 and e > w1:
                obs += 1
                break
    n_int = len(sub)
    t_u = np.empty(n_int); o_u = np.empty(n_int)
    t_e = np.empty(n_int); o_e = np.empty(n_int)
    for k, (s, e) in enumerate(zip(sub["s"], sub["e"])):
        ilen = e - s
        t_u[k] = clen - ilen + 1
        o_u[k] = sum(max(0, min(w2, clen - ilen + 1) - max(w1, 0)) for w1, w2 in wins)
        t_e[k], o_e[k] = overlap_stats(expr_runs[c], ilen, wins)
    cnt_u = (RNG.random((NPERM, n_int)) < o_u / t_u).sum(axis=1)
    cnt_e = (RNG.random((NPERM, n_int)) < o_e / t_e).sum(axis=1)
    p_u = (int((cnt_u >= obs).sum()) + 1) / (NPERM + 1)
    p_e = (int((cnt_e >= obs).sum()) + 1) / (NPERM + 1)
    pooled_u += cnt_u; pooled_e += cnt_e; pooled_obs += obs
    rows.append(dict(Chr=c, N_intervals=n_int, Observed=obs,
                     Null_mean_uniform=round(float(cnt_u.mean()), 1),
                     Null_mean_ExprRestricted=round(float(cnt_e.mean()), 1),
                     Empirical_p_ExprRestricted=p_e,
                     Availability_ratio=round(float(cnt_e.mean()) / max(cnt_u.mean(), 1e-9), 3)))
    print(f"Chr{c:2d} obs={obs:4d} | uniform={cnt_u.mean():6.1f} | ExprRestricted={cnt_e.mean():6.1f} p={p_e:.4f}")

res = pd.DataFrame(rows)
res.to_csv(GLM + "/Fig3C_terminal_enrichment_ExprRestricted.csv", index=False, encoding="utf-8-sig")
print("\nsaved:", GLM + "/Fig3C_terminal_enrichment_ExprRestricted.csv")
print("ExprRestricted significant chromosomes (p<=0.05) =", int((res['Empirical_p_ExprRestricted'] <= 0.05).sum()),
      "| (p<=0.005) =", int((res['Empirical_p_ExprRestricted'] <= 0.005).sum()))
pooled_p_e = (int((pooled_e >= pooled_obs).sum()) + 1) / (NPERM + 1)
pooled_p_u = (int((pooled_u >= pooled_obs).sum()) + 1) / (NPERM + 1)
print(f"POOLED: observed {pooled_obs} | uniform null mean {pooled_u.mean():.1f} p={pooled_p_u:.4f}"
      f" | ExprRestricted null mean {pooled_e.mean():.1f} p={pooled_p_e:.4f}")
