# -*- coding: utf-8 -*-
"""R2-Major10 follow-up: non-CDS-restricted terminal-enrichment permutation test.

Background: the published terminal-enrichment permutation (Fig3C_terminal_enrichment_permutation.csv) placed NCP
intervals uniformly along whole chromosomes, without controlling for gene density / non-coding space
availability - chromosome termini are gene-poor and non-coding-rich, so NCPs (which must sit in
non-coding sequence) may cluster there naturally. This script narrows the null sampling space from whole
chromosomes to the non-CDS complement (intergenic+intron+UTR space). All other conventions (distal-10%
windows per arm, length-preserving, within-chromosome, 2,000 permutations, empirical p=(>=obs+1)/2001)
match the original analysis, which is first recalibrated by reproducing the published numbers under the uniform null.

Inputs (read-only):
  00raw_data/GWHBISF00000000.gff        (CDS intervals -> define the non-coding space)
  00raw_data/Eu_genome.fasta            (chromosome lengths)
  reproduction/data/supplementary_table4_high_confidence.csv (NCP coordinates, S4)
  figure3/chrom_count.csv               (GWHBISF <-> ChrN mapping)
  figure3/figure3c_seq_len.txt          (chromosome lengths)
  output/SupFig5_terminal_enrichment/Fig3C_terminal_enrichment_permutation.csv (original result, for calibration)
Outputs:
  output/SupFig5_terminal_enrichment/Fig3C_terminal_enrichment_permutation_NCrestricted.csv
"""
import re
import numpy as np
import pandas as pd
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
RAW = REJ + r"/00raw_data"
FIG3 = REJ + r"/02figure/figure3"
DATA = OUTROOT + r"/_base_data_chain"
OUT = REJ + r"/02figure/figure6"  # placeholder; actual output goes to <out>/SupFig5_terminal_enrichment
GLM = OUTROOT + r"/SupFig5_terminal_enrichment"

NPERM = 2000
RNG = np.random.default_rng(20260906)

# ---------- 1. Chromosome lengths and mapping ----------
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

# ---------- 2. NCP intervals (S4 deduplicated, Chr1-17 only) ----------
s4 = pd.read_csv(DATA + "/supplementary_table4_high_confidence.csv")
s4 = s4[[cmap.get(c) in L for c in s4["Chr"]]].copy()
s4["chrN"] = [cmap[c] for c in s4["Chr"]]
s4["s"] = s4["Start"].astype(int) - 1          # 0-based
s4["e"] = s4["End"].astype(int)                 # half-open
ncp = s4.drop_duplicates(["chrN", "s", "e"])[["chrN", "s", "e"]]
print("Chromosomal NCP intervals:", len(ncp), "(expect 5,482)")

# ---------- 3. GFF CDS -> non-coding complement ----------
cds = {c: [] for c in L}
with open(RAW + "/GWHBISF00000000.gff", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#") or "\tCDS\t" not in line:
            continue
        fs = line.rstrip("\n").split("\t")
        c = cmap.get(fs[0])
        if c:
            cds[c].append((int(fs[3]) - 1, int(fs[4])))
nc_runs = {}   # chr -> list of (start, len) non-coding runs
for c, clen in L.items():
    iv = sorted(cds[c])
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    runs, pos = [], 0
    for a, b in merged:
        if a > pos:
            runs.append((pos, a - pos))
        pos = max(pos, b)
    if pos < clen:
        runs.append((pos, clen - pos))
    nc_runs[c] = runs
print("GFF CDS parsed; example Chr1 non-coding runs =", len(nc_runs[1]))

# ---------- 4. Permutation statistic helpers ----------
def overlap_stats(runs, ilen, windows):
    """Return (t, o): total valid start positions / valid starts overlapping the windows.
    A start p is valid when [p, p+ilen) lies fully inside one run. windows=[(w1,w2),...] are the overlap windows."""
    t = 0
    o = 0
    for rs, rl in runs:
        if rl < ilen:
            continue
        t += rl - ilen + 1
        lo_valid, hi_valid = rs, rs + rl - ilen        # closed interval of valid starts
        for w1, w2 in windows:
            lo = w1 - ilen + 1
            hi = w2 - 1
            lo2, hi2 = max(lo_valid, lo), min(hi_valid, hi)
            if hi2 >= lo2:
                o += hi2 - lo2 + 1
    return t, o

def windows_for(clen):
    w = clen * 10 // 100
    return [(0, w), (clen - w, clen)]               # distal 10% of each arm (both ends)

# ---------- 5. Per-chromosome computation under both nulls ----------
rows = []
pooled_u = np.zeros(NPERM, dtype=np.int64)
pooled_n = np.zeros(NPERM, dtype=np.int64)
pooled_obs_u = pooled_obs_n = 0
for c in sorted(L):
    clen = L[c]
    wins = windows_for(clen)
    sub = ncp[ncp["chrN"] == c]
    n = len(sub)
    if n == 0:
        continue
    # Observed
    obs = 0
    for s, e in zip(sub["s"], sub["e"]):
        ilen = e - s
        for w1, w2 in wins:
            if s < w2 and e > w1:
                obs += 1
                break
    # per-interval (valid starts t, overlapping starts o)
    t_u = np.empty(n); o_u = np.empty(n)      # uniform: whole chromosome valid
    t_n = np.empty(n); o_n = np.empty(n)      # noncoding-restricted
    for k, (s, e) in enumerate(zip(sub["s"], sub["e"])):
        ilen = e - s
        t_u[k] = clen - ilen + 1
        o_u[k] = sum(max(0, min(w2, clen - ilen + 1) - max(w1, 0)) for w1, w2 in wins)
        t_n[k], o_n[k] = overlap_stats(nc_runs[c], ilen, wins)
    f_u = o_u / t_u
    f_n = o_n / t_n
    # Bernoulli sum = permutation distribution equivalent to independent uniform sampling
    cnt_u = (RNG.random((NPERM, n)) < f_u).sum(axis=1)
    cnt_n = (RNG.random((NPERM, n)) < f_n).sum(axis=1)
    p_u = (int((cnt_u >= obs).sum()) + 1) / (NPERM + 1)
    p_n = (int((cnt_n >= obs).sum()) + 1) / (NPERM + 1)
    pooled_u += cnt_u; pooled_n += cnt_n
    pooled_obs_u += obs; pooled_obs_n += obs
    rows.append(dict(Chr=c, N_intervals=n, Chr_len=clen, Observed=obs,
                     Null_mean_uniform=round(float(cnt_u.mean()), 1),
                     Empirical_p_uniform=p_u,
                     Null_mean_NCrestricted=round(float(cnt_n.mean()), 1),
                     Empirical_p_NCrestricted=p_n,
                     Availability_ratio=round(float(cnt_n.mean()) / max(cnt_u.mean(), 1e-9), 3)))
    print(f"Chr{c:2d} N={n:4d} obs={obs:4d} | uniform null={cnt_u.mean():6.1f} p={p_u:.4f} | "
          f"NC-restricted null={cnt_n.mean():6.1f} p={p_n:.4f}")

res = pd.DataFrame(rows)
res.to_csv(GLM + "/Fig3C_terminal_enrichment_permutation_NCrestricted.csv", index=False, encoding="utf-8-sig")
print("\nsaved:", GLM + "/Fig3C_terminal_enrichment_permutation_NCrestricted.csv")
print("\nSummary: uniform-significant chromosomes =", int((res['Empirical_p_uniform'] <= 0.005).sum()),
      "| NC-restricted-significant chromosomes =", int((res['Empirical_p_NCrestricted'] <= 0.005).sum()))
pooled_p_u = (int((pooled_u >= pooled_obs_u).sum()) + 1) / (NPERM + 1)
pooled_p_n = (int((pooled_n >= pooled_obs_n).sum()) + 1) / (NPERM + 1)
print(f"POOLED (17 chromosomes): observed {pooled_obs_n} | uniform null mean {pooled_u.mean():.1f} p={pooled_p_u:.4f}"
      f" | NC-restricted null mean {pooled_n.mean():.1f} p={pooled_p_n:.4f}")
