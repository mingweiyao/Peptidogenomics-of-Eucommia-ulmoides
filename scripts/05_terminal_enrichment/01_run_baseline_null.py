# -*- coding: utf-8 -*-
"""SupFig5-A: baseline terminal-enrichment permutation test (chromosome-matched, length-preserving, uniform null).

Method (main-text SupFig5 panel A / response letter R2-Major10; verified: the observed counts on all
17 chromosomes match the Observed column of the original Fig3C_terminal_enrichment_permutation.csv exactly):
  Observed = number of NCP intervals overlapping either distal-10% window (>=1 bp overlap counts);
  Null = 2,000 permutations: each NCP re-placed uniformly at random within its chromosome, length
  preserved (Bernoulli-equivalent method); empirical p = (#perms >= obs + 1) / 2001.
Result (archive): all 17 chromosomes significant (p<=0.005; observed 66-157 vs null 46-92).

Inputs: output/_base_data_chain/supplementary_table4_high_confidence.csv (S4),
      REJ/02figure/figure3/figure3c_seq_len.txt、GEN13/01figure/figure3/chrom_count.csv
Output: output/SupFig5_terminal_enrichment/Fig3C_terminal_enrichment_permutation.csv
      (column names identical to the original archive; read by 06_plot_terminal_enrichment.py)
Usage: py 01_run_baseline_null.py --out ./output --rej ./data/rejected_20260204 --gen13 ./data/figure_middleware
"""
import re
import numpy as np
import pandas as pd
import os
import os as _os
import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal", help="upstream archive root")
_ap.add_argument("--gen13", default=r"./data/figure_middleware", help="figure middleware root")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))


REJ = os.path.abspath(os.path.expanduser(args.rej))
GEN13 = os.path.abspath(os.path.expanduser(args.gen13))
COMMON = OUTROOT + "/_base_data_chain"
OUT = OUTROOT + "/SupFig5_terminal_enrichment"
os.makedirs(OUT, exist_ok=True)

NPERM = 2000
RNG = np.random.default_rng(20260906)

L = {}
for line in open(os.path.join(args.rej, '02figure/figure3/figure3c_seq_len.txt'), encoding="utf-8"):
    parts = line.split()
    if len(parts) == 3 and parts[0].isdigit():
        L[int(parts[0])] = int(parts[2])
cmap = {}
for line in open(os.path.join(args.gen13, '01figure/figure3/chrom_count.csv'), encoding="utf-8"):
    m = re.match(r"(GWHBISF\d+),(Chr\d+),", line)
    if m:
        cmap[m.group(1)] = int(m.group(2)[3:])

s4 = pd.read_csv(f"{COMMON}/supplementary_table4_high_confidence.csv", low_memory=False)
s4 = s4[[cmap.get(c) in L for c in s4["Chr"]]].copy()
s4["chrN"] = [cmap[c] for c in s4["Chr"]]
s4["s"] = s4["Start"].astype(int) - 1
s4["e"] = s4["End"].astype(int)
ncp = s4.drop_duplicates(["chrN", "s", "e"])
print("Chromosomal NCP intervals:", len(ncp), "(expect 5,482)")

rows = []
for c in sorted(L):
    clen = L[c]
    w = clen * 10 // 100
    wins = [(0, w), (clen - w, clen)]
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
    # per-interval fraction of valid start positions (uniform over the whole chromosome)
    ilen = (sub["e"] - sub["s"]).to_numpy()
    f_u = np.empty(len(ilen))
    for k, il in enumerate(ilen):
        t = clen - il + 1                       # total valid start positions
        o = 0                                   # valid starts overlapping a terminal window
        for w1, w2 in wins:
            lo, hi = max(w1 - il + 1, 0), min(w2 - 1, clen - il)
            if hi >= lo:
                o += hi - lo + 1
        f_u[k] = o / t
    cnt = (RNG.random((NPERM, len(ilen))) < f_u).sum(axis=1)
    p = (int((cnt >= obs).sum()) + 1) / (NPERM + 1)
    rows.append({"Chr": c, "Observed_in_terminal10pct": obs,
                 "Null_mean": round(float(cnt.mean()), 1), "Empirical_p_gt": round(p, 4)})
    print(f"Chr{c:2d} N={n:4d} obs={obs:4d} null={cnt.mean():6.1f} p={p:.4f}")

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/Fig3C_terminal_enrichment_permutation.csv", index=False)
print("\nSignificant chromosomes (p<=0.005):", int((res["Empirical_p_gt"] <= 0.005).sum()), "/ 17 (archive: all 17)")
print("saved:", OUT + "/Fig3C_terminal_enrichment_permutation.csv")
