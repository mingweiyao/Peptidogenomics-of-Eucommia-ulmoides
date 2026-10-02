# -*- coding: utf-8 -*-
"""R2-Major10 follow-up step 3 (v2): interval-level GC-exact-matched terminal-enrichment permutation.

Null: each NCP interval (length preserved) is re-placed into non-coding space such that the L-mer GC
content at the placed site deviates from the interval's own GC by <= tolerance (rejection sampling;
tolerance starts at 0.03 and is relaxed stepwise, falling back to the closest candidate). All other
conventions match step 1 (distal-10% windows per arm, within chromosome, 2,000 permutations, p=(>=obs+1)/2001).
"""
import numpy as np
import pandas as pd
import re
from scipy import stats
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

s4 = pd.read_csv(DATA + "/supplementary_table4_high_confidence.csv")
s4 = s4[[cmap.get(c) in L for c in s4["Chr"]]].copy()
s4["chrN"] = [cmap[c] for c in s4["Chr"]]
s4["s"] = s4["Start"].astype(int) - 1
s4["e"] = s4["End"].astype(int)
ncp = s4.drop_duplicates(["chrN", "s", "e"])[["chrN", "s", "e"]]

# GFF CDS -> non-coding runs
cds = {c: [] for c in L}
with open(RAW + "/GWHBISF00000000.gff", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#") or "\tCDS\t" not in line:
            continue
        fs = line.rstrip("\n").split("\t")
        c = cmap.get(fs[0])
        if c:
            cds[c].append((int(fs[3]) - 1, int(fs[4])))
nc_runs = {}
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

# FASTA -> per-chromosome base-level cumulative GC counts
gc_cum = {}
fasta = open(RAW + "/Eu_genome.fasta", "rb")
cur, buf = None, bytearray()
def flush(cur, buf):
    if cur is None or cur not in L:
        return
    arr = np.frombuffer(bytes(buf), dtype=np.uint8)
    isgc = ((arr == ord("G")) | (arr == ord("C"))).astype(np.int64)
    gc_cum[cur] = np.concatenate(([0], np.cumsum(isgc))).astype(np.int64)
while True:
    chunk = fasta.read(1 << 24)
    if not chunk:
        break
    while True:
        i = chunk.find(b">")
        if i < 0:
            break
        if cur is not None:
            buf.extend(chunk[:i]); flush(cur, buf)
        nl = chunk.find(b"\n", i)
        name = chunk[i+1:nl].split()[0].decode()
        cur = cmap.get(name)
        chunk = chunk[nl+1:]
        buf = bytearray()
    if cur is not None:
        buf.extend(chunk.replace(b"\n", b""))
fasta.close()
flush(cur, buf)
print("FASTA cumulative GC done")

def windows_for(clen):
    w = clen * 10 // 100
    return [(0, w), (clen - w, clen)]

TOL_SCHEDULE = [0.03, 0.05, 0.08, 0.12, 0.20]
BATCH = 50_000
MAXBATCH = 80

rows = []
pooled = np.zeros(NPERM, dtype=np.int64)
pooled_obs = 0
tol_report = []
for c in sorted(L):
    clen = L[c]
    wins = windows_for(clen)
    sub = ncp[ncp["chrN"] == c]
    n = len(sub)
    if n == 0:
        continue
    gcum = gc_cum[c]
    runs = nc_runs[c]
    rs = np.array([a for a, b in runs], dtype=np.int64)
    rl = np.array([b for a, b in runs], dtype=np.int64)
    keep = rl > 0
    rs, rl = rs[keep], rl[keep]
    cum = np.cumsum(rl)
    total = cum[-1]

    s_arr = sub["s"].values.astype(np.int64)
    e_arr = sub["e"].values.astype(np.int64)
    Larr = e_arr - s_arr
    gobs = (gcum[e_arr] - gcum[s_arr]) / Larr

    obs = 0
    for s, e in zip(s_arr, e_arr):
        for w1, w2 in wins:
            if s < w2 and e > w1:
                obs += 1
                break

    cnt = np.zeros(NPERM, dtype=np.int64)
    placed_gc_all = []
    tol_used_all = []
    for k in range(n):
        Li = int(Larr[k])
        gi = float(gobs[k])
        pp = None
        for tol in TOL_SCHEDULE:
            got = 0
            acc_pos, acc_d = [], []
            for _ in range(MAXBATCH):
                r = RNG.random(BATCH) * total
                ridx = np.searchsorted(cum, r, side="right")
                ridx = ridx[ridx < len(rs)]
                pos = rs[ridx] + (RNG.random(ridx.size) * (rl[ridx] - Li + 1)).astype(np.int64)
                ok = (rl[ridx] >= Li) & (pos + Li <= rs[ridx] + rl[ridx])
                pos = pos[ok]
                if pos.size == 0:
                    continue
                g = (gcum[pos + Li] - gcum[pos]) / Li
                d = np.abs(g - gi)
                acc = d <= tol
                if acc.any():
                    acc_pos.append(pos[acc]); acc_d.append(d[acc])
                    got += int(acc.sum())
                if got >= NPERM:
                    break
            if got >= NPERM:
                pp = np.concatenate(acc_pos)[:NPERM].astype(np.int64)
                tol_used_all.append(tol)
                break
        if pp is None:
            # fallback: take the candidate closest to target at the tightest tolerance
            if acc_pos:
                allpos = np.concatenate(acc_pos); alld = np.concatenate(acc_d)
                pp = allpos[np.argsort(alld)[:NPERM]].astype(np.int64)
            else:
                pp = np.zeros(0, dtype=np.int64)
            tol_used_all.append(TOL_SCHEDULE[-1])
        placed_gc_all.append((gcum[pp + Li] - gcum[pp]) / Li)
        hit = np.zeros(len(pp), dtype=np.int64)
        for w1, w2 in wins:
            hit += ((pp < w2) & (pp + Li > w1)).astype(np.int64)
        add = np.zeros(NPERM, dtype=np.int64)
        add[: len(hit)] = hit
        cnt += add
    p = (int((cnt >= obs).sum()) + 1) / (NPERM + 1)
    pooled += cnt
    pooled_obs += obs
    # QC: placed GC vs observed GC
    allpgc = np.concatenate(placed_gc_all)
    ks = stats.ks_2samp(allpgc, gobs).pvalue
    rows.append(dict(Chr=c, N_intervals=n, Observed=obs,
                     Null_mean_GCexact=round(float(cnt.mean()), 1),
                     Empirical_p_GCexact=p,
                     Placed_GC_mean=round(float(allpgc.mean()), 4),
                     Observed_GC_mean=round(float(gobs.mean()), 4),
                     GCmatch_KS_p=round(float(ks), 4),
                     Max_tol=float(max(tol_used_all))))
    print(f"Chr{c:2d} obs={obs:4d} | GC-exact null={cnt.mean():7.1f} p={p:.4f} | "
          f"placedGC={allpgc.mean():.3f} obsGC={gobs.mean():.3f} KS_p={ks:.3f} | "
          f"max tolerance={max(tol_used_all)}")

res = pd.DataFrame(rows)
res.to_csv(GLM + "/Fig3C_terminal_enrichment_GCexact.csv", index=False, encoding="utf-8-sig")
pooled_p = (int((pooled >= pooled_obs).sum()) + 1) / (NPERM + 1)
print("\nsaved:", GLM + "/Fig3C_terminal_enrichment_GCexact.csv")
print(f"GC-exact significant chromosomes (p<=0.05) = {int((res['Empirical_p_GCexact'] <= 0.05).sum())}"
      f" | (p<=0.005) = {int((res['Empirical_p_GCexact'] <= 0.005).sum())}")
print(f"POOLED: observed {pooled_obs} | GC-exact null mean {pooled.mean():.1f} p={pooled_p:.4f}")
