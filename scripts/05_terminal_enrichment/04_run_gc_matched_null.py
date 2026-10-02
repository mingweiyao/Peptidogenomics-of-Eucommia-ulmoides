# -*- coding: utf-8 -*-
"""R2-Major10 follow-up step 3: GC-stratified-matched terminal-enrichment permutation.

Null: intervals (length preserved) may only land in non-coding space, and the GC content of the
hosting 10-kb bin must fall in the same GC decile as the observed interval's regional GC
(regional-GC stratified). All other conventions match step 1 (distal-10% windows per arm, within
chromosome, 2,000 permutations, empirical p=(>=obs+1)/2001). Observed/expected and window definitions follow the calibrated version.
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
_ap.add_argument("--rej", default=r"./data/rejected_20260204", help="upstream archive root")
_ap.add_argument("--gen13", default=r"./data/figure_middleware", help="figure middleware root")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

ROOT = os.path.abspath(os.path.expanduser(args.rej))
REJ = os.path.abspath(os.path.expanduser(args.rej))
RAW = REJ + r"/00raw_data"
FIG3 = REJ + r"/02figure/figure3"
DATA = OUTROOT + r"/_base_data_chain"
GLM = OUTROOT + r"/SupFig5_terminal_enrichment"
BIN = 10_000
NPERM = 2000
RNG = np.random.default_rng(20260906)

# Chromosome lengths + mapping
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

# NCP intervals
s4 = pd.read_csv(DATA + "/supplementary_table4_high_confidence.csv")
s4 = s4[[cmap.get(c) in L for c in s4["Chr"]]].copy()
s4["chrN"] = [cmap[c] for c in s4["Chr"]]
s4["s"] = s4["Start"].astype(int) - 1
s4["e"] = s4["End"].astype(int)
ncp = s4.drop_duplicates(["chrN", "s", "e"])[["chrN", "s", "e"]]

# GFF CDS -> non-coding runs (same as step 1)
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

# Genomic GC: read FASTA per chromosome, count GC per 10-kb bin; keep base-level cumulative GC (for QC)
bin_gc = {}     # chr -> np.array(GC count per bin)
bin_n = {}      # chr -> np.array(non-N length per bin, for GC%)
gc_cum = {}     # chr -> int32 cumulative GC counts (base level, QC)
fasta = open(RAW + "/Eu_genome.fasta", "rb")
cur = None
buf = bytearray()
def flush(cur, buf):
    if cur is None or cur not in L:
        return
    arr = np.frombuffer(bytes(buf), dtype=np.uint8)
    isgc = ((arr == ord("G")) | (arr == ord("C"))).astype(np.int64)
    isn = ((arr == ord("N")) | (arr == ord("n"))).astype(np.int64)
    nb = len(arr) // BIN + 1
    cgc = np.concatenate(([0], np.cumsum(isgc)))
    cnn = np.concatenate(([0], np.cumsum(isn)))
    g = np.array([cgc[min((b+1)*BIN, len(arr))] - cgc[b*BIN] for b in range(nb)])
    nn = np.array([cnn[min((b+1)*BIN, len(arr))] - cnn[b*BIN] for b in range(nb)])
    bin_gc[cur] = g
    bin_n[cur] = (np.minimum((np.arange(nb)+1)*BIN, len(arr)) - np.arange(nb)*BIN) - nn
    gc_cum[cur] = np.concatenate(([0], cgc)).astype(np.int64)
while True:
    chunk = fasta.read(1 << 24)
    if not chunk:
        break
    while True:
        i = chunk.find(b">")
        if i < 0:
            break
        if cur is not None:
            buf.extend(chunk[:i])
            flush(cur, buf)
        nl = chunk.find(b"\n", i)
        name = chunk[i+1:nl].split()[0].decode()
        cur = cmap.get(name)          # unmapped (scaffold) -> None, dropped
        chunk = chunk[nl+1:]
        buf = bytearray()
    if cur is not None:
        buf.extend(chunk.replace(b"\n", b""))
fasta.close()
flush(cur, buf)
print("GC bins done")

# GC% decile edges over all non-coding bins genome-wide
gc_all = []
for c in L:
    g, nn = bin_gc[c], bin_n[c]
    ok = nn > 0
    # count bins intersecting non-coding runs (lenient: non-N length>0 and intersecting any run)
    gc_all.append((g[ok] / nn[ok]))
gc_all = np.concatenate(gc_all)
edges = np.quantile(gc_all, np.arange(1, 10) / 10)
print("GC decile edges:", edges.round(4))

# Per chromosome: pieces (non-coding runs split by bin), observed regional-GC deciles, permutations
def windows_for(clen):
    w = clen * 10 // 100
    return [(0, w), (clen - w, clen)]

rows = []
cnt_list, obs_total, ntot_total = [], 0, 0
for c in sorted(L):
    clen = L[c]
    wins = windows_for(clen)
    sub = ncp[ncp["chrN"] == c]
    n = len(sub)
    if n == 0:
        continue
    g, nn = bin_gc[c], bin_n[c]
    ngc = np.divide(g, np.maximum(nn, 1))
    # pieces: intersect non-coding runs with bin edges
    piece_bin, piece_start, piece_len = [], [], []
    for rs, rl in nc_runs[c]:
        p = rs
        end = rs + rl
        while p < end:
            b = min(p // BIN, len(ngc) - 1)
            e = min((b + 1) * BIN, end)
            piece_bin.append(int(b))
            piece_start.append(p)
            piece_len.append(e - p)
            p = e
    piece_bin = np.array(piece_bin); piece_start = np.array(piece_start, dtype=np.int64)
    piece_len = np.array(piece_len, dtype=np.int64)
    # bin -> decile
    gcpc = np.divide(g, np.maximum(nn, 1))
    bdec = np.searchsorted(edges, gcpc, side="right")   # GC decile 0..9 per bin
    piece_dec = np.minimum(bdec[np.minimum(piece_bin, len(gcpc) - 1)], 9)  # GC decile per piece (capped at 9)
    # Observed: coordinates, regional-GC decile, terminal counts
    raw_dec = np.searchsorted(edges, gcpc[np.minimum(sub["s"].values // BIN, len(gcpc) - 1)], side="right")
    n_dec10 = int((raw_dec == 10).sum())
    print(f"    Diagnostic: {n_dec10}/{n} NCP start-bin GC exceeded the decile cap (dec=10); merged into the top group")
    dec_of = np.minimum(raw_dec, 9)
    obs = 0
    for s, e in zip(sub["s"], sub["e"]):
        for w1, w2 in wins:
            if s < w2 and e > w1:
                obs += 1
                break
    # Permute within GC deciles (length preserved; piece weight approximated by piece_len)
    cnt = np.zeros(NPERM, dtype=np.int64)
    Larr = (sub["e"] - sub["s"]).values.astype(np.int64)
    perm_pos = np.empty((NPERM, n), dtype=np.int64)   # keep placements for GC-match QC
    for d in range(10):
        sel = np.where(dec_of == d)[0]
        if len(sel) == 0:
            continue
        pm = piece_dec == d
        zs, zl = piece_start[pm], piece_len[pm]
        keep = zl > 0
        zs, zl = zs[keep], zl[keep]
        Li = Larr[sel]
        zi = RNG.integers(0, len(zs), (NPERM, len(sel)))
        zs2 = zs[zi]; zl2 = zl[zi]
        span = np.maximum(zl2 - Li[None, :] + 1, 1)
        pos = zs2 + (RNG.random((NPERM, len(sel))) * span).astype(np.int64)
        perm_pos[:, sel] = pos
        hit = np.zeros((NPERM, len(sel)), dtype=bool)
        for w1, w2 in wins:
            hit |= (pos < w2) & (pos + Li[None, :] > w1)
        cnt += hit.sum(axis=1)
    p = (int((cnt >= obs).sum()) + 1) / (NPERM + 1)
    cnt_list.append(cnt.copy()); obs_total += obs; ntot_total += n
    # QC: placed GC vs observed GC (first 500 permutations, all intervals, base-level cumulative GC)
    gcum = gc_cum[c]
    idx = RNG.choice(NPERM, size=min(500, NPERM), replace=False)
    pgc = (gcum[np.minimum(perm_pos[idx] + Larr[None, :], clen)] - gcum[perm_pos[idx]]) / Larr[None, :]
    ogc = (gcum[np.minimum(sub["s"].values + Larr, clen)] - gcum[sub["s"].values]) / Larr
    ks_p = stats.ks_2samp(pgc.ravel(), ogc).pvalue
    rows.append(dict(Chr=c, N_intervals=n, Observed=obs,
                     Null_mean_GCmatched=round(float(cnt.mean()), 1),
                     Empirical_p_GCmatched=p,
                     Placed_GC_mean=round(float(pgc.mean()), 4),
                     Observed_GC_mean=round(float(ogc.mean()), 4),
                     GCmatch_KS_p=round(float(ks_p), 4)))
    print(f"Chr{c:2d} obs={obs:4d} | GC-matched null={cnt.mean():6.1f} p={p:.4f}")

res = pd.DataFrame(rows)
pooled_cnt = np.sum(cnt_list, axis=0)
pooled_p = (int((pooled_cnt >= obs_total).sum()) + 1) / (NPERM + 1)
print(f"POOLED: observed {obs_total} | GC-matched null mean {pooled_cnt.mean():.1f} p={pooled_p:.4f}")
res.loc[len(res)] = dict(Chr='Pooled', N_intervals=ntot_total, Observed=obs_total,
                         Null_mean_GCmatched=round(float(pooled_cnt.mean()), 1),
                         Empirical_p_GCmatched=pooled_p,
                         Placed_GC_mean=np.nan, Observed_GC_mean=np.nan, GCmatch_KS_p=np.nan)
res.to_csv(GLM + "/Fig3C_terminal_enrichment_GCmatched.csv", index=False, encoding="utf-8-sig")
print("saved:", GLM + "/Fig3C_terminal_enrichment_GCmatched.csv")
print("GC-matched significant chromosomes (p<=0.05) =", int((res['Empirical_p_GCmatched'] <= 0.05).sum()),
      "| (p<=0.005) =", int((res['Empirical_p_GCmatched'] <= 0.005).sum()))
print("GC-match QC: placed vs observed GC mean (KS p>0.05 treated as consistent distributions)")
print(res[["Chr", "Placed_GC_mean", "Observed_GC_mean", "GCmatch_KS_p"]].to_string(index=False))
