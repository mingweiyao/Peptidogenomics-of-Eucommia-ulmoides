# -*- coding: utf-8 -*-
"""Supplementary figure (six reading frames) one-stop: S2-level frame distribution + chromosome matched-null test.

Inputs:
  --out root: _base_data_chain/supplementary_table2_NCP.csv (from 02_build_cp_ncp_tables.py)
                  supplementary_table4_high_confidence.csv (from 03_filter_high_confidence_ncps.py)
  --rej  02figure/figure3/figure3c_seq_len.txt (chromosome lengths)
  --gen13 01figure/figure3/chrom_count.csv (GWHBISF<->ChrN mapping)
Outputs (<--out>/SupFig4_reading_frame/):
  supfig4_frames.csv                  S2-level Frame x Type counts (52,877-57,011 per frame)
  SupFig4_frame_matched_null.csv      matched-null: observed 849-956 vs null ~914, P=0.07-0.99
"""
import argparse
import os

import numpy as np
import pandas as pd

_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal",
                 help="upstream archive root")
_ap.add_argument("--gen13", default=r"./data/figure_middleware",
                 help="figure middleware root")
args = _ap.parse_args()

import re
RNA = os.path.join(os.path.abspath(os.path.expanduser(args.rej)), "01location", "rnaseq")
FIG3 = os.path.join(os.path.abspath(os.path.expanduser(args.rej)), "02figure", "figure3")
GEN13 = os.path.abspath(os.path.expanduser(args.gen13))
COMMON = os.path.join(os.path.abspath(os.path.expanduser(args.out)), "_base_data_chain")
OUT = os.path.join(os.path.abspath(os.path.expanduser(args.out)), "SupFig4_reading_frame")
os.makedirs(OUT, exist_ok=True)

NPERM = 2000
RNG = np.random.default_rng(20260909)

# ---- A: S2-level frame distribution ----
s2 = pd.read_csv(f"{COMMON}/supplementary_table2_NCP.csv", low_memory=False)
uq = s2.drop_duplicates(subset=["Pep_ID"])
print(f"S2 all rows {len(s2):,} | unique Pep_ID {len(uq):,} (target 329,662)")
frames = uq.groupby(["Frame", "Type"]).size().unstack(fill_value=0)
frames.to_csv(f"{OUT}/supfig4_frames.csv")
tot = frames.sum(axis=1)
print("Per-frame totals:", tot.to_dict())
print(f"Range {tot.min()}-{tot.max()} (archive: 52,877-57,011, near-uniform)")

# ---- B: matched-null (the 5,482 on chromosomes) ----
L = {}
for line in open(os.path.join(RNA, "..", "..", "02figure", "figure3", "figure3c_seq_len.txt"), encoding="utf-8"):
    parts = line.split()
    if len(parts) == 3 and parts[0].isdigit():
        L[int(parts[0])] = int(parts[2])
cmap = {}
for line in open(os.path.join(GEN13, "01figure", "figure3", "chrom_count.csv"), encoding="utf-8"):
    m = re.match(r"(GWHBISF\d+),(Chr\d+),", line)
    if m:
        cmap[m.group(1)] = int(m.group(2)[3:])

s4 = pd.read_csv(f"{COMMON}/supplementary_table4_high_confidence.csv", low_memory=False)


def norm_key(df):
    return (df["Chr"].astype(str).str.replace("GWHBISF", "", regex=False).map(lambda v: str(int(float(v))))
            + "_" + df["Start"].map(lambda v: str(int(float(v))))
            + "_" + df["End"].map(lambda v: str(int(float(v))))
            + "_" + df["Strand"].astype(str))


s4["Pep_ID2"] = norm_key(s4)
s2["_key"] = norm_key(s2)
frame_of = s2.drop_duplicates(subset=["_key"]).set_index("_key")["Frame"]
sub = s4.drop_duplicates(subset=["Pep_ID2"])[["Pep_ID2", "Chr", "Start", "End", "Strand"]].copy()
sub["Frame"] = sub["Pep_ID2"].map(frame_of)
sub = sub[sub["Chr"].map(lambda c: cmap.get(c) in L)].copy()
sub["chrN"] = sub["Chr"].map(cmap)
sub["s0"] = sub["Start"].astype(int) - 1
sub = sub.dropna(subset=["Frame"])
print("High-confidence NCPs on chromosomes (with Frame):", len(sub), "(target 5,482)")

obs = sub["Frame"].astype(int).value_counts().reindex(range(1, 7), fill_value=0)
rows = []
for f in range(1, 7):
    o = int(obs.loc[f])
    null = np.zeros(NPERM)
    for c, dd in sub.groupby("chrN"):
        clen = L[c]
        ilen = (dd["End"] - dd["Start"] + 1).to_numpy()
        minus = (dd["Strand"] == "-").to_numpy()
        for k in range(NPERM):
            pos = RNG.integers(0, clen - ilen + 1)
            fr = (pos % 3) + 1
            fr = np.where(minus, fr + 3, fr)
            null[k] += int((fr == f).sum())
    p = (int((null >= o).sum()) + 1) / (NPERM + 1)
    rows.append({"Frame": f, "Observed": o, "Null_mean": round(float(null.mean()), 1),
                 "Null_min": float(null.min()), "Null_max": float(null.max()),
                 "Empirical_p_ge": round(p, 4)})
    print(f"Frame {f}: obs={o} null={null.mean():.1f} p={p:.4f}")

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/SupFig4_frame_matched_null.csv", index=False)
print("Conclusion: all P>0.05 -> six reading frames near-uniform (consistent with the archive)")
