# -*- coding: utf-8 -*-
"""Low-complexity flagging (R2-Major3 / R2-Minor8): flags low-complexity sequences among the 5,684 high-confidence candidates.

Rule (response-letter convention): entropy < 2.5 bits or top-residue fraction >= 50% -> low_complexity = True.
  entropy_bits = -sum p_i * log2(p_i) over the 20 amino-acid frequencies
  top1_frac = fraction of the most frequent residue
Verification (archived result): 37.72% of candidates flagged; PSM median = 1 in both flagged and unflagged classes.

Inputs: output/_base_data_chain/supplementary_table4_with_recomputed_props.csv (from 04_compute_peptide_properties.py;
      supplementary_table4_high_confidence.csv also works)
Output: output/18_revision_supplements/Fig_low_complexity_analysis.csv
      （Pep_ID2/top1_frac/entropy_bits/low_complexity + PSM）
Usage: py 05_flag_low_complexity_sequences.py --out ./output
"""
import pandas as pd
import numpy as np
import os
from collections import Counter
from math import log2
import os as _os
import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

OUT = OUTROOT + "/18_revision_supplements"
os.makedirs(OUT, exist_ok=True)

SRC = OUTROOT + "/_base_data_chain/supplementary_table4_with_recomputed_props.csv"
if not os.path.exists(SRC):
    SRC = OUTROOT + "/_base_data_chain/supplementary_table4_high_confidence.csv"

s4 = pd.read_csv(SRC, low_memory=False)
u = s4.drop_duplicates(subset=["Pep_ID2"] if "Pep_ID2" in s4.columns else ["Pep_ID"]).copy()
print("Unique candidates:", len(u), "(target 5,684)")

def stats_row(seq):
    seq = str(seq)
    c = Counter(seq)
    n = len(seq)
    top1 = max(c.values()) / n
    ent = -sum((v / n) * log2(v / n) for v in c.values())
    return top1, ent

u[["top1_frac", "entropy_bits"]] = u["Sequence"].apply(lambda s: pd.Series(stats_row(s)))
# Same convention as the original analysis: round to 3 decimals first, then threshold
u["top1_frac"] = u["top1_frac"].round(3)
u["entropy_bits"] = u["entropy_bits"].round(3)
u["low_complexity"] = (u["entropy_bits"] < 2.5) | (u["top1_frac"] >= 0.5)

cols = ["Pep_ID", "Sequence", "PSM", "top1_frac", "entropy_bits", "low_complexity"]
if "Pep_ID2" in u.columns:
    cols = ["Pep_ID2", "Sequence", "PSM", "top1_frac", "entropy_bits", "low_complexity"]
res = u[cols].copy()
res.to_csv(f"{OUT}/Fig_low_complexity_analysis.csv", index=False, encoding="utf-8-sig")

flag = res["low_complexity"]
print(f"Flagged low-complexity: {flag.sum()} ({100*flag.mean():.2f}%) [archive: 37.72%]")
for g, d in res.groupby(flag):
    print(f"  low_complexity={g}: n={len(d)}, PSM median = {d['PSM'].median()} (archive: 1 in both classes)")
print("saved:", OUT + "/Fig_low_complexity_analysis.csv")
