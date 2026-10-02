# -*- coding: utf-8 -*-
"""Append recomputed physicochemical property columns to S4 (high-confidence NCPs).

Note: 02_build_cp_ncp_tables.py already writes *_rec recomputed columns for S1/S2 in its own output
(and persists the property cache props_cache.pkl). This script only handles S4: values come from the
cache (S4 sequences are a subset of S2, so nearly all hit; missing ones are computed fresh with Biopython),

producing supplementary_table4_with_recomputed_props.csv.
  S4 unique：>15aa 1,159 (20.39%)、pI<7 5,019 (88.30%)、GRAVY<0 3,855（Biopython；
  the main text uses the author-confirmed 3,853 (67.79%) - do not treat as a correction), MW mean 1,233.14 / <2500 98.75%.

Inputs: $REPRO_OUT/_base_data_chain/supplementary_table4_high_confidence.csv (from 03_filter_high_confidence_ncps.py),
      $REPRO_OUT/_source_cache/props_cache.pkl (from 02_build_cp_ncp_tables.py; computed fresh if absent)
Output: $REPRO_OUT/_base_data_chain/supplementary_table4_with_recomputed_props.csv
Usage: py 04_compute_peptide_properties.py --out ./output   (~1-2 min with cache)
"""
import os
import pickle
import pandas as pd
import numpy as np

import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))
DATA = OUTROOT + "/_source_cache"
OUT = OUTROOT + "/_base_data_chain"

s4 = pd.read_csv(f"{OUT}/supplementary_table4_high_confidence.csv", low_memory=False)
s4["Sequence"] = s4["Sequence"].astype(str)

MW, PI, GR, LN = {}, {}, {}, {}
if os.path.exists(f"{DATA}/props_cache.pkl"):
    with open(f"{DATA}/props_cache.pkl", "rb") as f:
        _c = pickle.load(f)
    LN.update(_c["LN"]); MW.update(_c["MW"]); PI.update(_c["PI"]); GR.update(_c["GR"])
    print(f"Property cache loaded: {len(MW):,} entries")
else:
    print("[WARN] props_cache.pkl not found (product of 02_build_cp_ncp_tables.py); computing S4 sequences fresh (a few minutes)")

todo = [s for s in pd.unique(s4["Sequence"]) if s not in MW]
if todo:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    print(f"Computing missing sequences fresh: {len(todo):,}")
    for i, seq in enumerate(todo):
        a = ProteinAnalysis(seq)
        MW[seq] = round(a.molecular_weight(), 4)
        PI[seq] = round(a.isoelectric_point(), 6)
        GR[seq] = round(a.gravy(), 6)
        LN[seq] = len(seq)
    if os.path.exists(f"{DATA}/props_cache.pkl"):
        with open(f"{DATA}/props_cache.pkl", "rb") as f:
            _c = pickle.load(f)
        LN.update(_c["LN"]); MW.update(_c["MW"]); PI.update(_c["PI"]); GR.update(_c["GR"])
    with open(f"{DATA}/props_cache.pkl", "wb") as f:
        pickle.dump({"MW": MW, "PI": PI, "GR": GR, "LN": LN}, f)

s4["Length_rec"] = s4["Sequence"].map(LN)
s4["MW_rec"] = s4["Sequence"].map(MW)
s4["pI_rec"] = s4["Sequence"].map(PI)
s4["GRAVY_rec"] = s4["Sequence"].map(GR)
s4.to_csv(f"{OUT}/supplementary_table4_with_recomputed_props.csv", index=False)

u4 = s4.drop_duplicates(subset=["Pep_ID"])
gt15 = (u4["Length_rec"] > 15).sum()
pi7 = (u4["pI_rec"] < 7).sum()
gv = (u4["GRAVY_rec"] < 0).sum()
print(f"\nS4 high-confidence (recomputed properties): >15aa = {gt15} ({100*gt15/len(u4):.2f}%) [manuscript 1,159 (20.39%)]")
print(f"pI<7 = {pi7} ({100*pi7/len(u4):.2f}%) [manuscript 5,019 (88.30%)]")
print(f"GRAVY<0 = {gv} ({100*gv/len(u4):.2f}%) [Biopython; main text uses the author-confirmed 3,853 (67.79%)]")
print(f"MW mean = {u4['MW_rec'].mean():.2f}, <2500 = {100*(u4['MW_rec']<2500).mean():.2f}%")
print("DONE → supplementary_table4_with_recomputed_props.csv")
