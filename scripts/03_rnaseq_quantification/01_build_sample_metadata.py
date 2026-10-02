# -*- coding: utf-8 -*-
"""Supplementary table S3: metadata for the 122 transcriptome samples.

Rule: take the group sheet of Total_rna_seq.xlsx (Sample/Tissues/Group; 122 samples, 41 groups)
      left-joined to Sheet1 (NCBI SRA metadata, 189 Run-level rows), keeping 12 fields.
Verification: 122 rows; row-identical to the archived supplementary_table3_metadata.csv (recons/reproduction/data).

Input: REJ/01location/rnaseq/Total_rna_seq.xlsx (sheets: Sheet1 + group)
Output: output/S3_sample_metadata/supplementary_table3_metadata.csv
Usage: py 01_build_sample_metadata.py --out ./output --rej ./data/rejected_20260204
"""
import pandas as pd
import os
import os as _os
import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal", help="upstream archive root")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

REJ = os.path.abspath(os.path.expanduser(args.rej))
OUT = OUTROOT + "/S3_sample_metadata"
os.makedirs(OUT, exist_ok=True)

P = os.path.join(args.rej, '01location/rnaseq/Total_rna_seq.xlsx')
grp = pd.read_excel(P, sheet_name="group")
s1 = pd.read_excel(P, sheet_name="Sheet1")
print(f"group: {len(grp)} samples / {grp['Group'].nunique()} groups | Sheet1: {len(s1)} Runs")

meta = grp.merge(
    s1[["Run", "Experiment", "BioProject", "LibraryStrategy", "LibrarySelection",
        "LibraryLayout", "Platform", "Model", "CenterName", "SampleName"]],
    left_on="Sample", right_on="Run", how="left")
meta = meta.rename(columns={"Sample": "Sample"}).drop(columns=["Run"])
cols = ["Sample", "Tissues", "Group", "Experiment", "BioProject", "LibraryStrategy",
        "LibrarySelection", "LibraryLayout", "Platform", "Model", "CenterName", "SampleName"]
meta = meta[cols]
print(f"S3: {len(meta)} rows x {len(cols)} fields, missing cells {int(meta.isna().sum().sum())}")
meta.to_csv(f"{OUT}/supplementary_table3_metadata.csv", index=False, encoding="utf-8-sig")

# Compare with the archived version
ref = "./data/response_letter_reference/supplementary_table3_metadata.csv"
if os.path.exists(ref):
    old = pd.read_csv(ref)
    same = (old[cols].astype(str).reset_index(drop=True) == meta[cols].astype(str).reset_index(drop=True)).all().all()
    print("Matches archived supplementary_table3_metadata.csv:", bool(same))
print("saved:", OUT + "/supplementary_table3_metadata.csv")
