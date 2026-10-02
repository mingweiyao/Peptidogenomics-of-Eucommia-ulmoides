# -*- coding: utf-8 -*-
"""Transcriptomic filtering chain (second half) -> supplementary table S4 (high-confidence NCPs).

223,870 ->CV(groups with CV<=0.5 >= 95%)-> 148,266 (manuscript 148,268; 4 tied IDs are historically unstable, +/-2)
        ->CPM(all replicates >= 1 in at least one (Tissues,Group) group)-> 5,920 (exact)
        ->intersect S2(proteins==1)-> 5,684 unique (exact) -> expanded per tissue to 5,896 rows = supplementary table S4 (exact)

Inputs: output/_base_data_chain/gff_223870.pkl (from 01_gff_cleanup/01_clean_gff_coordinates.py), supplementary_table2_NCP.csv (from 02_build_cp_ncp_tables.py),
      REJ/01location/rnaseq/03ouput/finally_groupwise_cv.csv, total_combined_matrix_cpm.csv (manuscript-convention matrix),
      REJ/01location/rnaseq/Total_rna_seq.xlsx (group sheet, 41 groups)
Output: output/_base_data_chain/supplementary_table4_high_confidence.csv
Usage: py 03_filter_high_confidence_ncps.py --out ./output --rej ./data/rejected_20260204   (reads CSVs of several hundred MB; needs GBs of RAM; ~2-5 min)
"""
import pandas as pd
import numpy as np
import os as _os
import argparse
import os

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal", help="upstream archive root (rejected_20260204)")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))

REJ = os.path.abspath(os.path.expanduser(args.rej))
OUT = OUTROOT + "/_base_data_chain"

df3 = pd.read_pickle(f"{OUT}/gff_223870.pkl")
cv = pd.read_csv(f"{REJ}/01location/rnaseq/03ouput/finally_groupwise_cv.csv", low_memory=False)
cpm = pd.read_csv(f"{REJ}/01location/rnaseq/03ouput/total_combined_matrix_cpm.csv", low_memory=False)


def mkid(chrom, start, end, strand):
    num = str(int(str(chrom).replace("GWHBISF", "")))
    return f"{num}_{start}_{end}_{strand}"


ids_legacy = set(mkid(r.chrom, r.start, r.end, r.strand) for r in df3.itertuples())
gene_ids = set(cv["GeneID"].astype(str))
inter = ids_legacy & gene_ids
print(f"Intersection of legacy-cleaned IDs with CV-matrix GeneIDs: {len(inter)} / 223,870")

# CV filter: groups with CV<=0.5 must be >= 95% (41 groups -> at least 39)
group_cols = [c for c in cv.columns if c != "GeneID"]
print("Number of groups:", len(group_cols))
cvv = cv[group_cols].apply(pd.to_numeric, errors="coerce")
n_pass = (cvv <= 0.5).sum(axis=1)
need = int(np.ceil(0.95 * len(group_cols)))
in_legacy = cv["GeneID"].astype(str).isin(ids_legacy)
mask_cv = (n_pass >= need) & in_legacy
cv_pass = cv.loc[mask_cv, "GeneID"]
print(f"CV filter: groups >= {need}: {mask_cv.sum()} of 223,870 pass (target 148,268)")

# CPM filter: all replicates CPM>=1 in at least one group
grp = pd.read_excel(f"{REJ}/01location/rnaseq/Total_rna_seq.xlsx", sheet_name="group")
sample_groups = grp.groupby("Group")["Sample"].apply(list).to_dict()
print("Number of groups (group sheet):", len(sample_groups))
mat = cpm.set_index("GeneID")
cv_set = set(cv_pass.astype(str))
sub = mat[mat.index.astype(str).isin(cv_set)]
print("Rows passing CV and present in matrix:", len(sub))
groups_arr = {g: [s for s in samples if s in sub.columns] for g, samples in sample_groups.items()}
data = sub.to_numpy()
cols_idx = {s: i for i, s in enumerate(sub.columns)}
pass_mask = np.zeros(len(sub), bool)
row_ids = sub.index.astype(str).to_numpy()
for g, samples in groups_arr.items():
    idxs = [cols_idx[s] for s in samples if s in cols_idx]
    if not idxs:
        continue
    ok = (data[:, idxs] >= 1).all(axis=1)
    pass_mask |= ok
print(f"After CPM filter: {pass_mask.sum()} (target 5,920)")
final_ids = set(row_ids[pass_mask])

# Intersect with S2 (proteins==1 NCPs) -> 5,684 unique
s2 = pd.read_csv(f"{OUT}/supplementary_table2_NCP.csv", low_memory=False)
s2["Start"] = s2["Start"].astype("int64")
s2["End"] = s2["End"].astype("int64")
chrom_num = s2["Chr"].astype(str).str.replace("GWHBISF", "", regex=False).map(lambda v: str(int(v)))
pid = chrom_num + "_" + s2["Start"].astype(str) + "_" + s2["End"].astype(str) + "_" + s2["Strand"].astype(str)
s2 = s2.assign(Pep_ID2=pid)
final_set = s2[s2["Pep_ID2"].isin(final_ids)]
print(f"Intersection with proteins==1 NCPs: unique {final_set['Pep_ID2'].nunique()} (target 5,684)")
final_set.to_csv(f"{OUT}/supplementary_table4_high_confidence.csv", index=False)
print("Wrote supplementary_table4_high_confidence.csv")
