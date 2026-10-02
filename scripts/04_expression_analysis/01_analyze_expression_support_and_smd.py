# -*- coding: utf-8 -*-
"""Main figure 4 one-stop: every spatiotemporal-expression number (panels A/B/C/D/E).

Inputs:
  --out root: _base_data_chain/supplementary_table4_high_confidence.csv (from 03_filter_high_confidence_ncps.py)
  --rej upstream archive root (rejected_20260204):
         01location/rnaseq/03ouput/finally_expressed_sp_cpm.csv (universe matrix, panels B/C)
         01location/rnaseq/Total_rna_seq.xlsx (sample grouping)
         02figure/figure4/rnaseq_veen_analysis.xlsx (original Venn workbook for panel B, verification)
         02figure/figure4/SMD_leaf_stage.xlsx / SMD_leaf_located.xlsx (pipeline outputs for panels D/E)
Outputs (<--out>/Fig4_expression_smd/):
  fig4A_ms_tissue_counts.csv   Panel A: MS counts per tissue (2,958/1,647/1,012/279)
  fig4B_veen_reconciliation.csv Panel B: 15 Venn regions reconciled (all OK, intersection 4,960)
  fig4C_detection_distribution.csv Panel C: detection breadth across all 122 samples (2,880)
  fig4DE_SMD_recomputed.csv    Panels D/E: SMD recomputed and cross-checked
On-screen verification: A four tissue counts, B 15 regions all OK, C 2,880, D/E four numbers (32.04/69.14/32.17/89.74)
"""
import argparse
import os

import numpy as np
import pandas as pd

_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal",
                 help="upstream archive root")
args = _ap.parse_args()

REJ = os.path.abspath(os.path.expanduser(args.rej))
OUT = os.path.join(os.path.abspath(os.path.expanduser(args.out)), "Fig4_expression_smd")
os.makedirs(OUT, exist_ok=True)
FIG4 = os.path.join(REJ, "02figure", "figure4")

# ---- Panel A: MS counts per tissue (S4 tissue-level rows, not deduplicated) ----
s4 = pd.read_csv(os.path.join(os.path.dirname(OUT), "_base_data_chain",
                              "supplementary_table4_high_confidence.csv"), low_memory=False)
cnt = s4["Material"].value_counts()
print("A counts per material (tissue level, total should be 5,896):")
print(cnt.to_string(), "| total", int(cnt.sum()))
print("union unique Pep_ID:", s4["Pep_ID2"].nunique(), "(target 5,684)")
cnt.to_csv(os.path.join(OUT, "fig4A_ms_tissue_counts.csv"), encoding="utf-8-sig")

# ---- Panel B: transcript-support Venn (settled manuscript-convention verification) ----
ROOT4 = os.path.join(REJ, "01location")
mat = pd.read_csv(os.path.join(ROOT4, "rnaseq", "03ouput", "finally_expressed_sp_cpm.csv"), index_col=0)
grp = pd.read_excel(os.path.join(ROOT4, "rnaseq", "Total_rna_seq.xlsx"), sheet_name="group")
ana = pd.read_excel(FIG4 + "/rnaseq_veen_analysis.xlsx")
FIX = {"Xylem": "Bark", "fruit": "Fruit", "leaf": "Leaf", "root": "Root"}
declared = {}
for c in ana.columns[1:]:
    v = ana[c].iloc[0]
    if pd.notna(v) and str(v) != "value":
        declared[" ".join(FIX.get(w, w) for w in c.split())] = int(v)
groups = {}
for t in ["root", "Xylem", "fruit", "leaf", "peel"]:
    groups[t] = grp[grp["Tissues"] == t].groupby("Group")["Sample"].apply(list).to_dict()
groups["Bark"] = groups.pop("Xylem")
groups["Leaf"] = groups.pop("leaf")
groups["Fruit"] = {**groups.pop("fruit"), **groups.pop("peel")}
groups["Root"] = groups.pop("root")   # normalise region name: root -> Root (matches the declared naming)
S = {}
for t, gg in groups.items():
    X = mat.loc[:, sum(gg.values(), [])]
    ok = np.zeros(len(X), bool)
    for g, ss in gg.items():
        ok |= (X[ss] >= 1).all(axis=1).values
    S[t] = set(X.index[ok])
    print(f"B {t:6s}: {sum(len(v) for v in gg.values()):3d} samples {len(gg):2d} groups -> expressed {len(S[t]):5d}")
ts = ["Bark", "Fruit", "Leaf", "Root"]
rows, allok = [], True
for name, dd in declared.items():
    inside = name.split()
    reg = set.intersection(*[S[t] for t in inside])
    for t in [t for t in ts if t not in inside]:
        reg -= S[t]
    okk = len(reg) == dd
    allok &= okk
    rows.append({"region": name, "declared": dd, "recomputed": len(reg), "match": okk})
print(f"B all 15 regions consistent? {allok} (intersection = "
      f"{len(set.intersection(*[S[t] for t in ts]))}, target 4,960)")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "fig4B_veen_reconciliation.csv"),
                          index=False, encoding="utf-8-sig")

# ---- Panel C: detection breadth (universe matrix) ----
ndet = (mat >= 1).sum(axis=1)
print(f"C detected in all 122 samples: {int((ndet == 122).sum())} (manuscript 2,880)")
pd.Series(ndet).value_counts().sort_index(ascending=False).to_csv(
    os.path.join(OUT, "fig4C_detection_distribution.csv"))

# ---- Panels D/E: SMD (pipeline output files as reference + independent recomputation) ----
for tag, fn, target in [
    ("developmental stage", "SMD_leaf_stage.xlsx", "MAX>10 1,657 (32.04%) | MIN<1 3,575 (69.14%)"),
    ("spatial position", "SMD_leaf_located.xlsx", "MAX>10 1,627 (32.17%) | MIN<1 4,539 (89.74%)"),
]:
    dd = pd.read_excel(os.path.join(FIG4, fn))
    mx, mn = dd["MAX"], dd["MIN"]
    print(f"D/E {tag}: {len(dd)} rows | MAX>10: {int((mx>10).sum())} ({100*(mx>10).mean():.2f}%) | "
          f"MIN<1: {int((mn<1).sum())} ({100*(mn<1).mean():.2f}%) | target: {target}")

# Independent recomputation cross-check (SMD_analysis.xlsx, a different normalisation convention)
SMD = os.path.join(FIG4, "SMD_analysis.xlsx")
id_df = pd.read_excel(SMD, sheet_name="ID")
group_df = pd.read_excel(SMD, sheet_name="Group")
q = pd.read_excel(SMD, sheet_name="CPM")
q["ID"] = q["ID"].astype(str)
fd = q[q["ID"].isin(id_df["ID"].astype(str))]
gmap = group_df.groupby("Group")["Sample"].apply(list).to_dict()
gs = {}
for gname, ss in gmap.items():
    x = fd[[c for c in ss if c in fd.columns]]
    gs[gname] = (x.mean(axis=1), x.std(axis=1))
smd = pd.DataFrame({"ID": fd["ID"].astype(str).values})
done = set()
for g1 in gs:
    for g2 in gs:
        if g2 in done or g1 == g2:
            continue
        pooled = np.sqrt((gs[g1][1] ** 2 + gs[g2][1] ** 2) / 2)
        smd[f"SMD_{g1}_vs_{g2}"] = np.abs(gs[g1][0] - gs[g2][0]) / pooled
    done.add(g1)
smdcols = [c for c in smd.columns if c.startswith("SMD_")]
mx, mn = smd[smdcols].max(axis=1), smd[smdcols].min(axis=1)
print(f"[cross-check convention] max SMD>10: {(mx>10).sum()} ({100*(mx>10).mean():.2f}%) | "
      f"min SMD<1: {(mn<1).sum()} ({100*(mn<1).mean():.2f}%)")
smd.to_csv(os.path.join(OUT, "fig4DE_SMD_recomputed.csv"), index=False)

print("\nOutputs ->", OUT)
