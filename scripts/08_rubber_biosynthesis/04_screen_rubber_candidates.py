# -*- coding: utf-8 -*-
"""Main figure 6 / supplementary table S6 one-stop: Pearson verification (301/85/32) -> partial correlation (82/47) -> S6A/B table split.

Inputs:
  --out root (results written to <out>/S6_rubber_candidates/)
  --rej  archive root: rubber_count_5_fpkm.txt, id.xlsx, candidates_scored.xlsx
  --pcr  submission package: supplementary table 6.xlsx (comparison baseline)
Outputs:
  S6_pairs_matrix_verified.csv       301 pairs r/p reproduced bit-exactly (anchor r=0.704027)
  id301_partial_recomputed.csv       partial correlation (tissue dummy-variable residual method; 82 pairs / 47 pairs)
  Table_S6B_pairs_301.csv / Table_S6A_unique_NCP_85.csv   S6 split tables (BH q<0.05)
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal",
                 help="upstream archive root")
_ap.add_argument("--pcr", default=r"./data/pcr_20260205",
                 help="submission package (comparison baseline)")
args = _ap.parse_args()

REJ = os.path.abspath(os.path.expanduser(args.rej))
PCR = os.path.abspath(os.path.expanduser(args.pcr))
OUT = os.path.join(os.path.abspath(os.path.expanduser(args.out)), "S6_rubber_candidates")
os.makedirs(OUT, exist_ok=True)
BASE = REJ + "/02figure/figure5"


def bh_qvalues(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# ---- Step 1: bit-exact reproduction of the 301 Pearson pairs ----
M = pd.read_csv(BASE + "/rubber/WGCNA/rubber_count_5_fpkm.txt", sep="\t", index_col=0)
idp = pd.read_excel(BASE + "/id.xlsx", sheet_name="Sheet1")
print("Matrix:", M.shape, "| reference pairs:", idp.shape)

a, b = "EuNCP_trans_4462", "evm.TU.Chr17.1382"
r0, p0 = stats.pearsonr(M.loc[a].astype(float), M.loc[b].astype(float))
print(f"Anchor {a} vs {b}: r={r0:.6f} p={p0:.6e} (expect 0.704027 / 3.394415e-03)")

rs, ps, miss = [], [], 0
for _, row in idp.iterrows():
    x, y = str(row["novel"]), str(row["annotated"])
    if x not in M.index or y not in M.index:
        miss += 1
        rs.append(np.nan)
        ps.append(np.nan)
        continue
    r, p = stats.pearsonr(M.loc[x].astype(float), M.loc[y].astype(float))
    rs.append(r)
    ps.append(p)
idp["r_recomputed"] = rs
idp["p_recomputed"] = ps
ok_r = int((idp["r_recomputed"] - idp["pearson_r"]).abs().lt(5e-7).sum())
print(f"301 pairs: missing rows {miss} | r exactly equal {ok_r}/301 | max |dr| = "
      f"{(idp['r_recomputed'] - idp['pearson_r']).abs().max():.3e}")
print(f"unique novel = {idp['novel'].nunique()} (target 85), "
      f"unique annotated = {idp['annotated'].nunique()} (target 32)")
idp.to_csv(f"{OUT}/S6_pairs_matrix_verified.csv", index=False, encoding="utf-8-sig")

# ---- Step 2: partial correlation (controlling 15x5 tissue dummies) ----
Z = np.zeros((15, 5))
for j in range(5):
    Z[j * 3:(j + 1) * 3, j] = 1.0


def resid(v):
    beta, _, _, _ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ beta


rows = []
for _, pair in idp.iterrows():
    x = M.loc[str(pair["novel"])].astype(float).values
    y = M.loc[str(pair["annotated"])].astype(float).values
    r_raw, p_raw = stats.pearsonr(x, y)
    rp = np.corrcoef(resid(x), resid(y))[0, 1]
    t = rp * np.sqrt(13 / (1 - rp * rp))
    p_part = 2 * stats.t.sf(abs(t), 13)
    rows.append({"ncp": pair["novel"], "gene": pair["annotated"],
                 "r_raw": r_raw, "p_raw": p_raw, "r_partial": rp, "p_partial": p_part})
res = pd.DataFrame(rows)
res["q_raw"] = bh_qvalues(res["p_raw"])
res["q_partial"] = bh_qvalues(res["p_partial"])
res.to_csv(f"{OUT}/id301_partial_recomputed.csv", index=False, encoding="utf-8-sig")
print(f"Partial correlation: p_partial<0.05 = {int((res['p_partial']<0.05).sum())} pairs / "
      f"{res.loc[res['p_partial']<0.05,'ncp'].nunique()} NCPs (target 82/44)")
print(f"         q_partial<0.05 = {int((res['q_partial']<0.05).sum())} pairs / "
      f"{res.loc[res['q_partial']<0.05,'ncp'].nunique()} NCPs (target 47/32)")

# ---- Step 3: S6A/B table split (compared with the submission) ----
cs = pd.read_excel(REJ + "/02figure/figure5/candidates_scored.xlsx")
mp = dict(zip(cs["ID"].astype(str), cs["pep_id"].astype(str)))
idp["Pep_ID"] = idp["novel"].astype(str).map(lambda x: mp.get(x, x))
s6 = pd.read_excel(PCR + "/supplementary table 6.xlsx", header=2)
ann = s6.drop_duplicates(subset=["Pep_ID"]).set_index("Pep_ID")
idp["q"] = bh_qvalues(idp["pearson_p"])
print(f"After BH q<0.05: {int((idp['q'] < 0.05).sum())}/301 (max q = {idp['q'].max():.4f})")

s6b = pd.DataFrame({"Pep_ID": idp["Pep_ID"], "Gene": idp["annotated"].astype(str),
                    "R": idp["pearson_r"], "P": idp["pearson_p"], "Q_BH": idp["q"]})
s6b["Pep_trans"] = s6b["Pep_ID"].map(ann["Pep_trans"] if "Pep_trans" in ann else pd.Series(dtype=object))
s6b["Exons"] = s6b["Pep_ID"].map(ann["Exons"] if "Exons" in ann else pd.Series(dtype=object))
s6b.insert(0, "Pair_ID", [f"S6B-{i+1}" for i in range(len(s6b))])
s6b.to_csv(f"{OUT}/Table_S6B_pairs_301.csv", index=False, encoding="utf-8-sig")

s6a = (s6b.sort_values("R", ascending=False).drop_duplicates(subset=["Pep_ID"]).reset_index(drop=True))
s6a.insert(0, "NCP_no", range(1, len(s6a) + 1))
s6a.to_csv(f"{OUT}/Table_S6A_unique_NCP_85.csv", index=False, encoding="utf-8-sig")

sub_keys = set(zip(s6["Pep_ID"].astype(str), s6["Rubber genes"].astype(str)))
my_keys = set(zip(s6b["Pep_ID"].astype(str), s6b["Gene"].astype(str)))
print(f"S6B {len(s6b)} pairs vs submitted 301 | key diff: ours only {len(my_keys - sub_keys)}, submitted only {len(sub_keys - my_keys)}")
print(f"S6A unique NCPs: {len(s6a)} (target 85)")
print("DONE")
