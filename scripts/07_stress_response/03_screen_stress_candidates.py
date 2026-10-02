# -*- coding: utf-8 -*-
"""Main figure 5 / supplementary table S5 one-stop: nominal DE (51/35) -> candidate intersection (35/27) -> comparison with the submitted tables.

Inputs:
  --out root (results written to <out>/S5_stress_candidates/)
  --rej  archive root: pep_deseq.xlsx, veen.xlsx, candidates_scored.xlsx, id_with_candidates_*.xlsx
  --pcr  submission package: supplementary table 5.xlsx (comparison baseline)
Outputs:
  DE_nominal_DT.csv / DE_nominal_ST.csv   nominal DE lists (51/35, name-identical to veen)
  S5_dt_final.csv / S5_st_final.csv       candidates 35/27 (Pep_ID mapped, fully matching the submitted S5)
"""
import argparse
import os

import pandas as pd

_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--rej", default=r"./data/rejected 20260204 horticultural plant journal",
                 help="upstream archive root")
_ap.add_argument("--pcr", default=r"./data/pcr_20260205",
                 help="submission package (comparison baseline)")
args = _ap.parse_args()

REJ = os.path.abspath(os.path.expanduser(args.rej))
PCR = os.path.abspath(os.path.expanduser(args.pcr))
OUT = os.path.join(os.path.abspath(os.path.expanduser(args.out)), "S5_stress_candidates")
os.makedirs(OUT, exist_ok=True)

# ---- Step 1: nominal DE lists (P<0.05 and |log2FC|>=1, union of both time points) ----
P = os.path.join(REJ, "02figure", "figure6", "deseq", "pep_deseq.xlsx")
sheets = {"DT48": "DT_48-vs-CK.DESeq2", "DT144": "DT_144-vs-CK.DESeq2",
          "ST48": "ST_48h-vs-CK.DESeq2", "ST144": "ST_144-vs-CK.DESeq2"}
sets = {}
for k, s in sheets.items():
    d = pd.read_excel(P, sheet_name=s)
    de = d[(d["PValue"] < 0.05) & (d["log2FC"].abs() >= 1)]
    up = int((de["log2FC"] > 0).sum())
    dn = int((de["log2FC"] < 0).sum())
    sets[k] = set(de["name"].astype(str))
    print(f"{k}: nominal DE {len(de)} (up {up} / down {dn})")

dt = sets["DT48"] | sets["DT144"]
st = sets["ST48"] | sets["ST144"]
print(f"DE NCP lists: DT {len(dt)} (target 51) | ST {len(st)} (target 35)")

veen = pd.read_excel(REJ + "/02figure/figure6/veen.xlsx")
vd = set(veen["dt_deseq"].dropna().astype(str))
vs = set(veen["st_deseq"].dropna().astype(str))
print(f"Name-identical to veen.xlsx dt_deseq: {dt == vd} | st_deseq: {st == vs}")

pd.DataFrame({"EuNCP_trans": sorted(dt)}).to_csv(f"{OUT}/DE_nominal_DT.csv", index=False)
pd.DataFrame({"EuNCP_trans": sorted(st)}).to_csv(f"{OUT}/DE_nominal_ST.csv", index=False)

# ---- Step 2: candidate intersection + Pep_ID mapping ----
dt_mod = set(veen["dt_module"].dropna().astype(str))
st_mod = set(veen["st_module"].dropna().astype(str))
print(f"Module high-TOM NCPs (archive): dt_module {len(dt_mod)} / st_module {len(st_mod)} (main text 231/181)")

dt_cand = dt & dt_mod
st_cand = st & st_mod
print(f"DE intersect modules: DT {len(dt_cand)} (target 35) | ST {len(st_cand)} (target 27)")

cs = pd.read_excel(REJ + "/02figure/figure5/candidates_scored.xlsx")
mp = dict(zip(cs["ID"].astype(str), cs["pep_id"].astype(str)))
for name, cand in [("dt", dt_cand), ("st", st_cand)]:
    pd.DataFrame({"EuNCP_trans": sorted(cand),
                  "Pep_ID": [mp.get(x, x) for x in sorted(cand)]}).to_csv(
        f"{OUT}/S5_{name}_final.csv", index=False, encoding="utf-8-sig")

# ---- Step 3: comparison with the submitted S5 ----
s5 = pd.read_excel(PCR + "/supplementary table 5.xlsx", header=2)
sub_dt = set(s5.loc[s5["Condition"].astype(str).str.contains("rought|DT", na=False), "Pep_ID"].astype(str))
sub_st = set(s5.loc[s5["Condition"].astype(str).str.contains("alt|ST", na=False), "Pep_ID"].astype(str))
mine_dt = {mp.get(x, x) for x in dt_cand}
mine_st = {mp.get(x, x) for x in st_cand}
print(f"Intersection with submitted S5 Drought {len(mine_dt & sub_dt)}/35 | Salt {len(mine_st & sub_st)}/27")
print("DONE")
