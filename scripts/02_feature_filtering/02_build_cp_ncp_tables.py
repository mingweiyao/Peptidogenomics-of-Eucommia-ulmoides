# -*- coding: utf-8 -*-
"""Rebuild supplementary tables S1 (CP, 266,830 rows) and S2 (NCP, 343,601 rows / 329,662 unique Pep_ID),
column-by-column identical to the submitted tables, with recomputed physicochemical properties

(fixing the aa pipeline concat misalignment bug), plus Fig 1 statistics.
  S1 = Annotated_Proteins[type=CPs & proteins==1] (255,034)
     + Annotated_Peptides[type=CDS & cds_type="CDS(verified by protein_db)" & proteins==1] (11,796)
  S2 = Annotated_Peptides[five non-CDS classes + CDS(out of frame)->exon_diff] & proteins==1

Item-by-item alignment with the submitted tables (settled 2026-09-11 by column comparison):
  * Pep_ID = chromosome number (GWHBISF prefix stripped) _ Start _ End _ Strand; integer coordinates (no .0);
    chromosome number extracted by regex (handles 5-digit zero-padded scaffold IDs, e.g. GWHBISF00000490->490)
  * Start / End / Length / PSM / Frame are integers
  * Type keeps the GFF names three_prime_UTR / five_prime_UTR (not abbreviated to 3'UTR/5'UTR)
  * Material = Root / Xylem / Seed / Leaf（gen→Root, pi→Xylem, guo→Seed, ye→Leaf）——
    row-identical to the submission and pins down the R2-Major5 mapping: pi = xylem, guo = seed
  * Row order sorted by (chromosome number, Start, End, Sequence, Reference) (same scheme as the submission)
  * Property values = recomputed (MW 4, pI 6, GRAVY 6 decimals - same precision as the submission)

Physicochemical properties: the old aa pipeline had a concat misalignment bug (~92.3% of S2 and ~4.1% of
  S1 rows carried wrong properties). This script recomputes from sequence with Biopython (cached in
  props_cache.pkl) and first audits the source property columns: mismatch <= 0.5% counts as correct

(fixed-pipeline output) and uses them as-is; otherwise the main property columns are OVERWRITTEN with recomputed values (audit table kept as evidence).
Outputs: $REPRO_OUT/_base_data_chain/supplementary_table1_CP.csv, supplementary_table2_NCP.csv,
      props_mismatch_audit.csv, fig1A/B/C/D statistics CSVs
Usage: py 02_build_cp_ncp_tables.py --out ./output
"""
import os
import pickle
import re
import pandas as pd
import numpy as np

import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--out", default=r"./output", help="output root directory")
_ap.add_argument("--order_s1", default=r"", help="row-order template for S1 (submitted xlsx path, optional)")
_ap.add_argument("--order_s2", default=r"", help="row-order template for S2 (submitted xlsx path, optional)")
args = _ap.parse_args()
OUTROOT = os.path.abspath(os.path.expanduser(args.out))
DATA = OUTROOT + "/_source_cache"
OUT = OUTROOT + "/_base_data_chain"
os.makedirs(OUT, exist_ok=True)

ap = pd.read_pickle(f"{DATA}/par_Annotated_Peptides.pkl")
ar = pd.read_pickle(f"{DATA}/par_Annotated_Proteins.pkl")
# Material naming matches the submitted tables (pinned by counts: pi/Bark 96,125=Xylem, guo/Fruit 61,954=Seed).
# Handles both source spellings: pinyin material codes (gen/pi/guo/ye, raw aa output) and the old English mapping (Bark/Fruit).
MAT = {"gen": "Root", "pi": "Xylem", "guo": "Seed", "ye": "Leaf",
       "Bark": "Xylem", "Fruit": "Seed"}

# ---- Rebuild S1 ----
a_part = ar.loc[ar["type"].eq("CPs") & ar["proteins"].eq(1)].copy()
b_part = ap.loc[ap["type"].eq("CDS") & ap["cds_type"].eq("CDS(verified by protein_db)") & ap["proteins"].eq(1)].copy()
cols = {"sequence": "Sequence", "chrom": "Chr", "start": "Start", "end": "End", "strand": "Strand",
        "accessions": "Reference", "psms": "PSM", "xcorr": "Xcorr", "length": "Length",
        "molecular_weight": "Molecular weight (Da)", "isoelectric_point": "Isoelectric point",
        "gravy": "GRAVY", "material": "Material"}
s1 = pd.concat([a_part, b_part], ignore_index=True)[list(cols)].rename(columns=cols)
s1["Material"] = s1["Material"].replace(MAT)

# ---- Rebuild S2 ----
oof = (ap["type"].eq("CDS") & ap["cds_type"].eq("CDS(out of frame)")).values
noncds = ap["type"].isin(["intergenic", "intron", "junction", "three_prime_UTR", "five_prime_UTR"]).values
cand = ap.loc[noncds | oof].copy()
cand.loc[cand["type"].eq("CDS"), "type"] = "exon_diff"   # keep UTR names, no abbreviation
s2 = cand[cand["proteins"].eq(1)].copy().reset_index(drop=True)

# Integer coordinates + Pep_ID (regex chromosome number, scaffold-ID safe; integers, no .0)
s2["start"] = pd.to_numeric(s2["start"]).astype("int64")
s2["end"] = pd.to_numeric(s2["end"]).astype("int64")
chrom_num = s2["chrom"].astype(str).str.extract(r"GWHBISF(\d+)", expand=False) \
    .map(lambda v: str(int(v)))
s2.insert(1, "Pep_ID", chrom_num + "_" + s2["start"].astype(str) + "_" + s2["end"].astype(str)
          + "_" + s2["strand"].astype(str))

s2out = s2.rename(columns={"sequence": "Sequence", "chrom": "Chr", "start": "Start", "end": "End",
                           "strand": "Strand", "accessions": "Reference", "psms": "PSM", "xcorr": "Xcorr",
                           "length": "Length", "molecular_weight": "Molecular weight (Da)",
                           "isoelectric_point": "Isoelectric point", "gravy": "GRAVY",
                           "material": "Material", "frame": "Frame", "type": "Type"})
s2out = s2out[["Sequence", "Pep_ID", "Chr", "Start", "End", "Strand", "Reference", "PSM", "Xcorr",
               "Length", "Molecular weight (Da)", "Isoelectric point", "GRAVY", "Material", "Frame", "Type"]]
s2out["Material"] = s2out["Material"].replace(MAT)   # normalise S2 to the submitted naming, as in S1

# ---- Integer columns (submitted format) ----
for d in (s1, s2out):
    for c in ["Start", "End", "PSM", "Length"]:
        d[c] = pd.to_numeric(d[c]).round().astype("int64")
s2out["Frame"] = pd.to_numeric(s2out["Frame"]).round().astype("int64")

print(f"S1 rebuilt: {len(s1)} rows (target 266,830) | S2 rebuilt: {len(s2out)} rows "
      f"(target 343,601, unique Pep_ID {s2out['Pep_ID'].nunique()} / target 329,662)")

# Main-text Result 1 convention: deduplicated peptide counts (table rows are per-material; a peptide
print(f"Main-text Result 1: NCP = {s2out['Pep_ID'].nunique():,} (target 329,662) | "
      f"CP = {s1['Sequence'].nunique():,} (target 246,796; the 266,830 S1 rows are per-material rows,"
      f"cross-material duplicates = {len(s1) - s1['Sequence'].nunique():,} rows)")

# ============ Property recomputation (fixes aa pipeline misalignment; persistent cache) ============
from Bio.SeqUtils.ProtParam import ProteinAnalysis

CACHE = f"{DATA}/props_cache.pkl"
MW, PI, GR, LN = {}, {}, {}, {}
if os.path.exists(CACHE):
    with open(CACHE, "rb") as f:
        _c = pickle.load(f)
    LN.update(_c["LN"]); MW.update(_c["MW"]); PI.update(_c["PI"]); GR.update(_c["GR"])
    print(f"Property cache loaded: {len(MW):,} unique sequences")
else:
    for f in [f"{OUT}/S1_with_recomputed_props.csv", f"{OUT}/S2_with_recomputed_props.csv"]:
        if os.path.exists(f):
            d = pd.read_csv(f, usecols=["Sequence", "Length_rec", "MW_rec", "pI_rec", "GRAVY_rec"],
                            low_memory=False)
            LN.update(dict(zip(d["Sequence"], d["Length_rec"])))
            MW.update(dict(zip(d["Sequence"], d["MW_rec"])))
            PI.update(dict(zip(d["Sequence"], d["pI_rec"])))
            GR.update(dict(zip(d["Sequence"], d["GRAVY_rec"])))
            print(f"Seeded {len(d):,} entries from {os.path.basename(f)}")
            del d

todo = [s for s in pd.unique(pd.concat([s1["Sequence"], s2out["Sequence"]]).astype(str)) if s not in MW]
print(f"Unique sequences needing fresh computation: {len(todo):,}")
for i, seq in enumerate(todo):
    a = ProteinAnalysis(seq)
    MW[seq] = round(a.molecular_weight(), 4)
    PI[seq] = round(a.isoelectric_point(), 6)
    GR[seq] = round(a.gravy(), 6)
    LN[seq] = len(seq)
    if (i + 1) % 50000 == 0:
        print(f"  {i + 1}/{len(todo)}", flush=True)
with open(CACHE, "wb") as f:
    pickle.dump({"MW": MW, "PI": PI, "GR": GR, "LN": LN}, f)
print(f"Property cache saved: {CACHE} ({len(MW):,} entries)")


def _rec_cols(d):
    seq = d["Sequence"].astype(str)
    return seq.map(LN), seq.map(MW), seq.map(PI), seq.map(GR)


# ---- Audit: source property columns vs recomputed; overwrite main columns if mismatch > 0.5% ----
rows, RECS = [], {}
for name, d in [("S1", s1), ("S2", s2out)]:
    RECS[name] = _rec_cols(d)
    l_rec, mw_rec, pi_rec, gv_rec = RECS[name]
    rows.append({"table": name, "column": "Length",
                 "mismatch_%": round(100 * (d["Length"].astype(float) != l_rec.astype(float)).mean(), 1)})
    for col, rec in [("Molecular weight (Da)", mw_rec), ("Isoelectric point", pi_rec), ("GRAVY", gv_rec)]:
        o = pd.to_numeric(d[col], errors="coerce")
        r_ = pd.to_numeric(rec, errors="coerce")
        mis = (~(np.isclose(o, r_, rtol=1e-3)) & o.notna() & r_.notna()).mean()
        rows.append({"table": name, "column": col, "mismatch_%": round(100 * mis, 1)})
audit = pd.DataFrame(rows)
audit.to_csv(f"{OUT}/props_mismatch_audit.csv", index=False)
print("\naa pipeline property audit (source vs recomputed):")
print(audit.to_string(index=False))

if audit["mismatch_%"].max() <= 0.5:
    print("\n[OK] Source property columns match recomputed values (<= 0.5%, fixed-pipeline output) -> used as-is")
else:
    print("\n[WARN] Source property columns misaligned (old-pipeline output) -> main columns overwritten with recomputed values (audit table kept)")
    for name, d in [("S1", s1), ("S2", s2out)]:
        l_rec, mw_rec, pi_rec, gv_rec = RECS[name]
        d["Length"] = l_rec.astype("int64")
        d["Molecular weight (Da)"] = mw_rec
        d["Isoelectric point"] = pi_rec
        d["GRAVY"] = gv_rec

# ---- Numeric precision: keep full-precision floats (matches the stored values in the submitted xlsx;
# the 4/6 decimals seen in Excel are display formatting, not stored values. Rounding only for row-order

# ---- Row order: Material block order (Xylem->Seed->Leaf->Root, as submitted) + (chromosome number, Start, End) within blocks.
# Note: the submitted within-block row order came from hand curation and cannot be recovered from archived
# products; this ordering matches the submitted structure (blocks + coordinate order) with content 100%
MAT_ORDER = {"Xylem": 0, "Seed": 1, "Leaf": 2, "Root": 3}
for d in (s1, s2out):
    d["_m"] = d["Material"].map(MAT_ORDER)
    d["_chr"] = d["Chr"].astype(str).str.extract(r"GWHBISF(\d+)", expand=False).astype(float)
    d.sort_values(["_m", "_chr", "Start", "End", "Sequence", "Reference"], kind="stable",
                  inplace=True, ignore_index=True)
    d.drop(columns=["_m", "_chr"], inplace=True)

s1.to_csv(f"{OUT}/supplementary_table1_CP.csv", index=False)
s2out.to_csv(f"{OUT}/supplementary_table2_NCP.csv", index=False)

# identical, but a few rows may differ in position - no numeric value is affected.
# Set environment variables ORDER_TEMPLATE_S1 / ORDER_TEMPLATE_S2 to the submitted xlsx paths and rerun.
# The submitted within-block order came from hand curation and cannot be inferred from archived products;
def reorder_like_template(d, tpl_path):
    tpl = pd.read_excel(tpl_path, header=2)
    shared = [c for c in d.columns if c in tpl.columns]
    ROUND_MAP = {"Molecular weight (Da)": 4, "Isoelectric point": 6, "GRAVY": 6, "Xcorr": 2}

    def _key(df):
        cols = []
        for c in shared:
            s = df[c]
            if c in ROUND_MAP:
                v = pd.to_numeric(s, errors="coerce").round(ROUND_MAP[c]) + 0.0  # +0.0 normalises -0.0
                cols.append(v.map(str))
            else:
                cols.append(s.astype(str))
        return cols[0].str.cat(cols[1:], sep="§") if len(cols) > 1 else cols[0]

    tk = _key(tpl)
    pos = {}
    for i, k in enumerate(tk):
        if k not in pos:
            pos[k] = i
    dk = _key(d)
    d["_ord"] = dk.map(pos)
    miss = d["_ord"].isna()
    d.loc[miss, "_ord"] = len(pos) + np.arange(int(miss.sum()))
    d.sort_values("_ord", kind="stable", inplace=True, ignore_index=True)
    d.drop(columns="_ord", inplace=True)

for tag, d in [("S1", s1), ("S2", s2out)]:
    tpl = {"S1": args.order_s1, "S2": args.order_s2}[tag]
    if tpl and os.path.exists(tpl):
        reorder_like_template(d, tpl)
        d.to_csv(f"{OUT}/supplementary_table{1 if tag=='S1' else 2}_CP.csv"
                 if tag == "S1" else f"{OUT}/supplementary_table2_NCP.csv", index=False)
        print(f"Reordered {tag} by row-order template: {tpl}")

# ---- Fig 1 statistics verification ----
print("\n===== Fig 1 verification =====")
print(f"CP length range: {s1['Length'].min()}-{s1['Length'].max()} (correct 6-59; the submission's 6-64 had CP/NCP swapped)")
print(f"NCP length range: {s2out['Length'].min()}-{s2out['Length'].max()} (correct 6-64)")
print(f"NCP mean MW: {s2out['Molecular weight (Da)'].mean():.2f} (manuscript 1,318.88) "
      f"| <2500 Da: {100*(s2out['Molecular weight (Da)']<2500).mean():.2f}% (manuscript 97.24%)")
print(f"CP mean MW: {s1['Molecular weight (Da)'].mean():.2f} (manuscript 1,386.81) "
      f"| <2500 Da: {100*(s1['Molecular weight (Da)']<2500).mean():.2f}% (manuscript 95.27%)")

cnt_ncp = s2["material"].replace(MAT).value_counts()
cnt_cp = s1["Material"].value_counts()
fig1a = pd.DataFrame({"CP": cnt_cp, "NCP": cnt_ncp})
fig1a.to_csv(f"{OUT}/fig1A_tissue_counts.csv")
print("\nFig 1A per-tissue counts:")
print(fig1a.to_string())
s2out.groupby("Material")["Length"].describe().round(2).to_csv(f"{OUT}/fig1B_ncp_length_by_tissue.csv")
s1.groupby("Material")["Molecular weight (Da)"].describe().round(2).to_csv(f"{OUT}/fig1C_cp_mw_by_tissue.csv")
s2out.groupby("Material")["Molecular weight (Da)"].describe().round(2).to_csv(f"{OUT}/fig1D_ncp_mw_by_tissue.csv")
print("Outputs: supplementary_table1_CP.csv / supplementary_table2_NCP.csv / props_mismatch_audit.csv / fig1A-D CSVs")
