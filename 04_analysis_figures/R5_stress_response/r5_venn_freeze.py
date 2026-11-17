from pathlib import Path
import pandas as pd, numpy as np

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
A = f"{B}/output/04_article/R5_stress_response"
D = f"{A}/data"

TREAT_TRAITS = ["DT_48", "DT_144", "ST_48", "ST_144"]

def sig(panel, comp):
    d = pd.read_csv(f"{MS}/01_deseq2/{panel}/{comp}_results.tsv", sep="\t")
    d = d.rename(columns={d.columns[0]: "tid"})
    return set(d.loc[(d["padj"] < 0.05) & (d["log2FoldChange"].abs() >= 1), "tid"])

S1 = set()
for c in ["DT_48_vs_CK", "DT_144_vs_CK", "ST_48_vs_CK", "ST_144_vs_CK"]:
    S1 |= sig("ST_DT", c)
S3 = sig("heat", "heat_vs_heat_CK")

me_c = pd.read_csv(f"{MS}/02_wgcna/ST_DT/ME_trait_correlation.tsv", sep="\t", index_col=0)
me_p = pd.read_csv(f"{MS}/02_wgcna/ST_DT/ME_trait_padj.tsv", sep="\t", index_col=0)
treat_mods = [m for m in me_c.index
              if any(me_p.loc[m, t] < 0.05 and abs(me_c.loc[m, t]) >= 0.6
                     for t in TREAT_TRAITS)]
gm = pd.read_csv(f"{MS}/02_wgcna/ST_DT/gene_module_assignment.tsv", sep="\t")
core = gm[(gm["module"].isin(treat_mods)) & (gm["kME_own_module"] >= 0.7)]
S2 = set(core["transcript_id"])

n12, n13, n23 = len(S1 & S2), len(S1 & S3), len(S2 & S3)
triple = S1 & S2 & S3
print(f"Treatment-associated modules({len(treat_mods)}): {sorted(treat_mods)}")
print(f"S1 DT/ST DE={len(S1)} | S2 WGCNA core={len(S2)} | S3 Heat DE={len(S3)}")
print(f"S1∩S2={n12} S1∩S3={n13} S2∩S3={n23} triple={len(triple)}")
r = {"DT/ST only": len(S1 - S2 - S3), "WGCNA only": len(S2 - S1 - S3),
     "Heat only": len(S3 - S1 - S2), "DT/ST∩WGCNA only": len((S1 & S2) - S3),
     "DT/ST∩Heat only": len((S1 & S3) - S2), "WGCNA∩Heat only": len((S2 & S3) - S1),
     "triple": len(triple)}
print("Seven regions:", r)
pd.DataFrame([{"exclusive_region": k, "n": v} for k, v in r.items()]).to_csv(
    f"{D}/R5_venn3_exclusive_counts.tsv", sep="	", index=False)
print("Seven-region output:", f"{D}/R5_venn3_exclusive_counts.tsv")


rk = pd.read_csv(f"{D}/R5_candidate_ranking.tsv", sep="\t").set_index("transcript_id")
mod = gm.set_index("transcript_id")
rows = []
for tid in sorted(triple):
    k = rk.loc[tid]
    m = mod.loc[tid]
    rows.append({
        "transcript_id": tid, "locus": k.get("locus", ""),
        "wgcna_module": m["module"], "kME": m["kME_own_module"],
        "n_ncp_peptides": int(k["n_ncp_peptides"]), "psms_sum": int(k["psms_sum"]),
        "psm_tier": k.get("psm_tier", "strong" if k["psms_sum"] >= 3 else "low"),
        "pass_display_st": bool(k["pass_display_st"]), "pass_display_ht": bool(k["pass_display_ht"]),
        "LFC_DT48": k["DT_48_LFC"], "LFC_DT144": k["DT_144_LFC"],
        "LFC_ST48": k["ST_48_LFC"], "LFC_ST144": k["ST_144_LFC"],
        "LFC_heat": k["heat_LFC"],
        "sig_pattern_dtst": k.get("sig_pattern", ""),
    })
venn = pd.DataFrame(rows)

deg_cols = ["DT_48_DEG", "DT_144_DEG", "ST_48_DEG", "ST_144_DEG", "heat_DEG"]
lfc_cols = ["DT_48_LFC", "DT_144_LFC", "ST_48_LFC", "ST_144_LFC", "heat_LFC"]
opp = []
for tid in venn["transcript_id"]:
    k = rk.loc[tid]
    best, bl = None, 0
    for dc, lc in zip(deg_cols[:4], lfc_cols[:4]):
        if k[dc] and abs(k[lc]) > abs(bl):
            best, bl = k[lc], k[lc]
    hl = k["heat_LFC"] if k["heat_DEG"] else np.nan
    opp.append(bool(pd.notna(bl) and pd.notna(hl) and np.sign(bl) != np.sign(hl)))
venn["dtst_vs_heat_opposite"] = opp
venn.to_csv(f"{D}/R5_venn3_sets.tsv", sep="\t", index=False)
print(f"\n28 candidates with strong PSM (>=3): {(venn['psms_sum'] >= 3).sum()} | "
      f"low: {(venn['psms_sum'] < 3).sum()}")
print("Not passing both projects' display filters:", venn.loc[~(venn["pass_display_st"] & venn["pass_display_ht"]),
                                     "transcript_id"].tolist())
print("Opposing DT/ST and heat response directions:", int(venn["dtst_vs_heat_opposite"].sum()))
print("Wrote:", f"{D}/R5_venn3_sets.tsv")
