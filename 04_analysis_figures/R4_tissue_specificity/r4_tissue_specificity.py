from pathlib import Path
import pandas as pd, numpy as np
from scipy import stats

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article/R4_tissue_specificity"
D = f"{A}/data"
sd = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/sample_design.tsv", sep="\t")
vtpm = pd.read_csv(f"{B}/output/02filter_trans/verified_transcript_TPM.tsv", sep="\t")
nt = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/ncp_to_transcript.tsv", sep="\t")

tissue_sd = sd[sd["panel"] == "tissue"].copy()
order = ["leaf", "central peel", "the edges of peels", "Xylem", "fruit"]
conds = tissue_sd["condition"].unique().tolist()
print("tissue panel conditions:", conds)
sample_of = {c: tissue_sd.loc[tissue_sd["condition"] == c, "sample"].tolist() for c in conds}
for c, s in sample_of.items():
    assert len(s) == 3, (c, s)

focus = set(nt["transcript_id"])
sub = vtpm[vtpm["transcript_id"].isin(focus)].copy()
assert len(sub) == 1741, len(sub)
samples_all = [s for c in order for s in sample_of[c]]
T = sub.set_index("transcript_id")[samples_all].astype(float)


det = pd.DataFrame(index=T.index)
grp_mean = pd.DataFrame(index=T.index)
for c in order:
    ss = sample_of[c]
    det[c] = (T[ss] >= 1.0).sum(axis=1) >= 2
    grp_mean[c] = T[ss].mean(axis=1)

# ---- Tau ----
gm = grp_mean.copy()
mx = gm.max(axis=1)
tau = pd.Series(np.nan, index=gm.index)
has = mx > 0
tau[has] = ((1 - gm.loc[has].div(mx[has], axis=0)).sum(axis=1)) / (len(order) - 1)
det_n = det.sum(axis=1)
dominant = gm.idxmax(axis=1).where(has)

res = pd.DataFrame({
    "transcript_id": gm.index,
    "locus": sub.set_index("transcript_id")["locus"],
    "class_code": sub.set_index("transcript_id")["class_code"],
    "n_peptides": sub.set_index("transcript_id")["n_peptides"],
    "tau": tau.round(4), "n_tissues_detected": det_n,
    "dominant_tissue": dominant,
    "detect_leaf": det["leaf"].astype(int), "detect_central_peel": det["central peel"].astype(int),
    "detect_edges_of_peels": det["the edges of peels"].astype(int),
    "detect_Xylem": det["Xylem"].astype(int), "detect_fruit": det["fruit"].astype(int),
    "mean_TPM_leaf": grp_mean["leaf"].round(3), "mean_TPM_central_peel": grp_mean["central peel"].round(3),
    "mean_TPM_edges_of_peels": grp_mean["the edges of peels"].round(3),
    "mean_TPM_Xylem": grp_mean["Xylem"].round(3), "mean_TPM_fruit": grp_mean["fruit"].round(3),
})
res["tissue_specific_group"] = np.where(
    (res["tau"] >= 0.8) & (res["n_tissues_detected"] >= 1), res["dominant_tissue"] + "_specific",
    np.where(res["n_tissues_detected"] >= len(order), "widespread",
             np.where(res["n_tissues_detected"] == 0, "not_detected", "intermediate")))
res["overall_mean_TPM"] = T.mean(axis=1).round(3)
res.to_csv(f"{D}/tissue_specificity_results.tsv", sep="\t", index=False)


combo = det.apply(lambda r: "|".join(c for c in order if r[c]) or "none", axis=1)
upset = combo.value_counts()
upset_df = pd.DataFrame({"detected_combination": upset.index, "n_transcripts": upset.values})
upset_df.to_csv(f"{D}/fig4A_detection_upset.tsv", sep="\t", index=False)
per_tissue = pd.DataFrame({"tissue": order,
                           "n_detected": [int(det[c].sum()) for c in order],
                           "pct_of_1741": [round(100*det[c].mean(),1) for c in order]})
per_tissue.to_csv(f"{D}/fig4A_per_tissue_detection.tsv", sep="\t", index=False)
print("Detection by tissue:\n", per_tissue.to_string(index=False))


res.loc[res["n_tissues_detected"]>0].to_csv(f"{D}/fig4B_tau_values.tsv", sep="\t", index=False)
tau_stats = {"n_expressing": int(has.sum()), "tau_median": round(float(tau[has].median()),3),
             "tau_ge_0.8": int(((tau>=0.8)&has).sum()), "tau_ge_0.9": int(((tau>=0.9)&has).sum()),
             "widespread_5tissues": int((res["tissue_specific_group"]=="widespread").sum())}
dom_counts = res.loc[res["tissue_specific_group"].str.endswith("_specific"),"dominant_tissue"].value_counts()
pd.DataFrame({"dominant_tissue": dom_counts.index, "n_specific_candidates": dom_counts.values}
             ).to_csv(f"{D}/fig4C_specific_candidates_by_tissue.tsv", sep="\t", index=False)
print("Tau:", tau_stats)
print("Tissue-specific candidates:\n", dom_counts.to_string())


e = res[has & res["class_code"].notna()].copy()
groups = [g["tau"].dropna().values for _, g in e.groupby("class_code") if len(g) >= 10]
h, p = stats.kruskal(*groups)
cls_summary = e.groupby("class_code").agg(n=("tau","size"), median_tau=("tau","median"),
                                          n_tau08=("tau", lambda s: int((s>=0.8).sum())))
cls_summary["pct_tau08"] = (100*cls_summary["n_tau08"]/cls_summary["n"]).round(1)
cls_summary.to_csv(f"{D}/fig4E_classcode_tau.tsv", sep="\t")
print(f"KW (class vs tau, n>=10 classes): H={h:.2f} p={p:.2e}")
print(cls_summary.to_string())


res[["transcript_id","tau","overall_mean_TPM"]].to_csv(f"{D}/figS7_tau_vs_expression.tsv", sep="\t", index=False)

with open(f"{A}/stats/R4_summary.md", "w", encoding="utf-8") as fh:
    fh.write(f"""# R4 Key results
- Dataset: T-focus 1,741; 15 tissue-panel samples (3 each for leaf, central peel, edges of peels, Xylem, and fruit)
- Detection by tissue (at least 2 of 3 replicates with TPM >= 1): {per_tissue.to_dict('records')}
- Tau ({tau_stats['n_expressing']} expressed transcripts): median {tau_stats['tau_median']}; >= 0.8: {tau_stats['tau_ge_0.8']}; >= 0.9: {tau_stats['tau_ge_0.9']}; detected in all five tissues: {tau_stats['widespread_5tissues']}
- Tissue-specific candidates (Tau >= 0.8 and detected), by dominant tissue: {dom_counts.to_dict()}
- Class code vs Tau Kruskal-Wallis: H={h:.2f}, p={p:.2e}; median Tau for class j: {cls_summary.loc['j','median_tau'] if 'j' in cls_summary.index else 'NA'}; for class u: {cls_summary.loc['u','median_tau'] if 'u' in cls_summary.index else 'NA'}
- Interpretation: detection and Tau cutoffs are analytical definitions; Tau uses mean TPM for five tissue groups; all-zero transcripts are excluded
""")
print("DONE R4")
