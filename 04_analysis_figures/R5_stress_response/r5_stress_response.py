from pathlib import Path
import pandas as pd, numpy as np, os
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article/R5_stress_response"
D = f"{A}/data"
DE = f"{B}/output/03_ms_trans_analysis/01_deseq2"
WG = f"{B}/output/03_ms_trans_analysis/02_wgcna"

def load_results(panel, comp):
    df = pd.read_csv(f"{DE}/{panel}/{comp}_results.tsv", sep="\t")
    idc = df.columns[0]
    df = df.rename(columns={idc: "transcript_id"})
    df["comparison"] = comp
    return df

def deg(df, padj=0.05, lfc=1.0):
    return df[(df["padj"] < padj) & (df["log2FoldChange"].abs() >= lfc)]


summ = pd.read_csv(f"{DE}/deseq2_summary.tsv", sep="\t").set_index("comparison")
checks = {}
for panel, comp in [("ST_DT","DT_48_vs_CK"),("ST_DT","DT_144_vs_CK"),
                    ("ST_DT","ST_48_vs_CK"),("ST_DT","ST_144_vs_CK"),
                    ("heat","heat_vs_heat_CK"),("MT_melatonin","drought_vs_CK"),
                    ("MT_melatonin","Melatonin_vs_CK"),("MT_melatonin","Melatonin_vs_drought")]:
    df = load_results(panel, comp)
    n = len(deg(df))
    target = int(summ.loc[comp, "n_sig"])
    checks[comp] = (n, target, n == target)
print("Threshold reconciliation (padj < 0.05 and |LFC| >= 1):", checks)
assert all(v[2] for v in checks.values()), "Threshold does not match v4 n_sig"


res = {c: load_results("ST_DT", c) for c in ["DT_48_vs_CK","DT_144_vs_CK","ST_48_vs_CK","ST_144_vs_CK"]}
degs = {c: set(deg(df)["transcript_id"]) for c, df in res.items()}
DT = degs["DT_48_vs_CK"] | degs["DT_144_vs_CK"]
ST = degs["ST_48_vs_CK"] | degs["ST_144_vs_CK"]
common, dt_only, st_only = DT & ST, DT - ST, ST - DT
union = DT | ST
counts = pd.DataFrame([
    {"set":"DT_union","n":len(DT)}, {"set":"ST_union","n":len(ST)},
    {"set":"common(DT∩ST)","n":len(common)}, {"set":"DT_specific","n":len(dt_only)},
    {"set":"ST_specific","n":len(st_only)}, {"set":"union","n":len(union)},
    {"set":"heat_DEG","n":len(deg(load_results("heat","heat_vs_heat_CK")))},
])
counts.to_csv(f"{D}/fig5A_DE_counts.tsv", sep="\t", index=False)
print(counts.to_string(index=False))


lfc_tab = pd.DataFrame({c: res[c].set_index("transcript_id")["log2FoldChange"] for c in res})
shared = lfc_tab.loc[lfc_tab.index.isin(union)]
sign = np.sign(shared)
sig = pd.DataFrame({c: shared.index.isin(degs[c]) for c in res}, index=shared.index)
pat = sign.where(sig, other=np.nan)
grp = pat.apply(lambda r: "|".join({("up" if v > 0 else "down") for v in r.dropna()}), axis=1)
consistent = pat.dropna(how="any").apply(lambda r: len(set(np.sign(r.dropna()))) == 1, axis=1)
dir_df = pd.DataFrame({"transcript_id": shared.index}).set_index("transcript_id")
dir_df["pattern"] = grp
dir_df["in_n_sig"] = sig.sum(axis=1)
dir_df["direction_consistent_across_4"] = consistent.reindex(dir_df.index).fillna(False)
dir_df = dir_df.join(lfc_tab.round(3))
dir_df.reset_index().to_csv(f"{D}/fig5B_direction_matrix.tsv", sep="\t", index=False)
print(f"Shared candidates with consistent direction across four contrasts: {int(consistent.sum())}")


vst = pd.read_csv(f"{WG}/input_vst/ST_DT_vst.tsv", sep="\t").set_index("transcript_id")
sd = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/sample_design.tsv", sep="\t")
st_samples = sd.loc[sd["panel"]=="ST_DT","sample"].tolist()
cond_of = sd.set_index("sample").loc[st_samples,"condition"].to_dict()
X = vst.loc[vst.index.isin(union), st_samples]
Z = (X.sub(X.mean(axis=1), axis=0)).div(X.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)
Ln = linkage(pdist(Z, "euclidean"), method="ward")
best_k, best_score = None, -1
for k in range(2, 9):
    lab = fcluster(Ln, k, criterion="maxclust")
    ss = 0.0
    for c in np.unique(lab):
        m = Z.values[lab == c]
        ss += ((m - m.mean(axis=0))**2).sum()
    within = ss / len(Z)
    if best_k is None or within < best_score * 0.92:  
        best_k, best_score = k, within
    else:
        break
lab = fcluster(Ln, best_k, criterion="maxclust")
clusters = pd.DataFrame({"transcript_id": Z.index, "cluster": lab})
prof = Z.copy(); prof["cluster"] = lab

long = prof.melt(id_vars="cluster", var_name="sample", value_name="z")
long["condition"] = long["sample"].map(cond_of)
cond_order = ["CK","DT_48","DT_144","ST_48","ST_144"]
cond_order = [c for c in cond_order if c in set(long["condition"])]
mean_prof = long.pivot_table(index="cluster", columns="condition", values="z", aggfunc="mean")
mean_prof["n"] = prof.groupby("cluster").size()
mean_prof.index = [f"C{c}" for c in mean_prof.index]

def name_row(r):
    ck = r["CK"]; d48, d144 = r.get("DT_48", ck), r.get("DT_144", ck)
    s48, s144 = r.get("ST_48", ck), r.get("ST_144", ck)
    dt_trend = "up" if d144 > ck and d48 > ck else ("down" if d144 < ck and d48 < ck else "mixed")
    st_trend = "up" if s144 > ck and s48 > ck else ("down" if s144 < ck and s48 < ck else "mixed")
    if dt_trend == st_trend == "up": return "co-upregulated"
    if dt_trend == st_trend == "down": return "co-downregulated"
    if dt_trend == "up" and st_trend == "mixed": return "drought_up"
    if dt_trend == "down" and st_trend == "mixed": return "drought_down"
    if st_trend == "up" and dt_trend == "mixed": return "salt_up"
    if st_trend == "down" and dt_trend == "mixed": return "salt_down"
    return "time_or_opposite_specific"
mean_prof["pattern"] = [name_row(mean_prof.loc[i]) for i in mean_prof.index]
mean_prof.to_csv(f"{D}/fig5C_cluster_mean_profiles.tsv", sep="\t")
cmap = dict(zip(mean_prof.index, mean_prof["pattern"]))
clusters["cluster_name"] = clusters["cluster"].map(lambda c: f"C{c}")
clusters["pattern"] = clusters["cluster"].map(lambda c: cmap[f"C{c}"])
print("Cluster k =", best_k, "\n", mean_prof[["n","pattern"]].to_string())





d_ck = load_results("MT_melatonin","drought_vs_CK").set_index("transcript_id")
m_ck = load_results("MT_melatonin","Melatonin_vs_CK").set_index("transcript_id")
m_d = load_results("MT_melatonin","Melatonin_vs_drought").set_index("transcript_id")
mt_sig = deg(m_d.reset_index())
rows = []
for tid in mt_sig["transcript_id"]:
    l_dck = d_ck.loc[tid, "log2FoldChange"] if tid in d_ck.index else np.nan
    p_dck = d_ck.loc[tid, "padj"] if tid in d_ck.index else 1.0
    l_mck = m_ck.loc[tid, "log2FoldChange"] if tid in m_ck.index else np.nan
    p_mck = m_ck.loc[tid, "padj"] if tid in m_ck.index else 1.0
    l_md = m_d.loc[tid, "log2FoldChange"]
    drought_sig = bool(np.isfinite(l_dck) and p_dck < 0.05 and abs(l_dck) >= 1)
    opp = np.isfinite(l_dck) and (np.sign(l_md) == -np.sign(l_dck))
    if drought_sig and opp:
        label = "recovery_confirmed"
    elif drought_sig:
        label = "drought_affected_not_reversed"
    else:
        label = "melatonin_specific_no_drought_response"
    rows.append({"transcript_id": tid, "MT_vs_drought_log2FC": round(l_md,3),
                 "drought_vs_CK_log2FC": round(l_dck,3) if np.isfinite(l_dck) else np.nan,
                 "MT_vs_CK_log2FC": round(l_mck,3) if np.isfinite(l_mck) else np.nan,
                 "drought_sig": drought_sig,
                 "MT_vs_CK_sig": bool(np.isfinite(l_mck) and p_mck < 0.05 and abs(l_mck) >= 1),
                 "class": label})
mt = pd.DataFrame(rows)
mt.to_csv(f"{D}/fig5F_melatonin_classification.tsv", sep="\t", index=False)
print("\nMelatonin categories (v2 strict):", mt["class"].value_counts().to_dict())
print("  MT_vs_CK also significant in recovery_confirmed (overshoot):",
      int(((mt["class"]=="recovery_confirmed") & mt["MT_vs_CK_sig"]).sum()))


me_c = pd.read_csv(f"{WG}/ST_DT/ME_trait_correlation.tsv", sep="\t", index_col=0)
me_p = pd.read_csv(f"{WG}/ST_DT/ME_trait_padj.tsv", sep="\t", index_col=0)
key = []
for trait in me_p.columns:
    for mod in me_p.index:
        if me_p.loc[mod, trait] < 0.05 and abs(me_c.loc[mod, trait]) >= 0.6:
            key.append({"module": mod, "trait": trait, "r": round(me_c.loc[mod, trait],3),
                        "padj": me_p.loc[mod, trait]})
keydf = pd.DataFrame(key)
gm = pd.read_csv(f"{WG}/ST_DT/gene_module_assignment.tsv", sep="\t")
msize = pd.read_csv(f"{WG}/ST_DT/module_sizes.tsv", sep="\t").set_index("module")["n_genes"]
core = gm[gm["kME_own_module"] >= 0.7]


keydf_out = keydf.merge(core.groupby("module").agg(core_members=("transcript_id", "size")),
                        left_on="module", right_index=True, how="left")
keydf_out["all_members"] = keydf_out["module"].map(msize)
keydf_out.to_csv(f"{D}/fig5D_key_module_trait_pairs.tsv", sep="\t", index=False)
keymods = set(keydf["module"])
core_out = core[core["module"].isin(keymods)]
core_out.to_csv(f"{D}/fig5E_key_module_core_members.tsv", sep="\t", index=False)
print("\nKey module-trait associations:", len(keydf_out), "; modules:", keymods)


base = pd.DataFrame({"transcript_id": sorted(union)})
base["DT_DEG"] = base["transcript_id"].isin(DT)
base["ST_DEG"] = base["transcript_id"].isin(ST)
base["in_common"] = base["transcript_id"].isin(common)
base["DT_specific"] = base["transcript_id"].isin(dt_only)
base["ST_specific"] = base["transcript_id"].isin(st_only)
base = base.merge(dir_df[["pattern"]], left_on="transcript_id", right_index=True, how="left")
base = base.merge(clusters[["transcript_id","cluster_name","pattern"]].rename(
    columns={"pattern":"cluster_pattern"}), on="transcript_id", how="left")
base = base.merge(mt[["transcript_id","class"]].rename(columns={"class":"melatonin_class"}),
                  on="transcript_id", how="left")
mod_of = gm.set_index("transcript_id")["module"]
base["ST_DT_module"] = base["transcript_id"].map(mod_of).fillna("grey/unassigned")
base.to_csv(f"{D}/stress_DEG_pattern.tsv", sep="\t", index=False)

with open(f"{A}/stats/R5_summary.md", "w", encoding="utf-8") as fh:
    fh.write(f"""# R5 Key results
- DEG cutoff reconciled against v4 across eight comparisons (padj < 0.05 and |LFC| >= 1)
- Drought DEG union {len(DT)} / salt DEG union {len(ST)} / shared {len(common)} / drought-specific {len(dt_only)} / salt-specific {len(st_only)}; heat DEG {counts.iloc[6]['n']}
- Fully concordant response direction (all four comparisons significant with the same sign): {int(consistent.sum())} transcripts
- Expression clustering (union of {len(union)} DEGs; ST_DT VST, 15 samples; Ward): k={best_k}; cluster sizes {prof.groupby('cluster').size().to_dict()}; patterns {cmap}
- Melatonin direct comparison: {len(mt)} DEGs; v2 strict classification (recovery requires a drought response and an opposing Melatonin_vs_drought response): {mt['class'].value_counts().to_dict()}; recovery_confirmed also significant in MT_vs_CK (overshoot): {int(((mt['class']=='recovery_confirmed') & mt['MT_vs_CK_sig']).sum())}
- Interpretation revised 2026-10-28: significant Melatonin_vs_drought without a significant drought baseline is melatonin_specific_no_drought_response, not recovery
- ST_DT key module-trait associations (FDR < 0.05, |r| >= 0.6): {len(keydf_out)} across modules {sorted(keymods)}
""")
print("DONE R5")
