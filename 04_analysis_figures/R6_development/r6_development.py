from pathlib import Path
import pandas as pd, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article/R6_development"
D = f"{A}/data"
DE = f"{B}/output/03_ms_trans_analysis/01_deseq2"
WG = f"{B}/output/03_ms_trans_analysis/02_wgcna"

def load(panel, comp):
    df = pd.read_csv(f"{DE}/{panel}/{comp}_results.tsv", sep="\t")
    return df.rename(columns={df.columns[0]: "transcript_id"})

def deg(df):
    return df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() >= 1)]

summ = pd.read_csv(f"{DE}/deseq2_summary.tsv", sep="\t").set_index("comparison")

dev_comps = ["growing_leaves_vs_leaf_buds","young_leaves_vs_leaf_buds","mature_leaves_vs_leaf_buds"]
dev_res = {c: load("leaf_dev", c) for c in dev_comps}
for c, df in dev_res.items():
    n = len(deg(df)); assert n == int(summ.loc[c,"n_sig"]), (c, n, summ.loc[c,"n_sig"])
dev_sets = {c: set(deg(df)["transcript_id"]) for c, df in dev_res.items()}
dev_union = set().union(*dev_sets.values())
dev_all3 = set.intersection(*dev_sets.values())
print(f"dev union {len(dev_union)}, all-3 {len(dev_all3)}")


pos_comps = ["close_to_tip_vs_top","middle_and_lower_vs_top","lower_vs_top","skin_vs_top"]
pos_res = {c: load("leaf_position", c) for c in pos_comps}
for c, df in pos_res.items():
    n = len(deg(df)); assert n == int(summ.loc[c,"n_sig"]), (c, n)
pos_sets = {c: set(deg(df)["transcript_id"]) for c, df in pos_res.items()}
pos_union = set().union(*pos_sets.values())

skin_specific = pos_sets["skin_vs_top"] - (pos_sets["close_to_tip_vs_top"]
                                           | pos_sets["middle_and_lower_vs_top"]
                                           | pos_sets["lower_vs_top"])
counts = pd.DataFrame([
    {"set":"dev_union(3 comps vs buds)","n":len(dev_union)},
    {"set":"dev_DEG_in_all_3","n":len(dev_all3)},
    {"set":"leafpos_union(4 comps vs top)","n":len(pos_union)},
    {"set":"dev∩leafpos","n":len(dev_union & pos_union)},
    {"set":"skin_DEG_total","n":len(pos_sets["skin_vs_top"])},
    {"set":"skin_DEG_shared_with_other_leafpos","n":len(pos_sets["skin_vs_top"] & (pos_sets["close_to_tip_vs_top"] | pos_sets["middle_and_lower_vs_top"] | pos_sets["lower_vs_top"]))},
    {"set":"skin_specific_DEG","n":len(skin_specific)},
])
counts.to_csv(f"{D}/fig6A_DE_counts.tsv", sep="\t", index=False)
print(counts.to_string(index=False))
with open(f"{D}/skin_specific_ids.txt","w") as fh:
    fh.write("\n".join(sorted(skin_specific)))


vst = pd.read_csv(f"{WG}/input_vst/leaf_dev_vst.tsv", sep="\t").set_index("transcript_id")
sd = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/sample_design.tsv", sep="\t")
dev_samples = sd.loc[sd["panel"]=="leaf_dev","sample"].tolist()
cond_of = sd.set_index("sample").loc[dev_samples,"condition"].to_dict()
print("leaf_dev conditions:", sorted(set(cond_of.values())))
X = vst.loc[vst.index.isin(dev_union), dev_samples]
Z = (X.sub(X.mean(axis=1), axis=0)).div(X.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)
Ln = linkage(pdist(Z, "euclidean"), method="ward")
best_k, prev = None, None
for k in range(2, 9):
    l = fcluster(Ln, k, criterion="maxclust")
    ss = sum(((Z.values[l==c] - Z.values[l==c].mean(axis=0))**2).sum() for c in np.unique(l)) / len(Z)
    if prev is None or ss < prev * 0.92:
        best_k, prev = k, ss
    else:
        break
lab = fcluster(Ln, best_k, criterion="maxclust")
prof = Z.copy(); prof["cluster"] = lab
long = prof.melt(id_vars="cluster", var_name="sample", value_name="z")
long["condition"] = long["sample"].map(cond_of)
stage_order = [c for c in ["leaf buds","growing leaves","young leaves","mature leaves"] if c in set(long["condition"])]
mp = long.pivot_table(index="cluster", columns="condition", values="z", aggfunc="mean")[stage_order]
mp["n"] = prof.groupby("cluster").size()
mp.index = [f"D{c}" for c in mp.index]
def name_row(r):
    v = r[stage_order].values
    if v[-1] > v[0] and all(v[i+1] >= v[i] - 0.15 for i in range(len(v)-1)): return "monotonic_up"
    if v[-1] < v[0] and all(v[i+1] <= v[i] + 0.15 for i in range(len(v)-1)): return "monotonic_down"
    if v[0] >= max(v): return "bud_high"
    if v[-1] >= max(v): return "mature_high"
    peak = int(np.argmax(v))
    return f"peak_{stage_order[peak].replace(' ','_')}"
mp["pattern"] = [name_row(mp.loc[i]) for i in mp.index]
mp.to_csv(f"{D}/fig6B_dev_cluster_mean_profiles.tsv", sep="\t")
cmap = dict(zip(mp.index, mp["pattern"]))
dev_cl = pd.DataFrame({"transcript_id": Z.index, "cluster_name": [f"D{c}" for c in lab]})
dev_cl["pattern"] = dev_cl["cluster_name"].map(cmap)
print("Development cluster k =", best_k, "\n", mp[["n","pattern"]].to_string())


pos_rows = []
for c, df in pos_res.items():
    g = deg(df)
    pos_rows.append({"comparison": c, "DEG": len(g), "up": int((g["log2FoldChange"]>0).sum()),
                     "down": int((g["log2FoldChange"]<0).sum())})
pos_sum = pd.DataFrame(pos_rows)
pos_sum.to_csv(f"{D}/fig6C_leafpos_summary.tsv", sep="\t", index=False)
print("\nLeaf position:\n", pos_sum.to_string(index=False))


r4 = pd.read_csv(f"{B}/output/04_article/R4_tissue_specificity/data/tissue_specificity_results.tsv", sep="\t")
r4m = r4.set_index("transcript_id")
inter = dev_union & pos_union
cross = pd.DataFrame({"transcript_id": sorted(inter)})
cross["dev_DEG_n"] = [sum(t in dev_sets[c] for c in dev_comps) for t in cross["transcript_id"]]
cross["leafpos_DEG_n"] = [sum(t in pos_sets[c] for c in pos_comps) for t in cross["transcript_id"]]
cross["dev_pattern"] = cross["transcript_id"].map(dev_cl.set_index("transcript_id")["pattern"]).fillna("-")
cross["tissue_specific"] = cross["transcript_id"].map(r4m["tissue_specific_group"]).fillna("not_in_R4")
cross["dominant_tissue"] = cross["transcript_id"].map(r4m["dominant_tissue"])
cross.to_csv(f"{D}/fig6D_cross_panel_intersection.tsv", sep="\t", index=False)
print(f"\nCross-panel overlap {len(inter)}; including tissue-specific candidates {cross['tissue_specific'].str.endswith('_specific').sum()}")


key_rows = []
for panel, tag in [("leaf_position","formal"), ("leaf_dev","exploratory")]:
    me_c = pd.read_csv(f"{WG}/{panel}/ME_trait_correlation.tsv", sep="\t", index_col=0)
    me_p = pd.read_csv(f"{WG}/{panel}/ME_trait_padj.tsv", sep="\t", index_col=0)
    for trait in me_p.columns:
        for mod in me_p.index:
            if me_p.loc[mod, trait] < 0.05 and abs(me_c.loc[mod, trait]) >= 0.6:
                key_rows.append({"panel": panel, "status": tag, "module": mod, "trait": trait,
                                 "r": round(me_c.loc[mod, trait],3), "padj": me_p.loc[mod, trait]})
keydf = pd.DataFrame(key_rows)
gm_pos = pd.read_csv(f"{WG}/leaf_position/gene_module_assignment.tsv", sep="\t")
gm_dev = pd.read_csv(f"{WG}/leaf_dev/gene_module_assignment.tsv", sep="\t")
allmod = pd.concat([gm_pos.assign(panel="leaf_position"), gm_dev.assign(panel="leaf_dev")])
sz = allmod.groupby(["panel","module"]).size().rename("module_size")
keydf = keydf.merge(sz, on=["panel","module"], how="left")
keydf.to_csv(f"{D}/fig6E_key_module_trait_pairs.tsv", sep="\t", index=False)
print("Key module-trait associations:", len(keydf), f"(Leaf position {sum(keydf.status=='formal')}, exploratory development {sum(keydf.status=='exploratory')})")


base = pd.DataFrame({"transcript_id": sorted(dev_union)})
base["dev_DEG_n"] = [sum(t in dev_sets[c] for c in dev_comps) for t in base["transcript_id"]]
for c in dev_comps:
    base[c] = base["transcript_id"].isin(dev_sets[c])
base = base.merge(dev_cl, on="transcript_id", how="left")
lfc_map = {}
for c, df in dev_res.items():
    lfc_map[c] = df.set_index("transcript_id")["log2FoldChange"].round(3)
for c, s in lfc_map.items():
    base[f"LFC_{c}"] = base["transcript_id"].map(s)
base["leafpos_DEG_n"] = base["transcript_id"].map(lambda t: sum(t in pos_sets[c] for c in pos_comps))
base["skin_specific"] = base["transcript_id"].isin(skin_specific)
base["module_dev"] = base["transcript_id"].map(gm_dev.set_index("transcript_id")["module"]).fillna("grey/unassigned")
base.to_csv(f"{D}/leaf_dev_pattern.tsv", sep="\t", index=False)

with open(f"{A}/stats/R6_summary.md", "w", encoding="utf-8") as fh:
    fh.write(f"""# R6 Key results
- DEG thresholds reconciled against v4 for all seven comparisons
- Development DEG union {len(dev_union)} (significant across all three comparisons: {len(dev_all3)}); leaf-position DEG union {len(pos_union)}; cross-panel overlap {len(inter)}
- **Skin-specific DEGs = {len(skin_specific)}** (skin_vs_top total {len(pos_sets['skin_vs_top'])}; {len(pos_sets['skin_vs_top'])-len(skin_specific)} shared with other leaf-position contrasts; earlier estimates treating all 345 as skin-specific were corrected)
- Developmental order (verified 2026-10-28): Li et al. 2019 (PRJNA559770; PMC6721751) defines L1 leaf buds -> L2 growing leaves (3 cm) -> L3 young fully expanded leaves -> L4 mature leaves (60 days after full expansion)
- Development expression clustering (12 VST samples, Ward): k={best_k}, patterns {cmap}
- Leaf-position comparisons: {pos_rows}
- Key module-trait associations (FDR < 0.05, |r| >= 0.6): leaf-position formal {sum(keydf.status=='formal')} + development exploratory {sum(keydf.status=='exploratory')}
- Tissue-specific candidates within the overlap: {int(cross['tissue_specific'].str.endswith('_specific').sum())}
- Interpretation: skin is a distinct tissue region, not the endpoint of a leaf-position gradient; development WGCNA at n=12 is exploratory
""")
print("DONE R6")
