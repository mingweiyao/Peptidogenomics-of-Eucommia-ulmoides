from pathlib import Path
import pandas as pd, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
A = f"{B}/output/04_article/R6_development"
D = f"{A}/data"
rng = np.random.default_rng(20261028)

mem = pd.read_csv(f"{D}/R6_sets_D_P_S_membership.tsv", sep="\t")
rk = pd.read_csv(f"{D}/R6_candidate_ranking.tsv", sep="\t").set_index("transcript_id")


dup = mem[(mem["in_D"]) & (mem["in_P"] | mem["in_S"])].copy()
dup = dup.set_index("transcript_id").join(
    rk[["pass_display_dev", "pass_display_pos", "psms_sum", "n_ncp_peptides",
        "all_low_complexity", "score_dev", "score_pos", "dev_pattern",
        "locus_total_isoforms", "locus_single_isoform"]])
pool = dup[(dup["pass_display_dev"] & dup["pass_display_pos"])
           & (dup["psms_sum"] >= 3) & (~dup["all_low_complexity"].astype(bool))]
print("Conservative pool:", len(pool), "| classes:",
      pool["evidence_class"].value_counts().reindex(["D∩P∩S", "D∩P", "D∩S"]).to_dict())
assert len(pool) == 71, f"Pool size {len(pool)} != 71"
pool.to_csv(f"{D}/R6_candidate_display_pool_71.tsv", sep="\t")


pool = pool.copy()
pool["_maxscore"] = pool[["score_dev", "score_pos"]].max(axis=1)
sel_rows = []
QUOTA = {"D∩P∩S": 12, "D∩P": 10, "D∩S": 8}
for cls, n in QUOTA.items():
    sub = pool[pool["evidence_class"] == cls].sort_values(
        ["psms_sum", "_maxscore", "locus_single_isoform"],
        ascending=[False, False, False], kind="stable")

    sub = sub.assign(_tid=sub.index)
    sub = sub.sort_values(["psms_sum", "_maxscore", "_tid"],
                          ascending=[False, False, True], kind="stable")
    sub = sub.drop(columns="_tid")
    take = list(sub.index[:n])
    sel_rows += take
    print(f"{cls}: pool {len(sub)} -> selected {len(take)}: {take}")
selE = pool.loc[sel_rows].copy()
selE["E_block_order"] = selE["evidence_class"].map({"D∩P∩S": 0, "D∩P": 1, "D∩S": 2})
selE.to_csv(f"{D}/R6_Fig7E_selected_candidates.tsv", sep="\t")


sd = pd.read_csv(f"{MS}/00_input/sample_design.tsv", sep="\t")
dev_samp = sd.loc[sd["panel"] == "leaf_dev", "sample"].tolist()
pos_samp = sd.loc[sd["panel"] == "leaf_position", "sample"].tolist()
cond_of = sd.set_index("sample")["condition"].to_dict()
usecols = ["transcript_id"] + dev_samp + pos_samp
tpm = pd.read_csv(f"{B}/output/02filter_trans/07new_transcript_quant/TPM_matrix.tsv",
                  sep="\t", usecols=usecols).set_index("transcript_id")
STAGES = ["leaf buds", "growing leaves", "young leaves", "mature leaves"]
POS4 = ["top", "close to tip", "middle and lower", "lower"]

def cond_mean_log2(tid, samples):
    g = {}
    for s in samples:
        g.setdefault(cond_of[s], []).append(np.log2(float(tpm.loc[tid, s]) + 1))
    return {k: float(np.mean(v)) for k, v in g.items()}

rows = []
for tid in sel_rows:
    dm = cond_mean_log2(tid, dev_samp)
    pm = cond_mean_log2(tid, pos_samp)
    dv = np.array([dm[c] for c in STAGES])
    pv = np.array([pm[c] for c in POS4])
    dz = (dv - dv.mean()) / (dv.std(ddof=0) if dv.std(ddof=0) > 0 else np.nan)
    pz = (pv - pv.mean()) / (pv.std(ddof=0) if pv.std(ddof=0) > 0 else np.nan)
    rows.append({"transcript_id": tid,
                 **{f"E1_{c.replace(' ', '_')}": v for c, v in zip(STAGES, dz)},
                 **{f"E2_{c.replace(' ', '_')}": v for c, v in zip(POS4, pz)}})
hm = pd.DataFrame(rows).set_index("transcript_id")
skin = pd.read_csv(f"{MS}/01_deseq2/leaf_position/skin_vs_top_results.tsv", sep="\t")
skin = skin.rename(columns={skin.columns[0]: "tid"}).set_index("tid")
hm["E3_skin_log2FC"] = skin.loc[sel_rows, "log2FoldChange"].round(2).values
hm["E3_skin_sig"] = ((skin.loc[sel_rows, "padj"] < 0.05)
                     & (skin.loc[sel_rows, "log2FoldChange"].abs() >= 1)).values
hm.to_csv(f"{D}/R6_Fig7E_heatmap_matrix.tsv", sep="\t")
bm = pd.DataFrame([{**{"transcript_id": tid},
                    **{f"dev_{c}": cond_mean_log2(tid, dev_samp)[c] for c in STAGES},
                    **{f"pos_{c}": cond_mean_log2(tid, pos_samp)[c] for c in POS4},
                    "skin_log2FC": skin.loc[tid, "log2FoldChange"]}
                   for tid in sel_rows]).set_index("transcript_id")
bm.to_csv(f"{D}/R6_Fig7E_block_means.tsv", sep="\t")
print("E matrix:", hm.shape)


P_ids = mem.loc[mem["in_P"], "transcript_id"]
X = rk.loc[P_ids, [f"tpm_mean_pos_{c}" for c in POS4]].astype(float)
X.columns = POS4
Zm = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1, ddof=0).replace(0, np.nan),
                                       axis=0).fillna(0)
def within_k(M, k):
    Ln = linkage(pdist(M.values, "euclidean"), method="ward")
    return fcluster(Ln, k, criterion="maxclust")
best_k, prev = None, None
for k in range(2, 7):
    lab = within_k(Zm, k)
    ss = sum(((Zm.values[lab == c] - Zm.values[lab == c].mean(0)) ** 2).sum()
             for c in np.unique(lab)) / len(Zm)
    if prev is None or ss < prev * 0.80:
        best_k, prev = k, ss
    else:
        break
lab_full = within_k(Zm, best_k)
J = []
for _ in range(50):
    idx = rng.choice(len(Zm), size=int(0.8 * len(Zm)), replace=False)
    lab_b = within_k(Zm.iloc[idx], best_k)
    dfb = pd.DataFrame({"full": lab_full[idx], "boot": lab_b})
    for c in np.unique(lab_full):
        members = set(dfb.index[dfb["full"] == c])
        best_j = 0
        for cb in np.unique(lab_b):
            b_members = set(dfb.index[dfb["boot"] == cb])
            inter = len(members & b_members)
            if inter:
                best_j = max(best_j, inter / len(members | b_members))
        J.append(best_j)
print(f"C1 leaf-position clusters k={best_k} bootstrap median Jaccard={np.median(J):.2f}")
prof = Zm.copy(); prof["cluster"] = lab_full
cprof = prof.groupby("cluster")[POS4].mean()
cprof["n"] = prof.groupby("cluster").size()
cprof.to_csv(f"{D}/R6_C1_position_patterns.tsv", sep="\t")
pat_out = pd.DataFrame({"transcript_id": Zm.index, "position_cluster": lab_full})
pat_out.to_csv(f"{D}/R6_C1_position_cluster_assignment.tsv", sep="\t", index=False)
print(cprof.round(2).to_string())
print("DONE 02")
