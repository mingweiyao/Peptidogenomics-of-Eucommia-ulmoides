from pathlib import Path
import sys
import pandas as pd, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts_common"))
from plot_style import *  # noqa
from sklearn.decomposition import PCA

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article"
D = f"{A}/R4_tissue_specificity/data"
DE = f"{B}/output/03_ms_trans_analysis/01_deseq2"
OUT = f"{A}/R4_tissue_specificity/figures"

sd = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/sample_design.tsv", sep="\t")
tsd = sd[sd["panel"] == "tissue"]
vtpm = pd.read_csv(f"{B}/output/02filter_trans/verified_transcript_TPM.tsv", sep="\t").set_index("transcript_id")
res = pd.read_csv(f"{D}/tissue_specificity_results.tsv", sep="\t")
nt = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/ncp_to_transcript.tsv", sep="\t")
focus = vtpm.loc[vtpm.index.isin(set(nt["transcript_id"]))]
samples = tsd["sample"].tolist()
cond_of = dict(zip(tsd["sample"], tsd["condition"]))
TORDER = ["leaf", "central peel", "the edges of peels", "Xylem", "fruit"]
TCOL = {"leaf": "#499894", "central peel": "#F28E2B", "the edges of peels": "#EDC948",
        "Xylem": "#B07AA1", "fruit": "#59A14F"}

L = np.log2(focus[samples].astype(float) + 1)
pca = PCA(n_components=2)
P = pca.fit_transform(L.T)

fig, axes = plt.subplots(2, 4, figsize=(10.5, 5.2))
# A PCA
ax = axes[0, 0]
for c in TORDER:
    ss = tsd.loc[tsd["condition"] == c, "sample"]
    idx = [samples.index(s) for s in ss]
    ax.scatter(P[idx, 0], P[idx, 1], color=TCOL[c], s=22, label=c)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.0f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.0f}%)")
ax.legend(frameon=False, fontsize=5.5)
panel(ax, "A")

ax = axes[0, 1]
data = [L[samples].values.flatten()]
for i, c in enumerate(TORDER):
    ss = tsd.loc[tsd["condition"] == c, "sample"]
    v = L[ss].values.flatten()
    vp = ax.violinplot([v], positions=[i], widths=0.8, showmedians=True)
    vp["bodies"][0].set_facecolor(TCOL[c]); vp["bodies"][0].set_alpha(0.6)
ax.set_xticks(range(5), ["leaf","c.peel","edges","Xylem","fruit"], fontsize=6.5)
ax.set_ylabel("log2(TPM+1)"); panel(ax, "B")

comps = ["central_peel_vs_leaf", "the_edges_of_peels_vs_leaf", "Xylem_vs_leaf", "fruit_vs_leaf"]
vol_axes = [axes[0,2], axes[0,3], axes[1,0], axes[1,1]]
for k, (comp, ax) in enumerate(zip(comps, vol_axes)):
    df = pd.read_csv(f"{DE}/tissue/{comp}_results.tsv", sep="\t")
    df = df.rename(columns={df.columns[0]: "tid"})
    sig = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() >= 1)
    ax.scatter(df.loc[~sig, "log2FoldChange"], -np.log10(df.loc[~sig, "padj"].clip(lower=1e-300)),
               s=2, c="#B0B0B0", alpha=0.5)
    up = sig & (df["log2FoldChange"] > 0); dn = sig & (df["log2FoldChange"] < 0)
    ax.scatter(df.loc[up, "log2FoldChange"], -np.log10(df.loc[up, "padj"]), s=3, c="#C00000")
    ax.scatter(df.loc[dn, "log2FoldChange"], -np.log10(df.loc[dn, "padj"]), s=3, c="#0070C0")
    ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.6)
    ax.set_title(comp.replace("_vs_leaf", " vs leaf"), fontsize=6.5)
    ax.set_xlabel("log2FC", fontsize=6.5)
    if k % 2 == 0: ax.set_ylabel("-log10 padj", fontsize=6.5)
    panel(ax, chr(ord("C")+k), dx=0.0, dy=1.10)

ax = axes[1, 2]
ss_sorted = [s for c in TORDER for s in tsd.loc[tsd["condition"] == c, "sample"]]
top = res.sort_values("tau", ascending=False).head(60)
M = np.log2(vtpm.loc[top["transcript_id"], ss_sorted].astype(float) + 1)
Mz = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)
im = ax.imshow(Mz.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(15), ["lf","lf","lf","cp","cp","cp","ed","ed","ed","Xy","Xy","Xy","fr","fr","fr"], fontsize=5)
ax.set_yticks([])
cb = fig.colorbar(im, ax=ax, fraction=0.04); cb.set_label("Row Z", fontsize=6)
panel(ax, "G")

ax = axes[1, 3]
tau = res["tau"]
thr = [0.7, 0.8, 0.9]
n_spec = [int(((tau >= t) & (res["n_tissues_detected"] >= 1)).sum()) for t in thr]
ax.bar([str(t) for t in thr], n_spec, color=["#9DC3E6", "#4E79A7", "#2E5F8A"])
for i, v in enumerate(n_spec):
    ax.text(i, v, str(v), ha="center", va="bottom", fontsize=7)
ax.set_xlabel("Tau threshold"); ax.set_ylabel("Candidates")
panel(ax, "H")
fig.tight_layout()
save(fig, OUT, "FigS6_S7")
print("FigS6+S7 done")
