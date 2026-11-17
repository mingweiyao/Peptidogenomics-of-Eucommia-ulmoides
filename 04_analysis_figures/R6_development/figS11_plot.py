from pathlib import Path
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))) + "/scripts_common")
from plot_style import plt, panel, save  # noqa: E402
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
A = f"{B}/output/04_article/R6_development"
D, OUT = f"{A}/data", f"{A}/figures"
os.makedirs(OUT, exist_ok=True)

UP_C, DN_C, NS_C = "#C0504D", "#4E79A7", "#BFBFBF"
STAGE = ["leaf buds", "growing leaves", "young leaves", "mature leaves"]
DCOL = {"D1": "#E15759", "D2": "#4E79A7", "D3": "#59A14F", "D4": "#F28E2B", "D5": "#B07AA1"}

sel = pd.read_csv(f"{D}/R6_candidate_ranking.tsv", sep="\t").set_index("transcript_id")
def _ids(col):
    s = sel[col]
    return list(sel[s.notna() & (s.astype(str) != "")].index)
dev_ids, pos_ids = _ids("selected_dev"), _ids("selected_pos")

def load(panel, comp):
    df = pd.read_csv(f"{MS}/01_deseq2/{panel}/{comp}_results.tsv", sep="\t")
    return df.rename(columns={df.columns[0]: "transcript_id"})

VOLS = [
    ("A", "leaf_dev", "growing_leaves_vs_leaf_buds", "Growing vs buds", dev_ids),
    ("B", "leaf_dev", "young_leaves_vs_leaf_buds", "Young vs buds", dev_ids),
    ("C", "leaf_dev", "mature_leaves_vs_leaf_buds", "Mature vs buds", dev_ids),
    ("D", "leaf_position", "close_to_tip_vs_top", "Close to tip vs top", pos_ids),
    ("E", "leaf_position", "middle_and_lower_vs_top", "Middle-lower vs top", pos_ids),
    ("F", "leaf_position", "lower_vs_top", "Lower vs top", pos_ids),
    ("G", "leaf_position", "skin_vs_top", "Skin vs top", pos_ids),
]

fig = plt.figure(figsize=(7.2, 9.6))
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.15],
                      hspace=0.42, top=0.965, bottom=0.045, left=0.075, right=0.975)
sub0 = gs[0].subgridspec(1, 3, wspace=0.38)
sub1 = gs[1].subgridspec(1, 4, wspace=0.42)

def volcano(ax, panel_dir, comp, title, hl, letter):
    d = load(panel_dir, comp)
    d = d[d["padj"].notna()].copy()
    y = -np.log10(np.clip(d["padj"].astype(float), 1e-30, 1.0))
    sig = (d["padj"] < 0.05) & (d["log2FoldChange"].abs() >= 1)
    ax.scatter(d.loc[~sig, "log2FoldChange"], y[~sig], s=3, c=NS_C, lw=0, alpha=0.75)
    for v, c in [(True, UP_C), (False, DN_C)]:
        m = sig & (d["log2FoldChange"] > 0 if v else d["log2FoldChange"] < 0)
        ax.scatter(d.loc[m, "log2FoldChange"], y[m], s=5, c=c, lw=0)
    h = d.assign(_y=y).set_index("transcript_id")
    for tid in hl:
        if tid in h.index and h.loc[tid, "padj"] < 0.05 and abs(h.loc[tid, "log2FoldChange"]) >= 1:
            ax.scatter(h.loc[tid, "log2FoldChange"], h.loc[tid, "_y"], s=20,
                       facecolor="#EE854A", edgecolor="black", lw=0.7, zorder=5)
    for v in (-1, 1):
        ax.axvline(v, color="0.4", lw=0.6, ls="--")
    ax.axhline(-np.log10(0.05), color="0.4", lw=0.6, ls="--")
    n_up = int((sig & (d["log2FoldChange"] > 0)).sum())
    n_dn = int((sig & (d["log2FoldChange"] < 0)).sum())
    ax.set_title(f"{title}\n{n_up + n_dn} DETs ({n_up} up / {n_dn} down)", fontsize=6.5)
    ax.set_xlabel("log$_2$FC vs reference", fontsize=6.5)
    ax.tick_params(labelsize=6)

for i, (letter, pdir, comp, title, hl) in enumerate(VOLS[:3]):
    ax = fig.add_subplot(sub0[0, i])
    volcano(ax, pdir, comp, title, hl, letter)
    if i == 0:
        ax.set_ylabel("$-\\log_{10}$($P_{adj}$)", fontsize=7)
    panel(ax, letter, dy=1.10)
for i, (letter, pdir, comp, title, hl) in enumerate(VOLS[3:]):
    ax = fig.add_subplot(sub1[0, i])
    volcano(ax, pdir, comp, title, hl, letter)
    if i == 0:
        ax.set_ylabel("$-\\log_{10}$($P_{adj}$)", fontsize=7)
    panel(ax, letter, dy=1.10)


axH = fig.add_subplot(gs[2])
pattern = pd.read_csv(f"{D}/leaf_dev_pattern.tsv", sep="\t")
union_ids = pattern["transcript_id"].tolist()
cl_of = dict(zip(pattern["transcript_id"], pattern["cluster_name"]))
sd = pd.read_csv(f"{MS}/00_input/sample_design.tsv", sep="\t")
dev_samp = sd.loc[sd["panel"] == "leaf_dev", "sample"].tolist()
cond_of = sd.set_index("sample").loc[dev_samp, "condition"].to_dict()
dev_samp = sorted(dev_samp, key=lambda s: STAGE.index(cond_of[s]))   
conds = [cond_of[s] for s in dev_samp]
vst = pd.read_csv(f"{MS}/02_wgcna/input_vst/leaf_dev_vst.tsv", sep="\t").set_index("transcript_id")
X = vst.loc[vst.index.isin(union_ids), dev_samp]
Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)
order = sorted(Z.index, key=lambda t: (cl_of.get(t, "D9"), Z.loc[t].mean()))
Zo = Z.loc[order]
im = axH.imshow(Zo.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
axH.set_yticks([])
axH.set_xticks(range(12), ["buds", "buds", "buds", "grow.", "grow.", "grow.",
                           "young", "young", "young", "mat.", "mat.", "mat."],
               fontsize=6, rotation=60)
for i, c in enumerate(STAGE):
    j = conds.index(c)
    axH.axvline(j - 0.5, color="black", lw=0.6)
for i, t in enumerate(Zo.index):
    axH.add_patch(mpatches.Rectangle((-1.8, i - 0.5), 1.4, 1,
                     color=DCOL.get(cl_of.get(t, "D9"), "0.8"), clip_on=False))
axH.set_ylabel("695 development DE transcripts (cluster-ordered)", fontsize=7)
cb = fig.colorbar(im, ax=axH, fraction=0.025, pad=0.03)
cb.set_label("Row Z (VST)", fontsize=6.5); cb.ax.tick_params(labelsize=5.5)
handles = [mpatches.Patch(color=DCOL[d], label=f"{d} {p}") for d, p in
           pattern.drop_duplicates("cluster_name").set_index("cluster_name")["pattern"]
           .reindex([f"D{i}" for i in range(1, 6)]).items()]
axH.legend(handles=handles, frameon=False, fontsize=5.6, loc="upper center",
           bbox_to_anchor=(0.5, -0.40), ncol=3, columnspacing=1.2, handlelength=1.4)
panel(axH, "H", dx=-0.115, dy=1.045)

save(fig, OUT, "FigS11")
print("saved FigS11 (revised)")
