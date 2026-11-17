from pathlib import Path
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))) + "/scripts_common")
from plot_style import plt, panel, save  # noqa: E402

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
A = f"{B}/output/04_article/R5_stress_response"
D, OUT = f"{A}/data", f"{A}/figures"
os.makedirs(OUT, exist_ok=True)

UP_C, DN_C, NS_C = "#C0504D", "#4E79A7", "#BFBFBF"
sd = pd.read_csv(f"{MS}/00_input/sample_design.tsv", sep="\t")
st_samples = sd.loc[sd["panel"] == "ST_DT", "sample"].tolist()
cond_of = sd.set_index("sample").loc[st_samples, "condition"].to_dict()

def load_results(comp):
    df = pd.read_csv(f"{MS}/01_deseq2/ST_DT/{comp}_results.tsv", sep="\t")
    return df.rename(columns={df.columns[0]: "transcript_id"})

fig = plt.figure(figsize=(7.2, 8.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.35, 1.0],
                      hspace=0.5, wspace=0.3, top=0.96, bottom=0.05,
                      left=0.09, right=0.97)

def volcano(ax, comp, title, letter, dy=1.06):
    d = load_results(comp)
    d = d[d["padj"].notna()].copy()
    y = -np.log10(np.clip(d["padj"].astype(float), 1e-30, 1.0))
    sig = (d["padj"] < 0.05) & (d["log2FoldChange"].abs() >= 1)
    ax.scatter(d.loc[~sig, "log2FoldChange"], y[~sig], s=3, c=NS_C, lw=0, alpha=0.75)
    for v, c in [(True, UP_C), (False, DN_C)]:
        m = sig & (d["log2FoldChange"] > 0 if v else d["log2FoldChange"] < 0)
        ax.scatter(d.loc[m, "log2FoldChange"], y[m], s=5, c=c, lw=0)
    for v in (-1, 1):
        ax.axvline(v, color="0.4", lw=0.6, ls="--")
    ax.axhline(-np.log10(0.05), color="0.4", lw=0.6, ls="--")
    n_up = int((sig & (d["log2FoldChange"] > 0)).sum())
    n_dn = int((sig & (d["log2FoldChange"] < 0)).sum())
    ax.set_title(f"{title}: {n_up + n_dn} DETs ({n_up} up / {n_dn} down)", fontsize=7)
    ax.set_xlabel("log$_2$FC vs CK", fontsize=7)
    ax.set_ylabel("$-\\log_{10}$($P_{adj}$)", fontsize=7)
    ax.tick_params(labelsize=6.5)
    panel(ax, letter, dy=dy)

volcano(fig.add_subplot(gs[0, 0]), "DT_48_vs_CK", "Drought 48 h vs CK", "A")
volcano(fig.add_subplot(gs[0, 1]), "ST_48_vs_CK", "Salt 48 h vs CK", "B")


COND_ORDER = ["CK", "DT_48", "DT_144", "ST_48", "ST_144"]
st_samples = sorted(st_samples, key=lambda s: COND_ORDER.index(cond_of[s]))
pattern = pd.read_csv(f"{D}/stress_DEG_pattern.tsv", sep="\t")
union_ids = pattern["transcript_id"].tolist()
cl_of = dict(zip(pattern["transcript_id"], pattern["cluster_name"]))
vst = pd.read_csv(f"{MS}/02_wgcna/input_vst/ST_DT_vst.tsv", sep="\t").set_index("transcript_id")
X = vst.loc[vst.index.isin(union_ids), st_samples]
Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)
order = sorted(Z.index, key=lambda t: (cl_of.get(t, "C9"), Z.loc[t].mean()))
Zo = Z.loc[order]
cl_col = {"C1": "#E15759", "C2": "#F28E2B", "C3": "#59A14F", "C4": "#4E79A7"}
import matplotlib.patches as mpatches
axC = fig.add_subplot(gs[1:, 0])
im = axC.imshow(Zo.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
axC.set_yticks([])
conds = [cond_of[s] for s in st_samples]
axC.set_xticks(range(15), conds, fontsize=5, rotation=60)
seen = set()
for i, t in enumerate(Zo.index):
    c = cl_of.get(t, "C9")
    axC.add_patch(mpatches.Rectangle((-1.7, i - 0.5), 1.3, 1, color=cl_col.get(c, "0.8"),
                                     clip_on=False))
axC.set_ylabel("393 DE transcripts (cluster-ordered)", fontsize=7)
cb = fig.colorbar(im, ax=axC, fraction=0.035, pad=0.06)
cb.set_label("Row Z (VST)", fontsize=6.5); cb.ax.tick_params(labelsize=5.5)
panel(axC, "C", dx=-0.13, dy=1.04)


axD = fig.add_subplot(gs[1, 1])
prof = pd.read_csv(f"{D}/fig5C_cluster_mean_profiles.tsv", sep="\t", index_col=0)
x_dt, x_st = [0, 1, 2], [3, 4, 5]
for c in prof.index:
    col = cl_col.get(c, "0.5")
    r = prof.loc[c]
    axD.plot(x_dt, [r["CK"], r["DT_48"], r["DT_144"]], color=col, lw=1.3, marker="o", ms=3,
             label=f"{c} {r['pattern']} (n={int(r['n'])})")
    axD.plot(x_st, [r["CK"], r["ST_48"], r["ST_144"]], color=col, lw=1.1, ls="--",
             marker="s", ms=2.6)
axD.axvline(2.5, color="0.85", lw=0.6)
axD.set_xticks([0, 1, 2, 3, 4, 5], ["CK", "DT 48", "DT 144", "CK", "ST 48", "ST 144"],
               fontsize=6.5)
axD.set_xlim(-0.25, 5.4)
axD.set_ylabel("Mean row Z (VST)", fontsize=7)
axD.legend(frameon=False, fontsize=5.4, loc="lower center", bbox_to_anchor=(0.5, -0.34),
           ncol=2, columnspacing=0.9, handlelength=1.6)
axD.text(0.985, 0.985, "solid: drought series\ndashed: salt series\n(cluster names from\nmean profile vs CK)",
         transform=axD.transAxes, fontsize=5.2, ha="right", va="top", color="0.3")
panel(axD, "D", dy=1.06)


gsE = gs[2, 1].subgridspec(1, 2, wspace=0.42)
axE1 = fig.add_subplot(gsE[0, 0])
cnt = pd.read_csv(f"{D}/fig5A_DE_counts.tsv", sep="\t").set_index("set")["n"]
vals = [int(cnt["DT_specific"]), int(cnt["common(DT∩ST)"]), int(cnt["ST_specific"])]
axE1.bar([0, 1, 2], vals, color=["#2E5F8A", "#9A6B9E", "#2F7E4D"], width=0.6)
for i, v in enumerate(vals):
    axE1.text(i, v + 4, str(v), ha="center", fontsize=6.5)
axE1.set_xticks([0, 1, 2], ["Drought\nspecific", "Common", "Salt\nspecific"], fontsize=6.2)
axE1.set_ylabel("DE transcripts", fontsize=7)
panel(axE1, "E", dx=-0.26)
axE2 = fig.add_subplot(gsE[0, 1])
lfc_d = load_results("DT_144_vs_CK").set_index("transcript_id")["log2FoldChange"]
lfc_s = load_results("ST_144_vs_CK").set_index("transcript_id")["log2FoldChange"]
common_ids = pattern.loc[pattern["in_common"] == True, "transcript_id"]
x = lfc_d.reindex(common_ids); y2 = lfc_s.reindex(common_ids)
same = int((np.sign(x) == np.sign(y2)).sum())
axE2.scatter(x, y2, s=8, c=np.where(np.sign(x) == np.sign(y2), "#9A6B9E", NS_C), lw=0)
axE2.axhline(0, color="black", lw=0.6); axE2.axvline(0, color="black", lw=0.6)
axE2.set_xlabel("log$_2$FC drought 144 h vs CK", fontsize=6.5)
axE2.set_ylabel("log$_2$FC salt 144 h vs CK", fontsize=6.5)
dirM = pd.read_csv(f"{D}/fig5B_direction_matrix.tsv", sep="\t")
n4 = int(dirM["direction_consistent_across_4"].sum())


axE2.set_title(f"sign agreement at 144 h: {same}/150; {n4} in all four", fontsize=5.8)
axE2.tick_params(labelsize=6)

save(fig, OUT, "FigS8")
print("saved FigS8 (revised)")
