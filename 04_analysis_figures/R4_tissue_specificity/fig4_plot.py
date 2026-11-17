from pathlib import Path
import sys
import pandas as pd, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts_common"))
from plot_style import *  # noqa
import matplotlib.pyplot as plt

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article"
D = f"{A}/R4_tissue_specificity/data"
OUT = f"{A}/R4_tissue_specificity/figures"
res = pd.read_csv(f"{D}/tissue_specificity_results.tsv", sep="\t")
upset = pd.read_csv(f"{D}/fig4A_detection_upset.tsv", sep="\t")
per_t = pd.read_csv(f"{D}/fig4A_per_tissue_detection.tsv", sep="\t")
spec_t = pd.read_csv(f"{D}/fig4C_specific_candidates_by_tissue.tsv", sep="\t")
cls_tau = pd.read_csv(f"{D}/fig4E_classcode_tau.tsv", sep="\t")

TORDER = ["leaf", "central peel", "the edges of peels", "Xylem", "fruit"]
TCOL = {"leaf": "#499894", "central peel": "#F28E2B", "the edges of peels": "#EDC948",
        "Xylem": "#B07AA1", "fruit": "#59A14F"}
SAMPLES = ["leaf", "central peel", "the edges of peels", "Xylem", "fruit"]
tpm_cols = {"leaf": ["mean_TPM_leaf"], "central peel": ["mean_TPM_central_peel"],
            "the edges of peels": ["mean_TPM_edges_of_peels"], "Xylem": ["mean_TPM_Xylem"],
            "fruit": ["mean_TPM_fruit"]}

fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.25], hspace=0.55, wspace=0.35)


ax = fig.add_subplot(gs[0, 0])
ax.bar(range(5), per_t["n_detected"], color=[TCOL[t] for t in per_t["tissue"]])
for i, v in enumerate(per_t["n_detected"]):
    ax.text(i, v, f"{v}", ha="center", va="bottom", fontsize=7)
ax.set_xticks(range(5), ["leaf", "central\npeel", "edges of\npeels", "Xylem", "fruit"], fontsize=7)
ax.set_ylabel("Transcripts detected\n(>=2/3 reps TPM>=1)")
ax.set_ylim(0, 1150); panel(ax, "A")
ax = fig.add_subplot(gs[0, 1])
top_u = upset[upset["detected_combination"] != "none"].head(12).sort_values("n_transcripts")
ax.barh(range(len(top_u)), top_u["n_transcripts"], color="#4E79A7")
labs = [c.replace("the edges of peels", "edge").replace("central peel", "cp").replace("Xylem", "Xy").replace("fruit", "fr").replace("leaf", "lf") for c in top_u["detected_combination"]]
ax.set_yticks(range(len(top_u)), labs, fontsize=5.5)
ax.set_xlabel("Transcripts"); panel(ax, "B")


ax = fig.add_subplot(gs[1, 0])
tau = res.loc[res["n_tissues_detected"] > 0, "tau"].dropna()
ax.hist(tau, bins=30, color="#4E79A7", edgecolor="white", lw=0.3)
ax.axvline(0.8, color="#C00000", ls="--", lw=1)
ax.text(0.805, ax.get_ylim()[1]*0.9, "tau = 0.8", color="#C00000", fontsize=7)
ax.set_xlabel("Tissue specificity index (tau)"); ax.set_ylabel("Transcripts")
panel(ax, "C")


ax = fig.add_subplot(gs[1, 1])
st = spec_t.set_index("dominant_tissue").reindex(TORDER).fillna(0)
ax.bar(range(5), st["n_specific_candidates"], color=[TCOL[t] for t in TORDER])
for i, v in enumerate(st["n_specific_candidates"]):
    ax.text(i, v, f"{int(v)}", ha="center", va="bottom", fontsize=7)
ax.set_xticks(range(5), ["leaf", "central\npeel", "edges of\npeels", "Xylem", "fruit"], fontsize=7)
ax.set_ylabel("Tissue-specific candidates\n(tau >= 0.8)")
panel(ax, "D")


ax = fig.add_subplot(gs[2, :])
sd = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/sample_design.tsv", sep="\t")
tsd = sd[sd["panel"] == "tissue"].set_index("sample")["condition"]
vtpm = pd.read_csv(f"{B}/output/02filter_trans/verified_transcript_TPM.tsv", sep="\t").set_index("transcript_id")
spec = res[res["tissue_specific_group"].str.endswith("_specific")].sort_values(
    ["tau", "overall_mean_TPM"], ascending=False).head(100)
samples_sorted = [s for t in TORDER for s in tsd[tsd == t].index]
M = np.log2(vtpm.loc[spec["transcript_id"], samples_sorted].astype(float) + 1)
Mz = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)
im = ax.imshow(Mz.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_yticks([]); ax.set_xticks(range(15), [tsd[s] for s in samples_sorted], fontsize=6, rotation=45)
for i, t in enumerate(TORDER):
    ax.axvline(i*3 - 0.5, color="black", lw=0.6)
ax.set_ylabel(f"Top 100 tissue-specific transcripts\n(tau sorted)")
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cb.set_label("Row Z-score\nlog2(TPM+1)", fontsize=7)
panel(ax, "E", dx=-0.06)
fig.tight_layout()


fig2, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
cl = cls_tau[cls_tau["n"] >= 10]
bp = a1.boxplot([res.loc[(res["class_code"] == c) & (res["n_tissues_detected"] > 0), "tau"].dropna()
                 for c in cl["class_code"]],
                tick_labels=cl["class_code"], patch_artist=True, showfliers=False, widths=0.6,
                medianprops=dict(color="black", lw=0.8))
for p in bp["boxes"]:
    p.set_facecolor("#BFD3E6"); p.set_edgecolor("black"); p.set_linewidth(0.6)
a1.set_xlabel("class code"); a1.set_ylabel("tau")
a2b = a2.twinx()
a2.bar(cl["class_code"], cl["pct_tau08"], color="#F28E2B")
a2.set_ylabel("% tau >= 0.8"); a2b.spines["right"].set_visible(True)
a2.set_xlabel("class code"); panel(a1, "F"); panel(a2, "G", dx=0.0)
fig2.tight_layout()


fig3, axes = plt.subplots(1, 5, figsize=(7.2, 2.4), sharey=False)
for ax, t in zip(axes, TORDER):
    cand = res[(res["tissue_specific_group"] == f"{t}_specific")].sort_values("tau", ascending=False)
    tid = cand.iloc[0]["transcript_id"]
    cols = {"leaf": ["mean_TPM_leaf"], "central peel": ["mean_TPM_central_peel"],
            "the edges of peels": ["mean_TPM_edges_of_peels"], "Xylem": ["mean_TPM_Xylem"],
            "fruit": ["mean_TPM_fruit"]}
    xs, ys, cs = [], [], []
    for t2 in TORDER:
        s3 = tsd[tsd == t2].index
        vals = vtpm.loc[tid, s3].astype(float).values
        xs += [t2]*3; ys += list(vals); cs += [TCOL[t2]]*3
    ax.scatter(range(len(xs)), ys, c=cs, s=10)
    ax.set_xticks([0, 3, 6, 9, 12], ["leaf", "c.peel", "edges", "Xylem", "fruit"], fontsize=5.5, rotation=45)
    ax.set_yscale("symlog")
    ax.set_title(f"{tid}\n{t} (tau={cand.iloc[0]['tau']:.2f})", fontsize=6)
axes[0].set_ylabel("TPM (symlog)")
panel(axes[0], "H", dx=-0.32, dy=1.12)
fig3.tight_layout()
save(fig, OUT, "Fig5A-E")
save(fig2, OUT, "Fig5F-G")
save(fig3, OUT, "Fig5H_examples")
print("Fig4 done")
