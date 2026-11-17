from pathlib import Path
import sys, re
import pandas as pd, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts_common"))
from plot_style import *  # noqa
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article"
D = f"{A}/R3_transcript_architecture/data"
OUT = f"{A}/R3_transcript_architecture/figures"
tf = pd.read_csv(f"{D}/transcript_features.tsv", sep="\t")
orf = pd.read_csv(f"{D}/fig3A_class_code_enrichment.tsv", sep="\t")
rel = pd.read_csv(f"{D}/peptide_transcript_relative_positions.tsv", sep="\t")
lf = pd.read_csv(f"{D}/locus_isoform_features.tsv", sep="\t")
pairs = pd.read_csv(f"{B}/output/02filter_trans/peptide_transcript_pairs.tsv", sep="\t")

codes = ["j","k","m","n","o","u"]
bg = tf[tf["is_novel"] & ~tf["NCP_compatible"] & tf["class_code"].isin(codes)]
cm = tf[tf["NCP_compatible"] & tf["class_code"].isin(codes)]

fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8))

for (col, lab, letter), ax in zip(
        [("spliced_length_nt", "Spliced length (nt, log10)", "A"),
         ("exon_count", "Exon count", "B"),
         ("mean_exon_length_nt", "Mean exon length (nt)", "C"),
         ("intron_count", "Intron count", "D")],
        [axes[0,0], axes[0,1], axes[0,2], axes[1,0]]):
    data = [np.log10(bg[col]) if col == "spliced_length_nt" else bg[col],
            np.log10(cm[col]) if col == "spliced_length_nt" else cm[col]]
    bp = ax.boxplot(data, tick_labels=["background", "NCP-comp."], patch_artist=True,
                    showfliers=False, widths=0.55, medianprops=dict(color="black", lw=0.8))
    for p, c in zip(bp["boxes"], ["#BFD3E6", "#F8CBAD"]):
        p.set_facecolor(c); p.set_edgecolor("black"); p.set_linewidth(0.6)
    ax.set_ylabel(lab, fontsize=7); panel(ax, letter, dy=1.12)


ax = axes[1, 1]
sub = orf[orf["n_compatible"] > 0].copy()
# logit CI
a = sub["n_compatible"].values
b = cm.shape[0] - a
c = sub["n_background"].values
d = bg.shape[0] - c
se = np.sqrt(1/np.maximum(a,0.5) + 1/np.maximum(b,0.5) + 1/np.maximum(c,0.5) + 1/np.maximum(d,0.5))
lo = np.exp(np.log(sub["OR"].clip(lower=1e-3)) - 1.96*se)
hi = np.exp(np.log(sub["OR"].clip(lower=1e-3)) + 1.96*se)
y = np.arange(len(sub))[::-1]
ax.errorbar(sub["OR"], y, xerr=[sub["OR"]-lo, hi-sub["OR"]], fmt="o", ms=3.5,
            color="#4E79A7", lw=1, capsize=2)
ax.axvline(1, color="grey", ls="--", lw=0.8)
ax.set_xscale("log")
ax.set_yticks(y, sub["class_code"], fontsize=7)
ax.set_xlabel("Odds ratio (NCP-compatible vs background)")
for yi, p in zip(y, sub["padj_BH"]):
    ax.text(0.02, yi, "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns")),
            transform=ax.get_yaxis_transform(), fontsize=6, va="center")
panel(ax, "E")


ax = axes[1, 2]
for i, code in enumerate(["j","k","m","n","o","u"]):
    v = rel.loc[rel["class_code"] == code, "relative_midpoint"].dropna()
    if len(v) < 5: continue
    parts = ax.violinplot([v], positions=[i], widths=0.75, showmedians=True)
    for pc in parts["bodies"]:
        pc.set_facecolor(CLASS8[i % 8]); pc.set_alpha(0.6)
ax.set_xticks(range(6), codes, fontsize=7)
ax.set_xlabel("class code"); ax.set_ylabel("Relative midpoint\non spliced RNA")
panel(ax, "F")
fig.tight_layout()
save(fig, OUT, "FigS4")


fig2, axes2 = plt.subplots(1, 3, figsize=(7.2, 2.8))

k_cnt = pairs.groupby("peptide")["transcript_id"].nunique().value_counts().sort_index()
axes2[0].bar(k_cnt.index, k_cnt.values, color="#4E79A7")
axes2[0].set_yscale("log")
axes2[0].set_xlabel("Transcripts per NCP", fontsize=7); axes2[0].set_ylabel("NCPs (log)")
panel(axes2[0], "A")

x = lf["total_assembled_isoforms"].clip(upper=20)
y = lf["NCP_compatible_isoforms"].clip(upper=10)
axes2[1].hexbin(x, y, gridsize=25, cmap="Blues", mincnt=1, extent=(0.5, 20.5, 0.5, 10.6))
axes2[1].plot([0.5, 10.6], [0.5, 10.6], ls="--", color="grey", lw=0.8)
axes2[1].set_xlabel("Total isoforms", fontsize=7); axes2[1].set_ylabel("Compatible isoforms", fontsize=7)
panel(axes2[1], "B")

exons = defaultdict(list)
nt = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/ncp_to_transcript.tsv", sep="\t")
want = {}
for code in ["j", "k", "o"]:
    want[code] = nt[nt["class_code"] == code]["transcript_id"].iloc[0]
with open(f"{B}/output/02filter_trans/novel.gff3") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9 and f[2] == "exon":
            t = re.search(r"Parent=([^;]+)", f[8]).group(1)
            if t in want.values():
                exons[t].append((int(f[3]), int(f[4])))
ax = axes2[2]
spans = {}
for code, tid in want.items():
    lst = sorted(exons[tid])
    spans[code] = (max(e for _, e in lst) - lst[0][0])/1e3
max_span = max(spans.values())
yy = 2.4
for code, tid in want.items():
    lst = sorted(exons[tid]); lo = lst[0][0]
    span = spans[code]
    for s, e in lst:
        ax.add_patch(mpatches.Rectangle(((s-lo)/1e3, yy-0.12), (e-s+1)/1e3, 0.24,
                     facecolor="#9DC3E6", edgecolor="black", lw=0.4))
    for r in pairs[pairs["transcript_id"] == tid].itertuples():
        w = max((r.end-r.start+1)/1e3, max_span*0.02)
        ax.add_patch(mpatches.Rectangle(((r.start-lo)/1e3, yy-0.18), w, 0.36,
                     facecolor="#C55A11", edgecolor="black", lw=0.4))
    ax.text(span, yy, f" {code}", va="center", fontsize=7)
    yy -= 1.0
ax.set_ylim(-0.7, 3.0); ax.set_yticks([]); ax.set_xlim(-0.5, max_span*1.08)
ax.set_xlabel("Position from start (kb)", fontsize=7)
panel(ax, "C")
fig2.tight_layout()
save(fig2, OUT, "FigS5")
print("FigS4+S5 done")
