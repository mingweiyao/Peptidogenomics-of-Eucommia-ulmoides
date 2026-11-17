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
pairs = pd.read_csv(f"{B}/output/02filter_trans/peptide_transcript_pairs.tsv", sep="\t")
psum = pd.read_csv(f"{B}/output/02filter_trans/peptide_summary.tsv", sep="\t")
rel = pd.read_csv(f"{D}/peptide_transcript_relative_positions.tsv", sep="\t")
lf = pd.read_csv(f"{D}/locus_isoform_features.tsv", sep="\t")
comp_tab = pd.read_csv(f"{D}/fig3A_class_code_composition.tsv", sep="\t", index_col=0)

fig = plt.figure(figsize=(7.2, 9.2))
gs = fig.add_gridspec(4, 2, height_ratios=[0.9, 1, 0.9, 1.15], hspace=0.55, wspace=0.32)


ax = fig.add_subplot(gs[0, :])
sets = ["T-bg_novel_only", "T-match", "T-focus(single)"]
pct = comp_tab[sets].div(comp_tab[sets].sum(axis=0), axis=1) * 100
main_codes = pct.max(axis=1).sort_values(ascending=False).head(8).index.tolist()
other = [c for c in pct.index if c not in main_codes]
plot_pct = pct.loc[main_codes].copy()
if other:
    orow = pct.loc[other].sum().to_frame().T
    orow.index = ["other"]
    plot_pct = pd.concat([plot_pct, orow])
colors = CLASS8 + ["#BFBFBF"]
bottom = np.zeros(3)
for i, code in enumerate(plot_pct.index):
    vals = plot_pct.loc[code, sets].values
    ax.bar([0, 1, 2], vals, bottom=bottom, color=colors[i % len(colors)], width=0.55,
           label=code if code != "other" else "other")
    for j, v in enumerate(vals):
        if v > 12:
            ax.text(j, bottom[j] + v/2, f"{code} {v:.0f}%", ha="center", va="center",
                    fontsize=5.8, color="white" if i < 6 else "black")
    bottom += vals
ax.set_xticks([0, 1, 2], ["Novel transcripts\n(background)", "NCP-compatible\n(16,595)",
                          "Single-transcript set\n(1,741)"])
ax.set_ylabel("Proportion (%)"); ax.set_ylim(0, 100)
ax.legend(frameon=False, fontsize=6, ncol=8, loc="upper center", bbox_to_anchor=(0.5, 1.28))
panel(ax, "A", dx=0.01)


bg = tf[tf["is_novel"] & ~tf["NCP_compatible"] & tf["class_code"].isin(["j","k","m","n","o","u"])]
cd = tf[tf["NCP_compatible"] & tf["class_code"].isin(["j","k","m","n","o","u"])]
codes = ["j","k","m","n","o","u"]
ax = fig.add_subplot(gs[1, 0])
data = [np.log10(bg.loc[bg["class_code"]==c, "spliced_length_nt"]) for c in codes]
bp = ax.boxplot(data, tick_labels=codes, patch_artist=True, showfliers=False, widths=0.6,
                medianprops=dict(color="black", lw=0.8))
for p in bp["boxes"]:
    p.set_facecolor("#BFD3E6"); p.set_edgecolor("black"); p.set_linewidth(0.6)
ax.set_ylabel("Spliced length (nt, log10)")
ax.set_xlabel("class code (background novel)")
panel(ax, "B")
ax = fig.add_subplot(gs[1, 1])
data = [bg.loc[bg["class_code"]==c, "exon_count"] for c in codes]
bp = ax.boxplot(data, tick_labels=codes, patch_artist=True, showfliers=False, widths=0.6,
                medianprops=dict(color="black", lw=0.8))
for p in bp["boxes"]:
    p.set_facecolor("#BFD3E6"); p.set_edgecolor("black"); p.set_linewidth(0.6)
ax.set_ylabel("Exon count")
panel(ax, "C")


ax = fig.add_subplot(gs[2, 0])
agg = pd.read_csv(f"{D}/fig3C_transcripts_per_NCP.tsv", sep="\t")
agg["bin"] = agg["bin"].astype(str).replace(">=10", "10+")
agg = agg.groupby("bin")["n_NCPs"].sum()
agg = agg.reindex([str(i) for i in range(1, 10)] + ["10+"]).dropna().reset_index()

ax.bar(agg["bin"], agg["n_NCPs"] - 1, bottom=1, color="#4E79A7")
ax.set_yscale("log")
ax.set_ylim(bottom=1)
ax.set_xlabel("Compatible transcripts per NCP")
ax.set_ylabel("NCPs (log)"); panel(ax, "D")
ax = fig.add_subplot(gs[2, 1])
agg2 = pd.read_csv(f"{D}/fig3C_NCPs_per_transcript.tsv", sep="\t")
agg2.columns = ["k", "n_transcripts"]
agg2["bin"] = agg2["k"].map(lambda k: str(int(k)) if k < 10 else "10+")
agg2 = agg2.groupby("bin")["n_transcripts"].sum()
agg2 = agg2.reindex([str(i) for i in range(1, 10)] + ["10+"]).dropna().reset_index()

ax.bar(agg2["bin"], agg2["n_transcripts"] - 1, bottom=1, color="#EE854A")
ax.set_yscale("log")
ax.set_ylim(bottom=1)
ax.set_xlabel("Distinct NCPs per transcript")
ax.set_ylabel("Transcripts (log)"); panel(ax, "E", dy=1.10)


ax = fig.add_subplot(gs[3, 0])
x = lf["total_assembled_isoforms"].clip(upper=20)
y = lf["NCP_compatible_isoforms"]
jx = x + np.random.default_rng(1).uniform(-0.25, 0.25, len(lf))
jy = y + np.random.default_rng(2).uniform(-0.25, 0.25, len(lf))
ax.hexbin(jx, jy, gridsize=28, cmap="Blues", mincnt=1, extent=(0.5, 20.5, 0.5, 10.5))
lim = 10.5
ax.plot([0.5, lim], [0.5, lim], ls="--", color="grey", lw=0.8)
ax.set_xticks(range(2, 21, 2)); ax.set_yticks(range(1, 11))
ax.set_xlabel("Total isoforms per locus (clipped 20)")
ax.set_ylabel("NCP-compatible isoforms"); panel(ax, "F")


ax = fig.add_subplot(gs[3, 1])
v = rel["relative_midpoint"].dropna()
kde = __import__("scipy.stats", fromlist=["gaussian_kde"]).gaussian_kde(v)
xs = np.linspace(0, 1, 200)
ax.plot(xs, kde(xs), color="#4E79A7", lw=1.2)
ax.fill_between(xs, kde(xs), color="#4E79A7", alpha=0.25)
ax.set_xlabel("Relative midpoint on spliced RNA (5'->3')")
ax.set_ylabel("Density"); panel(ax, "G")
fig.tight_layout()


fig2, (axl, axr) = plt.subplots(1, 2, figsize=(7.4, 3.2), gridspec_kw={"width_ratios": [1, 1.5]})
ec = rel["exon_category"].value_counts(normalize=True)*100
order_ec = ["single_exon", "first_exon", "internal_exon", "last_exon"]
vals = [ec.get(k, 0) for k in order_ec]
axl.bar(range(4), vals, color=["#4E79A7", "#F28E2B", "#9C755F", "#59A14F"])
for i, v in enumerate(vals):
    axl.text(i, v, f"{v:.1f}%", ha="center", va="bottom", fontsize=7)
axl.set_xticks(range(4), ["single", "first", "internal", "last"], fontsize=6, rotation=35)
for t in axl.get_xticklabels():
    t.set_ha("right")
axl.set_xlabel("exon category", fontsize=7)
axl.set_ylabel("Pairs (%)"); panel(axl, "H")


import importlib
exons = defaultdict(list); txinfo = {}
with open(f"{B}/output/02filter_trans/novel.gff3") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        if f[2] == "exon":
            tid = re.search(r"Parent=([^;]+)", f[8]).group(1)
            exons[tid].append((int(f[3]), int(f[4])))
ex_y = 2.4
nt = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/ncp_to_transcript.tsv", sep="\t")
u_single = nt[nt["class_code"] == "u"]["transcript_id"].iloc[0]
j_single = nt[nt["class_code"] == "j"]["transcript_id"].iloc[0]
multi = lf.sort_values("NCP_compatible_isoforms", ascending=False).iloc[0]["gene_id"]
j_multi = pairs[pairs["locus"] == multi]["transcript_id"].iloc[0]
for tag, tid in [(f"multi-isoform locus {multi}", j_multi),
                 (f"single-match j ({j_single})", j_single),
                 (f"single-match u ({u_single})", u_single)]:
    lst = sorted(exons[tid])
    lo = lst[0][0]
    span = (max(e for _, e in lst) - lo) / 1e3
    for s, e in lst:
        axr.add_patch(mpatches.Rectangle(((s - lo)/1e3, ex_y - 0.12), (e - s + 1)/1e3, 0.24,
                      facecolor="#9DC3E6", edgecolor="black", lw=0.5))
    axr.plot([0, span], [ex_y, ex_y], color="black", lw=0.6)
    pr = pairs[pairs["transcript_id"] == tid]
    for r in pr.itertuples():
        w = max((r.end - r.start + 1)/1e3, span*0.02)   
        axr.add_patch(mpatches.Rectangle(((r.start - lo)/1e3, ex_y - 0.18), w, 0.36,
                      facecolor="#C55A11", edgecolor="black", lw=0.5, alpha=0.85))
    axr.text(span, ex_y, f"  {tag}", va="center", fontsize=6.5)
    ex_y -= 1.0
axr.set_ylim(-0.9, 3.1); axr.set_yticks([])
axr.set_xlabel("Position from transcript start (kb)")
axr.set_xlim(-0.5, span * 1.35 + 1)
orange = mpatches.Patch(fc="#C55A11", label="NCP span"); blue = mpatches.Patch(fc="#9DC3E6", label="exon")
axr.legend(handles=[blue, orange], frameon=False, fontsize=6.5, loc="lower right")
panel(axr, "I")
fig2.tight_layout()
save(fig, OUT, "Fig4A-G")
save(fig2, OUT, "Fig4H-I")
print("Fig3 done")
