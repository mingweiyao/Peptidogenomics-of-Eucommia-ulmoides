from pathlib import Path
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))) + "/scripts_common")
from plot_style import plt, panel, save  # noqa: E402
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
A = f"{B}/output/04_article/R6_development"
D, OUT = f"{A}/data", f"{A}/figures"
os.makedirs(OUT, exist_ok=True)

UP_C, DN_C = "#C0504D", "#4E79A7"
DEV_C, POS_C, SKIN_C = "#2F5233", "#4E79A7", "#B07AA1"
EVID_C = {"D∩P∩S": "#9A6B9E", "D∩P": "#4E79A7", "D∩S": "#F28E2B"}
PAT_C = {"D1": "#E15759", "D2": "#4E79A7", "D3": "#59A14F", "D4": "#F28E2B", "D5": "#B07AA1"}

summ = pd.read_csv(f"{MS}/01_deseq2/deseq2_summary.tsv", sep="\t").set_index("comparison")
hm = pd.read_csv(f"{D}/R6_Fig7E_heatmap_matrix.tsv", sep="\t").set_index("transcript_id")
selE = pd.read_csv(f"{D}/R6_Fig7E_selected_candidates.tsv", sep="\t").set_index("transcript_id")
order = (selE.assign(o=selE["E_block_order"])
         .sort_values(["o", "psms_sum", "_maxscore"],
                      ascending=[True, False, False], kind="stable").index.tolist())
hm = hm.loc[order]
selE = selE.loc[order]

fig = plt.figure(figsize=(7.2, 11.9))
gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 0.92, 1.95],
                      hspace=0.52, wspace=0.42, top=0.97, bottom=0.085,
                      left=0.09, right=0.975)


axA = fig.add_subplot(gs[0, :2])
COMPS = [("growing_leaves_vs_leaf_buds", "Grow\nvs buds", "dev"),
         ("young_leaves_vs_leaf_buds", "Young\nvs buds", "dev"),
         ("mature_leaves_vs_leaf_buds", "Mature\nvs buds", "dev"),
         ("close_to_tip_vs_top", "Tip\nvs top", "pos"),
         ("middle_and_lower_vs_top", "Mid-low\nvs top", "pos"),
         ("lower_vs_top", "Lower\nvs top", "pos"),
         ("skin_vs_top", "Skin\nvs top", "skin")]
for i, (c, lab, g) in enumerate(COMPS):
    r = summ.loc[c]
    x = i + (0.28 if i >= 3 else 0) + (0.28 if i == 6 else 0)
    col = {"dev": DEV_C, "pos": POS_C, "skin": SKIN_C}[g]
    axA.bar(x, int(r["n_up"]), width=0.62, color=col, lw=0)
    axA.bar(x, -int(r["n_down"]), width=0.62, color="white", edgecolor=col, lw=0.9)
    axA.text(x, int(r["n_up"]) + 10, str(int(r["n_up"])), ha="center", fontsize=5.2)
    axA.text(x, -int(r["n_down"]) - 24, str(int(r["n_down"])), ha="center", fontsize=5.2)
axA.axhline(0, color="black", lw=0.7)
axA.axvline(2.90, color="0.6", lw=0.7, ls=":")
axA.axvline(5.98, color="0.75", lw=0.5, ls=":")
axA.set_xticks([0, 1, 2, 3.28, 4.28, 5.28, 6.84],
               [c[1] for c in COMPS], fontsize=5.6)
axA.set_ylim(-370, 340)
axA.set_ylabel("DE transcripts (DET)")
axA.text(1.0, 315, "Development\nPRJNA559770 (n = 12)", fontsize=5.4, ha="center",
         color=DEV_C)
axA.text(4.28, 315, "Leaf position/region\nPRJNA878775 (n = 15)", fontsize=5.4,
         ha="center", color=POS_C)
axA.text(0.015, 0.02, "filled: up\nopen: down\npadj<0.05\n|log$_2$FC|≥1",
         transform=axA.transAxes, fontsize=5.0, ha="left", va="bottom",
         bbox=dict(fc="white", ec="0.75", lw=0.5, pad=1.0))
panel(axA, "A", dx=-0.13)


axB = fig.add_subplot(gs[0, 2:])
prof = pd.read_csv(f"{D}/fig6B_dev_cluster_mean_profiles.tsv", sep="\t", index_col=0)
STAGE_LAB = ["Buds", "Growing", "Young", "Mature"]
ST_COLS = ["leaf buds", "growing leaves", "young leaves", "mature leaves"]
for c in prof.index:
    axB.plot(range(4), prof.loc[c, ST_COLS], marker="o", ms=3, lw=1.3, color=PAT_C[c],
             label=f"{c} {prof.loc[c,'pattern']} (n={int(prof.loc[c,'n'])})")
axB.set_xticks(range(4), STAGE_LAB, fontsize=6.5)
axB.set_ylabel("Mean row Z (VST)", fontsize=7)
axB.legend(frameon=False, fontsize=5.0, loc="upper center",
           bbox_to_anchor=(0.5, -0.20), ncol=2, columnspacing=0.8,
           title="descriptive clusters (ward, k = 5)", title_fontsize=5.0)
panel(axB, "B", dx=-0.085)


subC = gs[1, :2].subgridspec(1, 2, wspace=0.45, width_ratios=[1.15, 1])
axC1 = fig.add_subplot(subC[0, 0])
cp = pd.read_csv(f"{D}/R6_C1_position_patterns.tsv", sep="\t", index_col=0)
POS_LAB = ["Top", "Tip", "Mid-low", "Lower"]
C1_COL = {1: "#4E79A7", 2: "#2F5233", 3: "#59A14F"}
for c in cp.index:
    axC1.plot(range(4), cp.loc[c, ["top", "close to tip", "middle and lower", "lower"]],
              marker="o", ms=3, lw=1.3, color=C1_COL[c],
              label=f"P{c} (n={int(cp.loc[c,'n'])})")
axC1.axhline(0, color="0.8", lw=0.5)
axC1.set_xticks(range(4), POS_LAB, fontsize=6)
axC1.set_ylabel("Mean row Z", fontsize=7)
axC1.set_ylim(-1.95, 2.15)   
axC1.legend(frameon=False, fontsize=5.4, loc="lower right")
axC1.text(0.985, 0.985, "458 position DETs\nposition-associated\n80% subsampling stability\nmedian J = 0.80 (50 iter.)",
          transform=axC1.transAxes, fontsize=5.0, ha="right", va="top", color="0.35")
panel(axC1, "C", dx=-0.20, dy=1.10)

axC2 = fig.add_subplot(subC[0, 1])
skin = pd.read_csv(f"{MS}/01_deseq2/leaf_position/skin_vs_top_results.tsv", sep="\t")
skin = skin.rename(columns={skin.columns[0]: "tid"})
ssig = skin[(skin["padj"] < 0.05) & (skin["log2FoldChange"].abs() >= 1)]
lf = ssig["log2FoldChange"]
axC2.hist(lf[lf > 0], bins=24, color=UP_C, alpha=0.85)
axC2.hist(lf[lf < 0], bins=24, color=DN_C, alpha=0.85)
axC2.axvline(0, color="black", lw=0.6)
axC2.set_xlabel("Skin vs Top log$_2$FC", fontsize=6.5)
axC2.set_ylabel("DET count", fontsize=6.5)
axC2.tick_params(labelsize=6)
axC2.set_title("Skin vs Top (region contrast, same BioProject)\n345 DETs (146 up / 199 down);\n"
               "195 shared with position DETs", fontsize=5.2, color="0.35")

# ============================ D UpSet ============================
up = pd.read_csv(f"{D}/R6_upset_exclusive_counts.tsv", sep="\t")
_ex = pd.read_csv(f"{D}/R6_upset_exclusive_counts.tsv", sep="	").set_index("exclusive_region")["n"]
REGIONS = [("D only", int(_ex["D only"]), (1, 0, 0)), ("D∩P-S", int(_ex["D∩P∖S"]), (1, 1, 0)),
           ("D∩P∩S", int(_ex["D∩P∩S"]), (1, 1, 1)), ("P only", int(_ex["P only"]), (0, 1, 0)),
           ("D∩S-P", int(_ex["D∩S∖P"]), (1, 0, 1)), ("S only", int(_ex["S only"]), (0, 0, 1)),
           ("P∩S-D", int(_ex["P∩S∖D"]), (0, 1, 1))]

_ex = pd.read_csv(f"{D}/R6_upset_exclusive_counts.tsv", sep="\t").set_index("exclusive_region")["n"]
axD = fig.add_subplot(gs[1, 2:])
axM = axD.inset_axes([0.32, 0.00, 0.66, 0.34])     
axH = axD.inset_axes([0.32, 0.44, 0.66, 0.50])     
axS = axD.inset_axes([0.00, 0.00, 0.26, 0.34])    
for i, (name, n, memb) in enumerate(REGIONS):
    x = i
    axH.bar(x, n, width=0.62, color="#4878CF", lw=0)
    axH.text(x, n + 10, str(n), ha="center", fontsize=6)
    for r, s in enumerate(["D", "P", "S"]):
        y = 2 - r
        on = memb[r]
        axM.plot(x, y, marker="s", ms=5,
                 mfc="#333333" if on else "white", mec="#333333", mew=0.7)
    if sum(memb) > 1:
        rs = [2 - r for r, v in enumerate(memb) if v]
        axM.plot([x, x], [min(rs), max(rs)], color="#333333", lw=0.8)
axH.set_xticks([])
axH.set_ylabel("exclusive DETs", fontsize=6)
axH.tick_params(labelsize=5.5)
axH.set_ylim(0, 380)
axM.set_xticks(range(len(REGIONS)), [r[0] for r in REGIONS], fontsize=5,
               rotation=45, ha="right", rotation_mode="anchor")
axM.set_yticks([2, 1, 0], ["D", "P", "S"], fontsize=6)
axM.set_xlim(-0.6, len(REGIONS) - 0.4); axM.set_ylim(-0.6, 2.6)
axM.tick_params(length=0)
for sp in axM.spines.values():
    sp.set_visible(False)
set_sizes = [("D", 695), ("P", 458), ("S", 345)]
for r, (nm, n) in enumerate(set_sizes):
    axS.barh(2 - r, n, color="#BFBFBF", height=0.6)
    axS.text(n + 12, 2 - r, str(n), va="center", fontsize=5.6)
axS.set_xlim(0, 860); axS.set_ylim(-0.6, 2.6)
axS.set_yticks([])
axS.set_xticks([])
axS.invert_xaxis()
axS.set_title("sets", fontsize=5.6, color="0.35")
axD.axis("off")
panel(axD, "D", dx=-0.045, dy=1.12)
axD.text(0.985, 1.10, "D dev (695) · P position, no skin (458) · S skin (345)\n"
         "P & S: same BioProject, different comparisons",
         transform=axD.transAxes, fontsize=5.0, ha="right", va="top", color="0.35")


E1C = ["E1_leaf_buds", "E1_growing_leaves", "E1_young_leaves", "E1_mature_leaves"]
E2C = ["E2_top", "E2_close_to_tip", "E2_middle_and_lower", "E2_lower"]
subE = gs[2, :].subgridspec(1, 30)
axL = fig.add_subplot(subE[0, 0:9])       # E1
axR = fig.add_subplot(subE[0, 10:19])     # E2
axK = fig.add_subplot(subE[0, 21:27])     # E3
cmap_z = LinearSegmentedColormap.from_list("z", ["#2166AC", "#F7F7F7", "#B2182B"])
cmap_l = LinearSegmentedColormap.from_list("l", ["#2166AC", "#F7F7F7", "#B2182B"])
M1 = hm[E1C].values.astype(float)
M2 = hm[E2C].values.astype(float)
im1 = axL.imshow(M1, cmap=cmap_z, vmin=-2, vmax=2, aspect="auto")
axR.imshow(M2, cmap=cmap_z, vmin=-2, vmax=2, aspect="auto")
lfc = hm["E3_skin_log2FC"].values.astype(float)
sig = hm["E3_skin_sig"].values.astype(bool)
lfc_plot = np.where(np.isnan(lfc), 0, lfc)
im3 = axK.imshow(lfc_plot.reshape(-1, 1), cmap=cmap_l, vmin=-4, vmax=4, aspect="auto")
axL.set_xticks(range(4), ["Buds", "Grow", "Young", "Mature"], fontsize=5.4)
axR.set_xticks(range(4), ["Top", "Tip", "Mid", "Lower"], fontsize=5.4)
axK.set_xticks([0], ["Skin\nlog$_2$FC"], fontsize=5.4)
axL.set_yticks(range(len(hm)), hm.index, fontsize=5.0)
axR.set_yticks([]); axK.set_yticks([])
axL.set_title("E1 Development\n(row Z, log$_2$TPM+1)", fontsize=5.6)
axR.set_title("E2 Leaf position\n(row Z, own project)", fontsize=5.6)
axK.set_title("E3 Skin vs Top\n(bold sig; grey = ns)", fontsize=5.6)
for i, (v, s) in enumerate(zip(lfc, sig)):
    if np.isnan(v):
        axK.text(0, i, "n.t.", ha="center", va="center", fontsize=4.8, color="0.5")
    else:
        grey = not bool(s)   
        axK.text(0, i, f"{v:+.1f}", ha="center", va="center", fontsize=4.8,
                 color="0.45" if grey else ("white" if abs(v) > 2.6 else "black"),
                 fontweight="bold" if s else "normal")

for i, tid in enumerate(hm.index):
    axL.add_patch(mpatches.Rectangle((-1.55, i - 0.5), 0.6, 1,
                     color=EVID_C[selE.loc[tid, "evidence_class"]], clip_on=False, lw=0))
    pat = selE.loc[tid, "dev_pattern"]
    dc = {"mature_high": PAT_C["D1"], "monotonic_up": PAT_C["D2"],
          "peak_young_leaves": PAT_C["D3"], "bud_high": PAT_C["D4"],
          "monotonic_down": PAT_C["D5"]}.get(pat, "0.8")
    axL.add_patch(mpatches.Rectangle((-0.85, i - 0.5), 0.6, 1, color=dc,
                     clip_on=False, lw=0))
axL.set_xlim(-1.7, 3.5)

ev_cls = selE["evidence_class"].values
ds_idx = [i for i, c in enumerate(ev_cls) if c == "D∩S"]
if ds_idx:
    y0d, y1d = min(ds_idx) - 0.5, max(ds_idx) + 0.5
    axR.add_patch(mpatches.Rectangle((-0.5, y0d), 4, y1d - y0d, fill=False,
                    hatch="////", edgecolor="0.45", lw=0.0, alpha=0.55, zorder=5))
    axR.text(3.55, (y0d + y1d) / 2, "P\nns", fontsize=4.6, ha="left", va="center",
             color="0.35")
    axR.set_xlim(-0.5, 4.6)
for i in range(1, len(ev_cls)):
    if ev_cls[i] != ev_cls[i - 1]:
        for axx in (axL, axR, axK):
            axx.axhline(i - 0.5, color="black", lw=0.9)

for i, tid in enumerate(hm.index):
    axK.text(1.25, i, str(int(selE.loc[tid, "psms_sum"])), fontsize=4.8,
             va="center", ha="left", color="#7A3E0B")
axK.set_xlim(-0.5, 2.0)
cb1 = fig.colorbar(im1, ax=axL, fraction=0.05, pad=0.02)
cb1.set_label("row Z", fontsize=5.5); cb1.ax.tick_params(labelsize=5)
cb3 = fig.colorbar(im3, ax=axK, fraction=0.10, pad=0.28, extend="both")
cb3.set_label("log$_2$FC (colours clipped at ±4)", fontsize=5.5); cb3.ax.tick_params(labelsize=5)
leg = [mpatches.Patch(color=EVID_C[k], label=k) for k in ["D∩P∩S", "D∩P", "D∩S"]]
fig.legend(handles=leg, loc="lower left", bbox_to_anchor=(0.09, 0.012),
           ncol=3, fontsize=5.4, frameon=False, title="evidence class (left bars)",
           title_fontsize=5.2)
fig.text(0.45, 0.020, "Row Z = within-project relative pattern, NOT statistical significance; hatched E2 rows\n"
         "(D∩S-P) did not meet the position-DE threshold. E3 colours clipped at ±4; printed values\n"
         "are raw DESeq2 log$_2$FC. Second bar = developmental pattern; right numbers = PSM count.",
         fontsize=5.0, color="grey")
panel(axL, "E", dx=-0.185, dy=1.16)

save(fig, OUT, "Fig7")
print("saved Fig7 (final design v2.0, A–E)")
