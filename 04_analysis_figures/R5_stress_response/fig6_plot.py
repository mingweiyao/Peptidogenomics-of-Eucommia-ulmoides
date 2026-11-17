from pathlib import Path
import os, pickle, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))) + "/scripts_common")
from plot_style import plt, panel, save  # noqa: E402
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram as scipy_dendrogram
from scipy.spatial.distance import pdist

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
A = f"{B}/output/04_article/R5_stress_response"
D, OUT = f"{A}/data", f"{A}/figures"
os.makedirs(OUT, exist_ok=True)

UP_C, DN_C, NS_C = "#C0504D", "#4E79A7", "#BFBFBF"
DT_SHADE, ST_SHADE = "#2E5F8A", "#2F7E4D"
HEAT_C = "#E15759"
MOD_COLORS = {"black": "#1F1F1F", "blue": "#4E79A7", "brown": "#8C510A",
              "cyan": "#00CED1", "green": "#59A14F", "greenyellow": "#BADC7C",
              "grey": "#D9D9D9", "grey60": "#8F8F8F", "lightcyan": "#D4F2F2",
              "magenta": "#C51B7D", "midnightblue": "#08306B", "pink": "#FBB4AE",
              "purple": "#7B3294", "red": "#C0504D", "salmon": "#FA8072",
              "tan": "#D2B48C", "turquoise": "#40E0D0", "yellow": "#F1C40F"}

def load_results(panel_dir, comp):
    df = pd.read_csv(f"{MS}/01_deseq2/{panel_dir}/{comp}_results.tsv", sep="\t")
    return df.rename(columns={df.columns[0]: "transcript_id"})

def _nl(p):
    return -np.log10(np.clip(p.astype(float), 1e-30, 1.0))

venn28 = pd.read_csv(f"{D}/R5_venn3_sets.tsv", sep="\t")
ids28 = set(venn28["transcript_id"])

fig = plt.figure(figsize=(7.2, 11.6))
gs = fig.add_gridspec(3, 4, height_ratios=[1.02, 1.18, 1.30],
                      hspace=0.52, wspace=0.42, top=0.965, bottom=0.05,
                      left=0.08, right=0.975)


axA = fig.add_subplot(gs[0, 0])
COMPS = [("DT_48_vs_CK", "DT\n48 h", DT_SHADE), ("DT_144_vs_CK", "DT\n144 h", DT_SHADE),
         ("ST_48_vs_CK", "ST\n48 h", ST_SHADE), ("ST_144_vs_CK", "ST\n144 h", ST_SHADE),
         ("heat_vs_heat_CK", "Heat", HEAT_C)]
for i, (comp, lab, col) in enumerate(COMPS):
    d = load_results("ST_DT" if comp.startswith(("DT", "ST")) else "heat", comp)
    s = d[(d["padj"] < 0.05) & (d["log2FoldChange"].abs() >= 1)]
    up, dn = int((s["log2FoldChange"] > 0).sum()), int((s["log2FoldChange"] < 0).sum())
    x = i + (0.32 if i == 4 else 0)
    axA.bar(x, up, width=0.62, color=col, lw=0)
    axA.bar(x, -dn, width=0.62, color="white", edgecolor=col, lw=0.9)
    axA.text(x, up + 8, str(up), ha="center", fontsize=5.8)
    axA.text(x, -dn - 20, str(dn), ha="center", fontsize=5.8)
axA.axhline(0, color="black", lw=0.7)
axA.axvline(3.82, color="0.6", lw=0.7, ls=":")
axA.set_xticks([0, 1, 2, 3, 4.32], ["DT\n48 h", "DT\n144 h", "ST\n48 h", "ST\n144 h", "Heat"],
               fontsize=6)
axA.set_ylim(-285, 165)
axA.set_ylabel("DE transcripts")
axA.text(1.5, 152, "PRJNA741324 (n = 15)", fontsize=5.6, ha="center", color="0.3")
axA.text(4.32, 152, "PRJNA1030736\n(n = 6)", fontsize=5.6, ha="center", color="0.3")
axA.text(0.02, 0.02, "filled: up\nopen: down\nDT/ST union 393\n(deduplicated)",
         transform=axA.transAxes, fontsize=5.4, ha="left", va="bottom",
         bbox=dict(fc="white", ec="0.75", lw=0.5, pad=1.2))
panel(axA, "A", dx=-0.26)


def volcano(ax, comp, panel_dir, title, letter=None):
    d = load_results(panel_dir, comp)
    d = d[d["padj"].notna()].copy()
    y = _nl(d["padj"])
    sig = (d["padj"] < 0.05) & (d["log2FoldChange"].abs() >= 1)
    ax.scatter(d.loc[~sig, "log2FoldChange"], y[~sig], s=3, c=NS_C, lw=0, alpha=0.75)
    for v, c in [(True, UP_C), (False, DN_C)]:
        m = sig & (d["log2FoldChange"] > 0 if v else d["log2FoldChange"] < 0)
        ax.scatter(d.loc[m, "log2FoldChange"], y[m], s=5, c=c, lw=0)
    h = d.assign(_y=y).set_index("transcript_id")
    hl = [t for t in ids28 if t in h.index and h.loc[t, "padj"] < 0.05
          and abs(h.loc[t, "log2FoldChange"]) >= 1]
    for t in hl:
        ax.scatter(h.loc[t, "log2FoldChange"], h.loc[t, "_y"], s=18,
                   facecolor="#EE854A", edgecolor="black", lw=0.6, zorder=5)
    for v in (-1, 1):
        ax.axvline(v, color="0.4", lw=0.6, ls="--")
    ax.axhline(-np.log10(0.05), color="0.4", lw=0.6, ls="--")
    n_up = int((sig & (d["log2FoldChange"] > 0)).sum())
    n_dn = int((sig & (d["log2FoldChange"] < 0)).sum())
    ax.set_title(f"{title}\n{n_up + n_dn} DETs ({n_up}/{n_dn})", fontsize=6.5)
    ax.set_xlabel("log$_2$FC vs control", fontsize=6.5)
    ax.tick_params(labelsize=6)
    if letter:
        panel(ax, letter, dy=1.12)

subB = gs[0, 1:].subgridspec(1, 3, wspace=0.40)
volcano(fig.add_subplot(subB[0, 0]), "DT_144_vs_CK", "ST_DT", "DT 144 h", "B")
fig.axes[-1].set_ylabel("$-\\log_{10}$($P_{adj}$)", fontsize=7)
volcano(fig.add_subplot(subB[0, 1]), "ST_144_vs_CK", "ST_DT", "ST 144 h")
volcano(fig.add_subplot(subB[0, 2]), "heat_vs_heat_CK", "heat", "Heat")





axC1 = fig.add_subplot(gs[1, :2])
axC1.axis("off")
axC1.set_title("Transcript dendrogram (WGCNA, β = 12, signed TOM; n = 1,092)",
               fontsize=6.5)
axC1.text(0.5, -0.015, "module colours: archived assignment, reproduced 100% by re-run",
          transform=axC1.transAxes, fontsize=5.0, color="0.4", ha="center", va="top")
panel(axC1, "C", dx=-0.075, dy=1.14)

axC2 = fig.add_subplot(gs[1, 2])
me_c = pd.read_csv(f"{MS}/02_wgcna/ST_DT/ME_trait_correlation.tsv", sep="\t", index_col=0)
me_p = pd.read_csv(f"{MS}/02_wgcna/ST_DT/ME_trait_padj.tsv", sep="\t", index_col=0)
stab = pd.read_csv(f"{MS}/02_wgcna/ST_DT_outlier_removed_beta12/association_stability_vs_15s.tsv",
                   sep="\t")



traits = ["DT_48", "DT_144", "ST_48", "ST_144"]
t_lab = ["DT\n48h", "DT\n144h", "ST\n48h", "ST\n144h"]
TREAT_TRAITS = traits
mods_treat = ["black", "brown", "greenyellow", "magenta", "pink", "red", "yellow"]
mods = mods_treat
R, P = me_c.loc[mods, traits], me_p.loc[mods, traits]
im = axC2.imshow(R.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
axC2.set_xticks(range(len(traits)), t_lab, fontsize=5, rotation=55,
                ha="right", rotation_mode="anchor")
ylabels = [f"{m} *" for m in mods]
axC2.set_yticks(range(len(mods)), ylabels, fontsize=5.5)
for i, m in enumerate(mods):
    axC2.get_yticklabels()[i].set_fontweight("bold")
sig_set = {(r.module_15s, r.trait): r.still_FDR_sig for _, r in stab.iterrows()}
for i, m in enumerate(mods):
    for j, t in enumerate(traits):
        key = P.iloc[i, j]
        if key < 0.05 and abs(R.iloc[i, j]) >= 0.6:
            kept = sig_set.get((m, t), "?") == "Y"
            txt_c = "white" if abs(R.iloc[i, j]) > 0.65 else "black"
            axC2.text(j, i + 0.10, f"{R.iloc[i,j]:.2f}", ha="center", va="center",
                      fontsize=4.8, color=txt_c, fontweight="bold")
            axC2.plot(j - 0.34, i - 0.30, marker="s", ms=2.4,
                      mfc="black" if kept else "none", mec=txt_c, mew=0.6, clip_on=False)
        else:

            axC2.text(j, i + 0.10, f"{R.iloc[i,j]:.2f}", ha="center", va="center",
                      fontsize=4.4, color="0.55")
axC2.set_title("Module–trait $r$ (exploratory; all values shown)\n"
               "* modules entering the panel-D set; bold = BH FDR<0.05 & |r|≥0.6\n"
               "(▪ kept / □ lost after outlier removal; CK/time traits: FigS9D)",
               fontsize=5.2, pad=3)
panel(axC2, "", dx=-0.30, dy=1.10)


axD = fig.add_subplot(gs[1, 3])
_ex = pd.read_csv(f"{D}/R5_venn3_exclusive_counts.tsv", sep="\t").set_index("exclusive_region")["n"]
R_V = {"A_only": int(_ex["DT/ST only"]), "B_only": int(_ex["WGCNA only"]),
       "C_only": int(_ex["Heat only"]), "AB": int(_ex["DT/ST∩WGCNA only"]),
       "AC": int(_ex["DT/ST∩Heat only"]), "BC": int(_ex["WGCNA∩Heat only"]),
       "ABC": int(_ex["triple"])}
r_v = 0.78
circ = [("A", (-0.46, 0.0), "#2E5F8A"), ("B", (0.46, 0.0), "#59A14F"),
        ("C", (0.0, 0.74), "#E15759")]
for _, (cx, cy), c in circ:
    axD.add_patch(mpatches.Circle((cx, cy), r_v, facecolor=c, alpha=0.32, lw=1.0,
                    edgecolor=c))
pos = {"A_only": (-0.88, -0.28), "B_only": (0.88, -0.28), "C_only": (0.0, 1.28),
       "AB": (0.0, -0.44), "AC": (-0.44, 0.50), "BC": (0.44, 0.50), "ABC": (0.0, 0.20)}
for k, (x, y) in pos.items():
    axD.text(x, y, str(R_V[k]), ha="center", va="center", fontsize=6.5)
axD.text(-1.28, 0.55, "DT/ST DE\n393", fontsize=5.8, ha="center", color="#2E5F8A")
axD.text(1.28, 0.55, "WGCNA core\n282", fontsize=5.8, ha="center", color="#4a7d3f")
axD.text(0.0, 1.72, "Heat DE\n206", fontsize=5.8, ha="center", color="#E15759")
axD.set_xlim(-1.85, 1.85); axD.set_ylim(-1.25, 2.0)
axD.set_aspect("equal"); axD.axis("off")
axD.text(1.42, -1.02, "not to scale", fontsize=5, color="0.5", ha="center")
panel(axD, "D", dx=-0.10, dy=1.16)





rk = pd.read_csv(f"{D}/R5_candidate_ranking.tsv", sep="\t").set_index("transcript_id")
v28 = venn28.set_index("transcript_id")
ST_COLS = ["tpm_mean_st_CK", "tpm_mean_st_DT_48", "tpm_mean_st_DT_144",
           "tpm_mean_st_ST_48", "tpm_mean_st_ST_144"]
HT_COLS = ["tpm_mean_ht_heat_CK", "tpm_mean_ht_heat"]
rows = []
for tid in venn28["transcript_id"]:
    k = rk.loc[tid]
    lv = np.log2(k[ST_COLS].astype(float) + 1)   
    hv = np.log2(k[HT_COLS].astype(float) + 1)   
    left = (lv - lv.mean()).values
    right = (hv - hv.mean()).values
    rows.append(list(left) + list(right))
M = pd.DataFrame(rows, columns=["CK", "DT48", "DT144", "ST48", "ST144", "HeatCK", "Heat"],
                 index=venn28["transcript_id"])
Ln = linkage(pdist(M.iloc[:, :5].values, "euclidean"), method="ward")
M = M.iloc[leaves_list(Ln)]
meta = v28.loc[M.index]

axE = fig.add_subplot(gs[2, :4])
cmap = LinearSegmentedColormap.from_list("dlt", ["#2166AC", "#F7F7F7", "#B2182B"])
data = M.values          
vlim = 2.5
im = axE.imshow(np.clip(data, -vlim, vlim), cmap=cmap, vmin=-vlim, vmax=vlim,
                aspect="auto")
axE.set_xticks(range(7), ["CK", "DT48", "DT144", "ST48", "ST144", "Heat\nCK", "Heat"],
               fontsize=6)
low_det = [not (meta.loc[t, "pass_display_st"] and meta.loc[t, "pass_display_ht"])
           for t in M.index]
y_labels = [f"{t} †" if ld else t for t, ld in zip(M.index, low_det)]
axE.set_yticks(range(len(M)), y_labels, fontsize=4.6)
for lbl, ld in zip(axE.get_yticklabels(), low_det):
    if ld:
        lbl.set_color("0.45")

axE.axvline(4.5, color="black", lw=1.2)
axE.text(2.0, -1.6, "PRJNA741324 (DT/ST)", fontsize=6, ha="center", color="0.25")
axE.text(5.5, -1.6, "PRJNA1030736 (Heat)", fontsize=6, ha="center", color="0.25")

x_mod, x_psm, x_opp = 7.15, 7.55, 7.95
for i, tid in enumerate(M.index):
    axE.add_patch(mpatches.Rectangle((x_mod, i - 0.5), 0.3, 1,
                     color=MOD_COLORS.get(meta.loc[tid, "wgcna_module"], "#D9D9D9"),
                     clip_on=False, lw=0))
    strong = meta.loc[tid, "psms_sum"] >= 3
    axE.plot(x_psm, i, marker="o", ms=3.0,
             mfc="black" if strong else "none", mec="black", mew=0.6, clip_on=False)
    if meta.loc[tid, "dtst_vs_heat_opposite"]:
        axE.plot(x_opp, i, marker="D", ms=2.6, color="#B2182B", clip_on=False)
    if not (meta.loc[tid, "pass_display_st"] and meta.loc[tid, "pass_display_ht"]):
        pass  
axE.set_xlim(-0.5, 8.4)
cb = fig.colorbar(im, ax=axE, fraction=0.018, pad=0.10, extend="both")
cb.set_label("Δ log$_2$(TPM+1) from project mean (clipped at ±2.5)", fontsize=6)
cb.ax.tick_params(labelsize=5)
axE.text(8.15, -1.6, "module | PSM | opp.", fontsize=5.2, ha="center", color="0.3")
axE.text(0.0, -0.30,
         "Colour = log$_2$(group mean TPM + 1) minus the transcript's own mean across groups of the same project\n"
         "(row-centred separately per project; not comparable across projects). Stress-associated NCP-linked\n"
         "transcripts — 28 triple-intersection candidates. ▪ module colour; ● PSM ≥ 3 (strong, 8) / ○ < 3 (low, 20);\n"
         "red diamond = opposite DT/ST-vs-Heat direction; † low detection (< 2/3 replicates TPM ≥ 1 in one project).",
         transform=axE.transAxes, fontsize=5.0, color="grey", va="top")
panel(axE, "E", dx=-0.062, dy=1.08)

save(fig, OUT, "Fig6")


import pymupdf
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
tb = fig.get_tightbbox(renderer)          
pos = axC1.get_position()
W_IN, H_IN = fig.get_figwidth(), fig.get_figheight()
pad_in = 0.1                              

def to_page(x_in, y_in):
    return ((x_in - (tb.x0 - pad_in)) * 72.0,
            ((tb.y1 + pad_in) - y_in) * 72.0)

x0, y1 = to_page(pos.x0 * W_IN, pos.y1 * H_IN)
x1, y0 = to_page(pos.x1 * W_IN, pos.y0 * H_IN)
doc = pymupdf.open(f"{OUT}/Fig6.pdf")
rdoc = pymupdf.open(f"{OUT}/fig6C1_dendro_R.pdf")
doc[0].show_pdf_page(pymupdf.Rect(x0 + 1, y1 + 1, x1 - 1, y0 - 1), rdoc, 0)
doc.save(f"{OUT}/Fig6.pdf", incremental=True, encryption=pymupdf.PDF_ENCRYPT_KEEP)

svg = doc[0].get_svg_image(matrix=pymupdf.Matrix(1, 1), text_as_path=False)
with open(f"{OUT}/Fig6.svg", "w", encoding="utf-8") as fh:
    fh.write(svg)
doc.close(); rdoc.close()
print(f"stamped R dendrogram into C1 slot ({x1 - x0:.0f}x{y0 - y1:.0f} pt); SVG regenerated with tree")

print("saved Fig6 (final design A–E; F reserved for qPCR)")
