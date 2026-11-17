from pathlib import Path
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))) + "/scripts_common")
from plot_style import plt, panel, save, STRESS  # noqa: E402

B = str(Path(__file__).resolve().parents[2])
MS = f"{B}/output/03_ms_trans_analysis"
WG = f"{MS}/02_wgcna"
A = f"{B}/output/04_article/R5_stress_response"
D, OUT = f"{A}/data", f"{A}/figures"
os.makedirs(OUT, exist_ok=True)

fig = plt.figure(figsize=(7.2, 7.6))
gs = fig.add_gridspec(2, 3, hspace=0.52, wspace=0.38, top=0.94, bottom=0.07,
                      left=0.09, right=0.97)


axA = fig.add_subplot(gs[0, 0])
st = pd.read_csv(f"{WG}/ST_DT/softthreshold.tsv", sep="\t")
beta = int(st["chosen_beta"].iloc[0])


signed_r2 = -np.sign(st["slope"]) * st["SFT.R.sq"]
axA.plot(st["Power"], signed_r2, "o-", ms=2.5, lw=1, color="#4E79A7", label="signed $R^2$")
axA.axhline(0.85, color="0.5", lw=0.6, ls="--")
axA.axvline(beta, color="#C0504D", lw=0.8)
axA.text(beta + 0.6, 0.02, f"β = {beta}", fontsize=6.5, color="#C0504D")
axA2 = axA.twinx()
axA2.plot(st["Power"], st["mean.k."], "s--", ms=2, lw=0.9, color="#F28E2B", label="mean $k$")
axA2.tick_params(labelsize=6, colors="#F28E2B")
axA2.text(0.98, 1.06, "orange dashed: mean $k$ (right axis)", transform=axA.transAxes,
          fontsize=5.2, color="#F28E2B", ha="right", va="bottom")
axA.set_xlabel("soft-threshold power", fontsize=7)
axA.set_ylabel("scale-free $R^2$ (signed)", fontsize=7)
axA.tick_params(labelsize=6)
panel(axA, "A", dx=-0.22)


axB = fig.add_subplot(gs[0, 1])
qc = pd.read_csv(f"{WG}/ST_DT/sample_qc.tsv", sep="\t").sort_values("Z_k")
cols = [STRESS.get(c, "#9A9A9A") for c in qc["condition"]]
cols = ["#C0504D" if o == "Y" else c for o, c in zip(qc["outlier"], cols)]
axB.bar(range(len(qc)), qc["Z_k"], color=cols, width=0.65)
axB.axhline(-2, color="0.4", lw=0.6, ls="--")
axB.text(0.98, 0.06, "|Z| > 2 → outlier", transform=axB.transAxes, fontsize=5.5,
         color="0.35", ha="right")
for i, r in qc.reset_index().iterrows():
    if r["outlier"] == "Y":
        axB.annotate(r["sample"], (i, r["Z_k"]), fontsize=5.5, ha="center",
                     va="top", color="#C0504D")
axB.set_xticks(range(len(qc)),
               [c.replace("SRR1497", "…") for c in qc["sample"]], rotation=90, fontsize=4.8)
axB.set_ylabel("", fontsize=7)
axB.text(0.5, 1.05, "sample connectivity $Z_k$", transform=axB.transAxes,
         fontsize=6.5, ha="center", va="bottom")
axB.tick_params(labelsize=6)
panel(axB, "B", dx=-0.20)


axC = fig.add_subplot(gs[0, 2])
ms = pd.read_csv(f"{WG}/ST_DT/module_sizes.tsv", sep="\t")
kd = pd.read_csv(f"{D}/fig5D_key_module_trait_pairs.tsv", sep="\t")
core_of = kd.groupby("module")["core_members"].max().reindex(ms["module"]).fillna(0)
all_of = kd.groupby("module")["all_members"].max().reindex(ms["module"]).fillna(0)
all_of = all_of.where(all_of > 0, ms.set_index("module")["n_genes"]).reindex(ms["module"])
ms = ms.assign(core=core_of.values, alln=all_of.values).sort_values("n_genes")
axC.barh(range(len(ms)), ms["alln"], color="#D9E1F2", height=0.62)
axC.barh(range(len(ms)), ms["core"], color="#4E79A7", height=0.62)
axC.set_yticks(range(len(ms)), ms["module"], fontsize=5.5)
axC.set_xlabel("transcripts in module", fontsize=7)
axC.text(0.97, 0.05, "filled: core (kME ≥ 0.7)", transform=axC.transAxes,
         fontsize=5.5, ha="right", color="0.3")
axC.tick_params(labelsize=6)
panel(axC, "C", dx=-0.26)


axD = fig.add_subplot(gs[1, 0])
me_c = pd.read_csv(f"{WG}/ST_DT/ME_trait_correlation.tsv", sep="\t", index_col=0)
me_p = pd.read_csv(f"{WG}/ST_DT/ME_trait_padj.tsv", sep="\t", index_col=0)
traits = ["CK", "DT_48", "DT_144", "ST_48", "ST_144", "drought_hours", "salt_hours"]
t_lab = ["CK", "DT\n48h", "DT\n144h", "ST\n48h", "ST\n144h", "DT\nhrs", "ST\nhrs"]
R, P = me_c[traits], me_p[traits]
im = axD.imshow(R.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
axD.set_xticks(range(len(traits)), t_lab, fontsize=4.8, rotation=55,
               ha="right", rotation_mode="anchor")
axD.set_yticks(range(len(R)), R.index, fontsize=5.2)
for i in range(R.shape[0]):
    for j in range(R.shape[1]):
        if P.iloc[i, j] < 0.05:
            star = "**" if P.iloc[i, j] < 0.01 else "*"
            axD.text(j, i - 0.22, f"{R.iloc[i,j]:.1f}", ha="center", va="center", fontsize=3.6)
            axD.text(j, i + 0.22, star, ha="center", va="center", fontsize=4.5)
        else:
            axD.text(j, i, f"{R.iloc[i,j]:.1f}", ha="center", va="center",
                     fontsize=3.6, color="0.45")
cb = fig.colorbar(im, ax=axD, fraction=0.045, pad=0.03)
cb.ax.tick_params(labelsize=5)
cb.ax.text(0.5, 1.10, "$r$", transform=cb.ax.transAxes, ha="center", fontsize=6)
axD.set_title("* BH FDR < 0.05;  |r| ≥ 0.6 → key pair", fontsize=6)
panel(axD, "D", dx=-0.24, dy=1.08)


axE = fig.add_subplot(gs[1, 1])
bs = pd.read_csv(f"{WG}/ST_DT/beta_sensitivity.tsv", sep="\t")
size12 = pd.read_csv(f"{WG}/ST_DT/module_sizes.tsv", sep="\t").set_index("module")["n_genes"]
mods = sorted(bs["module"].unique(), key=lambda m: -size12.get(m, 0))
for off, c, mk in [(10, "#4E79A7", "o"), (14, "#F28E2B", "s")]:
    sub = bs[bs["beta_used"] == off].set_index("module").reindex(mods)
    axE.scatter(range(len(mods)), sub["jaccard"], s=14, marker=mk, color=c,
                label=f"β = {off}", lw=0)
axE.axhline(0.6, color="0.5", lw=0.6, ls="--")
axE.text(0.02, 0.64, "Jaccard 0.6", transform=axE.get_xaxis_transform(), fontsize=5.5,
         color="0.4")
axE.set_xticks(range(len(mods)), mods, rotation=90, fontsize=5)
axE.set_ylabel("Jaccard vs β = 12 modules", fontsize=7)
axE.legend(frameon=False, fontsize=6, loc="upper right")
axE.tick_params(labelsize=6)
panel(axE, "E", dx=-0.22)


axF = fig.add_subplot(gs[1, 2])
stab = pd.read_csv(f"{WG}/ST_DT_outlier_removed_beta12/association_stability_vs_15s.tsv", sep="\t")
kept = stab["still_FDR_sig"] == "Y"
axF.scatter(stab.loc[kept, "r_15s"], stab.loc[kept, "r_14s"], s=16, c="#4E79A7", lw=0,
            label=f"kept ({kept.sum()})")
axF.scatter(stab.loc[~kept, "r_15s"], stab.loc[~kept, "r_14s"], s=16, facecolor="none",
            edgecolor="#C0504D", lw=0.8, label=f"lost ({(~kept).sum()})")
lim = 1.08
axF.plot([-lim, lim], [-lim, lim], ls="--", color="grey", lw=0.7)
axF.set_xlim(-lim, lim); axF.set_ylim(-lim, lim)
axF.set_xlabel("$r$ (15 samples, β = 12)", fontsize=6.5)
axF.set_ylabel("$r$ (14 samples, outlier removed)", fontsize=6.5)
axF.legend(frameon=False, fontsize=5.6, loc="upper left")
axF.tick_params(labelsize=6)
panel(axF, "F", dx=-0.22)

save(fig, OUT, "FigS9")
print("saved FigS9 (revised, WGCNA diagnostics)")
