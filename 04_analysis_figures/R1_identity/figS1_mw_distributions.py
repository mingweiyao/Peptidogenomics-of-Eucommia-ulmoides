from pathlib import Path
import sys
import pandas as pd, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts_common"))
from plot_style import *  # noqa
from scipy.stats import mannwhitneyu
from scipy.stats import gaussian_kde as kde

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article/R1_identity"
OUT = f"{A}/figures"
CP = pd.read_excel(f"{B}/output/01location/CP_peptides.xlsx", sheet_name="CP",
                   usecols=["material", "molecular_weight"])
NCP = pd.read_excel(f"{B}/output/01location/NCP_peptides.xlsx", sheet_name="NCP",
                    usecols=["material", "molecular_weight"])
CP["molecular_weight"] = pd.to_numeric(CP["molecular_weight"], errors="coerce")
NCP["molecular_weight"] = pd.to_numeric(NCP["molecular_weight"], errors="coerce")

ORDER = ["Root", "Xylem", "Seed", "Leaf"]
fig, axes = plt.subplots(2, 2, figsize=(6.4, 5.6))
stats_rows = []
for ax, tissue, letter in zip(axes.flat, ORDER, "ABCD"):
    a = NCP.loc[NCP["material"] == tissue, "molecular_weight"].dropna().values
    b = CP.loc[CP["material"] == tissue, "molecular_weight"].dropna().values
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    xs = np.linspace(400, 4200, 300)
    for v, c, lab in [(a, NCP_C, "NCP"), (b, CP_C, "CP")]:
        v = v[(v >= 400) & (v <= 4200)]
        ax.plot(xs, kde(v)(xs), color=c, lw=1.3, label=lab)
        ax.fill_between(xs, kde(v)(xs), color=c, alpha=0.18)
    ax.axvline(np.mean(a), color=NCP_C, ls="--", lw=0.8)
    ax.axvline(np.mean(b), color=CP_C, ls="--", lw=0.8)
    ax.set_xlim(400, 4200)
    ax.set_xlabel("Molecular weight (Da)", fontsize=7.5)
    if letter in ("A", "C"):
        ax.set_ylabel("Density")
    ax.set_title(tissue, fontsize=9)
    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
    ax.text(0.97, 0.96, f"Wilcoxon {sig}\nP = {p:.1e}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7)
    if letter == "A":
        ax.legend(frameon=False, fontsize=7, loc="upper left")
    panel(ax, letter)
    stats_rows.append({"tissue": tissue, "n_NCP": len(a), "n_CP": len(b),
                       "NCP_median": round(np.median(a), 1), "CP_median": round(np.median(b), 1),
                       "NCP_mean": round(np.mean(a), 1), "CP_mean": round(np.mean(b), 1),
                       "MW_U": u, "p_value": p})
fig.tight_layout()

stats = pd.DataFrame(stats_rows)
stats.to_csv(f"{A}/data/figS1_MW_wilcoxon.tsv", sep="\t", index=False)
for ext in ("pdf", "svg"):
    fig.savefig(f"{OUT}/FigS1_MW_distributions.{ext}")
plt.close(fig)
import pymupdf, io
from PIL import Image
doc = pymupdf.open(f"{OUT}/FigS1_MW_distributions.pdf")
pix = doc[0].get_pixmap(matrix=pymupdf.Matrix(600/72, 600/72))
img = Image.open(io.BytesIO(pix.tobytes("png")))
img.save(f"{OUT}/tiff_600/FigS1_MW_distributions.tiff", compression="tiff_deflate")
pix2 = doc[0].get_pixmap(matrix=pymupdf.Matrix(2.2, 2.2))
pix2.save(f"{B}/output/04_article/tmp/fig_png/R1_identity__figures__FigS1_MW_distributions.png")
doc.close()
print(stats.to_string(index=False))
print("FigS1 (our style) done")
