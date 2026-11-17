from pathlib import Path
import sys
import pandas as pd, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts_common"))
from plot_style import *  # noqa
from scipy.stats import gaussian_kde

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article"
D = f"{A}/R2_peptide_properties/data"
OUT = f"{A}/R2_peptide_properties/figures"
props = pd.read_csv(f"{D}/peptide_properties_all_ncp.tsv", sep="\t")
comp = pd.read_csv(f"{D}/figS2_source_composition.tsv", sep="\t", index_col=0)
m = props[props["matched"] == True]
u = props[props["matched"] == False]

fig, axes = plt.subplots(2, 3, figsize=(10.5, 4.8))
ax = axes[0, 0]
cats = comp.index.tolist()
x = np.arange(len(cats))
wd = 0.38
ax.barh(x - wd/2, comp["all_pct"]*100, wd, color="#B0B0B0", label="All NCPs")
ax.barh(x + wd/2, comp["matched_pct"]*100, wd, color=NCP_C, label="Transcript-supported")
ax.set_yticks(x, [c.replace("(CDS out-of-frame)", "").strip() for c in cats], fontsize=6.5)
ax.set_xlabel("Proportion (%)"); ax.legend(frameon=False, fontsize=7)
panel(ax, "A")

for (prop, lab, rng, letter), ax in zip(
        [("length", "Peptide length (aa)", (5, 40), "B"),
         ("molecular_weight", "Molecular weight (Da)", (400, 4000), "C"),
         ("isoelectric_point", "Theoretical pI", (3.5, 12), "D"),
         ("gravy", "GRAVY", (-2.5, 2.5), "E")],
        [axes[0,1], axes[0,2], axes[1,0], axes[1,1]]):
    xs = np.linspace(*rng, 200)
    for df, col, lab2 in [(u, "#B0B0B0", "unmatched"), (m, NCP_C, "transcript-supported")]:
        v = df[prop].astype(float).values
        v = v[(v >= rng[0]) & (v <= rng[1])]
        if len(v) < 10: continue
        ax.plot(xs, gaussian_kde(v)(xs), color=col, lw=1.2, label=lab2)
        ax.fill_between(xs, gaussian_kde(v)(xs), color=col, alpha=0.18)
    ax.set_xlabel(lab, fontsize=7.5); ax.set_ylabel("Density")
    panel(ax, letter)
axes[0,1].legend(frameon=False, fontsize=7)

ax = axes[1, 2]
ax.axis("off")
ax.text(0.0, 0.9, "Sourcing note", fontsize=8, fontweight="bold", va="top")
ax.text(0.0, 0.62,
        "Density shifts between transcript-supported\n(n = 7,619) and unmatched NCPs reflect\n"
        "the composition change of the screened set \u2014\nnot a biological enrichment of properties.\n\n"
        "Source composition (panel A) shows the same\nchange at the category level\n"
        "(CDS out-of-frame 6.6% \u2192 63.0%).",
        fontsize=7, va="top", transform=ax.transAxes, linespacing=1.5)
fig.tight_layout()
for ext in ("pdf", "svg"):
    fig.savefig(f"{OUT}/FigS2.{ext}")
plt.close(fig)
print("FigS2 (A-E) done")
