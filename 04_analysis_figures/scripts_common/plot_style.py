import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 110, "savefig.bbox": "tight",
})


CP_C = "#4878CF"      
NCP_C = "#EE854A"     
TISSUE = {"Root": "#4E79A7", "Bark": "#F28E2B", "Fruit": "#59A14F", "Leaf": "#E15759",

          "Xylem_MS": "#F28E2B", "Seed_MS": "#59A14F",
          "Xylem": "#B07AA1", "central peel": "#F28E2B", "the edges of peels": "#EDC948",
          "fruit": "#59A14F", "leaf": "#499894"}
STRESS = {"CK": "#9A9A9A", "DT_48": "#7FA6D0", "DT_144": "#2E5F8A",
          "ST_48": "#8FCB9B", "ST_144": "#2F7E4D", "heat": "#E15759", "Melatonin": "#B07AA1"}
CLASS8 = ["#4E79A7","#F28E2B","#59A14F","#E15759","#B07AA1","#76B7B2","#EDC948","#9C755F"]

def panel(ax, letter, dx=-0.12, dy=1.04):
    """Panel."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left")

def save(fig, outdir, name):
    import os
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "svg"):
        fig.savefig(f"{outdir}/{name}.{ext}")
    plt.close(fig)
    print(f"saved {outdir}/{name}.pdf/.svg")

def significance_mark(ax, x1, x2, y, h, p):
    if p < 0.001: s = "***"
    elif p < 0.01: s = "**"
    elif p < 0.05: s = "*"
    else: return
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.8, c="black")
    ax.text((x1+x2)/2, y+h, s, ha="center", va="bottom", fontsize=8)
