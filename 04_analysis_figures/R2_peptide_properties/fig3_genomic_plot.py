from pathlib import Path
import sys
import pickle
import pandas as pd, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts_common"))
from plot_style import *  # noqa
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article"
D = f"{A}/R2_peptide_properties/data"
OUT = f"{A}/R2_peptide_properties/figures"
props = pd.read_csv(f"{D}/peptide_properties_matched.tsv", sep="\t")
LENS = pickle.load(open(r"../../tmp/chr_lens.pkl", "rb"))

NAME = pickle.load(open(r"../../tmp/chr_name_map.pkl", "rb"))

SHORT = {"exon_diff(CDS out-of-frame)": "CDS-OF", "intergenic": "Intergenic",
         "intron": "Intron", "junction": "Junction",
         "3'UTR": "3'UTR", "5'UTR": "5'UTR"}
props["cat_short"] = props["source_cat"].map(SHORT)
cat_n = props["cat_short"].value_counts()
cat_order = cat_n.index.tolist()
CCOL = dict(zip(cat_order, CLASS8))
m = props.copy()
m["chr_suffix"] = pd.to_numeric(m["chrom"].str.extract(r"([0-9]+)$")[0])
mm = m[m["chr_suffix"].between(485, 501)].copy()
allp = pd.read_csv(f"{D}/peptide_properties_all_ncp.tsv", sep="\t")
allp["chr_suffix"] = pd.to_numeric(allp["chrom"].str.extract(r"([0-9]+)$")[0])
aa = allp[allp["chr_suffix"].between(485, 501)].copy()
chroms = sorted([c for c in LENS if c[-3:].isdigit() and 485 <= int(c[-3:]) <= 501],
                key=lambda c: int(NAME[c][3:]))   
labels = [NAME[c] for c in chroms]

fig = plt.figure(figsize=(11.0, 11.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.5, 1.0], hspace=0.38, wspace=0.30,
                      left=0.06, right=0.97, top=0.96, bottom=0.05)


ax = fig.add_subplot(gs[0, 0])
ax.set_aspect("equal")
fr_n = mm["frame"].value_counts().reindex([1, 2, 3, 4, 5, 6]).fillna(0)
fr_lab = ["+1", "+2", "+3", "-1", "-2", "-3"]
fr_col = ["#9DC3E6", "#6FA8DC", "#2E5F8A", "#F6B26B", "#E69138", "#C0504D"]

ax.pie(cat_n.reindex(cat_order).values, radius=1.0, colors=[CCOL[c] for c in cat_order],
       startangle=90, counterclock=False,
       wedgeprops=dict(width=0.26, edgecolor="white", lw=1.0))

w2, t2 = ax.pie(fr_n.values, radius=0.70, colors=fr_col, startangle=90, counterclock=False,
       wedgeprops=dict(width=0.30, edgecolor="white", lw=1.0))
for w, lab, v in zip(w2, fr_lab, fr_n.values):
    ang = np.deg2rad((w.theta1 + w.theta2) / 2)
    r = 0.55
    ax.text(r*np.cos(ang), r*np.sin(ang), f"{lab}\n{v:,}", ha="center", va="center",
            fontsize=6.0, color="white" if lab in ("+3", "-1", "-2", "-3") else "#333333")
ax.legend(handles=[plt.Rectangle((0,0),1,1,fc=CCOL[c]) for c in cat_order],
          labels=[f"{c}  {cat_n[c]:,} ({100*cat_n[c]/len(m):.1f}%)" for c in cat_order],
          frameon=False, fontsize=6.3, loc="center left", bbox_to_anchor=(0.98, 0.5))
ax.set_title("Inner: reading frame   |   Outer: genomic source", fontsize=7.5, color="#555555", pad=14)
panel(ax, "A", dx=-0.05)


ax = fig.add_subplot(gs[0, 1])
plus_p = []
for c in chroms:
    sub = mm[mm["chrom"] == c]
    plus_p.append(100*(sub["strand"] == "+").mean())
x = np.arange(len(chroms))
ax.bar(x, plus_p, 0.75, color="#2E5F8A", label="+ strand")
ax.bar(x, [100-p for p in plus_p], 0.75, bottom=plus_p, color="#C0504D", label="- strand")
ax.set_xticks(x, labels, fontsize=6.5)
ax.set_ylabel("Strand proportion (%)")
ax.set_ylim(0, 100)
ax.legend(frameon=False, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
panel(ax, "B")


ax = fig.add_subplot(gs[1, :])
y_of = {c: i for i, c in enumerate(chroms)}
for c in chroms:
    sub = mm[mm["chrom"] == c]
    y = y_of[c]
    maxpos = sub["end"].max() / 1e6
    ax.plot([0, maxpos], [y, y], color="#C8C8C8", lw=2.4, solid_capstyle="butt", zorder=1)
    for strand, col, sgn in [("+", "#2E5F8A", 1), ("-", "#C0504D", -1)]:
        ss = sub[sub["strand"] == strand]
        ax.vlines(ss["start"]/1e6, y, y + sgn*0.30, color=col, lw=0.5, zorder=2)
    ax.text(86.5, y, f"{len(sub):,}", va="center", fontsize=7.5, color="#333333")
ax.text(86.5, len(chroms) - 0.2, "NCPs", fontsize=8, fontweight="bold")
ax.set_yticks(range(len(chroms)), labels, fontsize=7.5)
ax.set_ylim(-0.7, len(chroms) - 0.3)
ax.invert_yaxis()
ax.set_xlim(-0.5, 86)
ax.set_xlabel("Genomic position (Mb)")

ax.set_ylabel("Chromosome")
ax.legend(handles=[plt.Line2D([], [], color="#2E5F8A", lw=2),
                   plt.Line2D([], [], color="#C0504D", lw=2)],
          labels=["+ strand", "- strand"], frameon=False, fontsize=7.5,
          loc="upper right", bbox_to_anchor=(0.92, 1.10), ncol=2)
panel(ax, "C", dx=-0.055)


ax = fig.add_subplot(gs[2, 0])
win = 1_000_000
nwin = max(int(np.ceil(LENS[c]/win)) for c in chroms)
M = np.full((len(chroms), nwin), np.nan)
for i, c in enumerate(chroms):
    sub = mm[mm["chrom"] == c]
    w = (sub["start"] // win).astype(int)
    maxw = int(np.ceil(LENS[c]/win)) - 1
    for wi, cnt in w.value_counts().items():
        if wi <= maxw:
            M[i, wi] = cnt
im = ax.imshow(M, aspect="auto", cmap="Blues", vmin=0)
ax.set_yticks(range(len(chroms)), labels, fontsize=6.5)
ax.set_xticks(range(0, nwin, 10), [str(v) for v in range(0, nwin, 10)], fontsize=6.5)
ax.set_xlabel("1 Mb window (start position, Mb; range = true chromosome length)")
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.01)
cb.set_label("NCPs per window", fontsize=7)
panel(ax, "D")


def region_codes(df):
    codes = np.full(len(df), -1, dtype=np.int8)
    pos_mid = (df["start"].values + df["end"].values) / 2
    for c in chroms:
        L = LENS[c]
        sel = np.where((df["chrom"] == c).values)[0]
        mid = pos_mid[sel]
        cd = np.ones(len(sel), dtype=np.int8)
        cd[mid <= 0.1*L] = 0
        cd[mid > 0.9*L] = 2
        codes[sel] = cd
    return codes

code_m = region_codes(mm)
code_a = region_codes(aa)
n_m = int((code_m >= 0).sum())
all_cnt = np.bincount(code_a[code_a >= 0], minlength=3).astype(float)
all_frac = all_cnt / all_cnt.sum()
obs_cnt = np.bincount(code_m[code_m >= 0], minlength=3).astype(float)
exp_cnt = all_frac * n_m
ratio_pool = obs_cnt / exp_cnt
rng = np.random.default_rng(20261028)
NPERM = 1000
code_a_valid = code_a[code_a >= 0]
null = np.zeros((NPERM, 3))
for i in range(NPERM):
    idx = rng.choice(len(code_a_valid), size=n_m, replace=False)
    null[i] = np.bincount(code_a_valid[idx], minlength=3).astype(float) / (all_frac * n_m)
p_emp = [(1 + int((null[:, k] >= ratio_pool[k]).sum())) / (NPERM + 1) for k in range(3)]
regions = ["first 10%", "middle 80%", "last 10%"]

ax = fig.add_subplot(gs[2, 1])
xs = np.arange(3)
for i, r in enumerate(regions):
    per_chr = []
    for c in chroms:
        sm = mm[mm["chrom"] == c]
        sa = aa[aa["chrom"] == c]
        cd_m = region_codes(sm); cd_a = region_codes(sa)
        k = regions.index(r)
        nm = int((cd_m == k).sum()); na = int((cd_a == k).sum())
        exp_chr = (na / max(len(sa), 1)) * len(sm)
        per_chr.append(nm / exp_chr if exp_chr > 0 else np.nan)
    pts = np.array([v for v in per_chr if np.isfinite(v)], dtype=float)
    jit = np.full(len(pts), i) + rng.uniform(-0.12, 0.12, len(pts))
    ax.scatter(jit, pts, s=9, color="#B0B0B0", alpha=0.7, zorder=2)
ax.scatter(xs, ratio_pool, s=55, color=NCP_C, edgecolor="black", lw=0.7, zorder=3)
for i, (v, p) in enumerate(zip(ratio_pool, p_emp)):
    ax.text(i + 0.16, v, f"{v:.2f}\n(p_perm={p:.3f})", fontsize=6.8, va="center")
ax.axhline(1, color="grey", ls="--", lw=0.8)
ax.set_xticks(xs, regions, fontsize=7.5)
ax.set_ylabel("Observed / expected\n(same formula for grey and orange)")
panel(ax, "E")
fig.savefig(f"{OUT}/Fig3.pdf")
fig.savefig(f"{OUT}/Fig3.svg")
plt.close(fig)
print("Figure 3 rendered (A ring | B strand | C chromosomes | D heatmap | E O/E)")
