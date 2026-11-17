from pathlib import Path
import pandas as pd, numpy as np, os
from scipy import stats

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article/R2_peptide_properties"
D = f"{A}/data"
NCP = pd.read_excel(f"{B}/output/01location/NCP_peptides.xlsx", sheet_name="NCP")
psum = pd.read_csv(f"{B}/output/02filter_trans/peptide_summary.tsv", sep="\t")
unmat = pd.read_csv(f"{B}/output/02filter_trans/unmatched_peptides.tsv", sep="\t")

CAT_MAP = {"intergenic": "intergenic", "intron": "intron", "junction": "junction",
           "three_prime_UTR": "3'UTR", "five_prime_UTR": "5'UTR",
           "CDS": "exon_diff(CDS out-of-frame)"}
PROP = ["length", "molecular_weight", "isoelectric_point", "gravy"]


pos = NCP.drop_duplicates(subset=["sequence","chrom","start","end","strand"]).copy().reset_index(drop=True)
assert len(pos) == pos["sequence"].nunique() == 103284, "Violated one-to-one sequence-to-position relationship"
for c in PROP:
    pos[c] = pd.to_numeric(pos[c], errors="coerce")
pos["source_cat"] = pos["type"].map(CAT_MAP)
assert pos["source_cat"].notna().all(), pos.loc[pos["source_cat"].isna(),"type"].unique()
matched = set(psum["peptide"])
pos["matched"] = pos["sequence"].isin(matched)
cols_out = ["sequence","chrom","start","end","strand","frame","type","source_cat",
            "material","psms","matched"] + PROP
pos[cols_out].to_csv(f"{D}/peptide_properties_all_ncp.tsv", sep="\t", index=False)
mpos = pos[pos["matched"]].copy()
mpos.drop(columns=["matched"]).to_csv(f"{D}/peptide_properties_matched.tsv", sep="\t", index=False)
print(f"Unique genomic positions: {len(pos):,}; matched: {len(mpos):,}")


funnel = pd.DataFrame([
    {"stage":"1 NCP screening","step":"QC-passed PSM records (four criteria, both libraries)","n":110282,"unit":"rows"},
    {"stage":"1 NCP screening","step":"genome-located rows","n":109058,"unit":"rows"},
    {"stage":"1 NCP screening","step":"NCP-classified rows (material-level)","n":105966,"unit":"rows"},
    {"stage":"1 NCP screening","step":"unique NCP sequences = unique positions (P-all)","n":103284,"unit":"sequences"},
    {"stage":"2 NCP-transcript matching","step":"NCPs with >=1 candidate transcript","n":7619,"unit":"sequences"},
    {"stage":"2 NCP-transcript matching","step":"verified compatible pairs (P-T)","n":29660,"unit":"pairs"},
    {"stage":"2 NCP-transcript matching","step":"compatible transcripts (T-match)","n":16595,"unit":"transcripts"},
    {"stage":"2 NCP-transcript matching","step":"compatible loci (L-match)","n":4466,"unit":"loci"},
    {"stage":"2 NCP-transcript matching","step":"single-transcript NCPs (P-single)","n":2519,"unit":"sequences"},
    {"stage":"2 NCP-transcript matching","step":"T-focus transcripts","n":1741,"unit":"transcripts"},
])
funnel.to_csv(f"{D}/fig2A_funnel.tsv", sep="\t", index=False)


b = mpos["source_cat"].value_counts()
fig2b = pd.DataFrame({"source_cat": b.index, "n": b.values})
fig2b["pct"] = (100*fig2b["n"]/fig2b["n"].sum()).round(2)
fig2b.to_csv(f"{D}/fig2B_source_categories_matched.tsv", sep="\t", index=False)
print("\nFig2B (matched):\n", fig2b.to_string(index=False))


rows = []
for cat, g in mpos.groupby("source_cat"):
    for p, lab in [("length","Peptide length (aa)"),("molecular_weight","MW (Da)"),
                   ("isoelectric_point","pI"),("gravy","GRAVY")]:
        v = g[p].astype(float)
        rows.append({"source_cat":cat,"property":lab,"n":len(v),"mean":round(v.mean(),3),
                     "median":round(v.median(),3),"q1":round(v.quantile(.25),3),
                     "q3":round(v.quantile(.75),3)})
cd = pd.DataFrame(rows)
cd.to_csv(f"{D}/fig2CDE_by_source_describe.tsv", sep="\t", index=False)
kw_rows = []
for p in PROP:
    groups = [g[p].astype(float).values for _, g in mpos.groupby("source_cat")]
    h, pk = stats.kruskal(*groups)
    kw_rows.append({"property":p,"KW_H":round(h,2),"p":pk,"n_groups":len(groups),
                    "n_total":sum(len(g) for g in groups)})
pd.DataFrame(kw_rows).to_csv(f"{D}/fig2CDE_kruskal.tsv", sep="\t", index=False)
print("\nKW:", pd.DataFrame(kw_rows).to_string(index=False))


ch = mpos.groupby("chrom").agg(n=("sequence","size"), plus=("strand", lambda s:(s=="+").sum()),
                               minus=("strand", lambda s:(s=="-").sum())).reset_index()
ch["chr_num"] = ch["chrom"].str.extract(r"GWHBISF(\d+)").astype(float)
ch = ch.sort_values("chr_num")
ch.to_csv(f"{D}/fig2F_chromosome_counts.tsv", sep="\t", index=False)

def win_density(df, tag):
    w = df.copy()
    w["win"] = (w["start"] // 1_000_000) * 1_000_000
    g = w.groupby(["chrom","win"]).size().reset_index(name=f"n_{tag}")
    return g
wd = win_density(mpos,"matched").merge(win_density(pos,"all_ncp"), on=["chrom","win"], how="outer").fillna(0)
wd["chr_num"] = wd["chrom"].str.extract(r"GWHBISF(\d+)").astype(float)
wd = wd.sort_values(["chr_num","win"])
wd.to_csv(f"{D}/fig2F_window_density_1Mb.tsv.gz", sep="\t", index=False)


ks_rows = []
for p in PROP:
    a = mpos[p].astype(float); b = pos.loc[~pos["matched"], p].astype(float)
    ks = stats.ks_2samp(a, b)
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    ks_rows.append({"property":p,"matched_mean":round(a.mean(),3),"unmatched_mean":round(b.mean(),3),
                    "KS_D":round(ks.statistic,4),"KS_p":ks.pvalue,"MW_p":mw.pvalue})
s2 = pd.DataFrame(ks_rows)
s2.to_csv(f"{D}/figS2_matched_vs_unmatched.tsv", sep="\t", index=False)
comp_all = pos["source_cat"].value_counts(normalize=True).rename("all_pct")
comp_m = mpos["source_cat"].value_counts(normalize=True).rename("matched_pct")
comp = pd.concat([comp_all, comp_m], axis=1).round(4)
comp.to_csv(f"{D}/figS2_source_composition.tsv", sep="\t")
print("\nS2 KS:\n", s2.to_string(index=False))
print("\nS2 composition:\n", comp.to_string())


fr = pd.crosstab(pos["frame"], pos["matched"])
fr.to_csv(f"{D}/figS3A_frame_distribution.tsv", sep="\t")
st = pd.crosstab(pos["strand"], pos["matched"], normalize="columns").round(4)
st.to_csv(f"{D}/figS3A_strand_distribution.tsv", sep="\t")
print("\nFrame dist (all):\n", pos["frame"].value_counts().to_string())


with open(f"{A}/stats/R2_summary.md", "w", encoding="utf-8") as fh:
    fh.write(f"""# R2 Key results (generated 2026-10-28)

- P-all = {len(pos):,} (unique sequences equal unique positions; 1:1 confirmed; 105,966 NCP tissue-level records)
- P-match = {len(mpos):,}; unmatched unique peptides: {pos['matched'].eq(False).sum():,}
- Fig2B genomic source composition for matched peptides: {'; '.join(f"{r.source_cat} {r.n:,} ({r.pct}%)" for r in fig2b.itertuples())}
- Genomic source composition for all NCPs (comparison): {pos['source_cat'].value_counts().to_dict()}
- Fig2C-E Kruskal-Wallis tests: {', '.join(f"{r.property} H={r.KW_H} p={r.p:.2e}" for r in pd.DataFrame(kw_rows).itertuples())}
- S2 KS tests (matched vs unmatched): {', '.join(f"{r.property} D={r.KS_D} p={r.KS_p:.2e}" for r in s2.itertuples())}
- Reading frame distribution (all): {pos['frame'].value_counts().to_dict()}
- Number of chromosomes with matched peptides: {mpos['chrom'].nunique()}
- Selection summary: 105,966 tissue-level records -> 103,284 unique -> 7,619 with candidates -> 29,660 pairs -> 7,619 peptides / 16,595 transcripts / 4,466 loci -> 2,519 single matches
""")
print("\nDONE R2 ->", D)
