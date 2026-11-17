from pathlib import Path
import pandas as pd, numpy as np, re, os
from collections import defaultdict
from scipy import stats

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article/R3_transcript_architecture"
D = f"{A}/data"
pairs = pd.read_csv(f"{B}/output/02filter_trans/peptide_transcript_pairs.tsv", sep="\t")
psum = pd.read_csv(f"{B}/output/02filter_trans/peptide_summary.tsv", sep="\t")
cls = pd.read_csv(f"{B}/output/02filter_trans/07new_transcript_quant/classification.tsv", sep="\t")
novel_ids = set(l.strip() for l in open(f"{B}/output/02filter_trans/07new_transcript_quant/novel_transcripts.txt") if l.strip())
nt = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/ncp_to_transcript.tsv", sep="\t")


tx, exons = {}, defaultdict(list)
with open(f"{B}/output/02filter_trans/novel.gff3", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        if f[2] in ("transcript", "mRNA"):
            tid = re.search(r"ID=([^;]+)", f[8]).group(1)
            m_g = re.search(r'geneID="?([^";]+)"?', f[8])
            tx[tid] = (f[0], f[6], m_g.group(1) if m_g else tid.rsplit(".",1)[0])
        elif f[2] == "exon":
            tid = re.search(r"Parent=([^;]+)", f[8]).group(1)
            exons[tid].append((int(f[3]), int(f[4])))
print(f"gff3: {len(tx):,} transcripts, {sum(len(v) for v in exons.values()):,} exons")

rows = []
spliced_len_of, exon_list_of = {}, {}
for tid, (chrom, strand, gene) in tx.items():
    lst = sorted(exons.get(tid, []))
    lens = [e - s + 1 for s, e in lst]
    sl = sum(lens)
    spliced_len_of[tid] = sl
    exon_list_of[tid] = (lst, strand)
    rows.append({"transcript_id": tid, "gene_id": gene, "chrom": chrom, "strand": strand,
                 "spliced_length_nt": sl, "exon_count": len(lst),
                 "intron_count": max(0, len(lst)-1),
                 "mean_exon_length_nt": round(np.mean(lens), 1) if lens else 0})
tf = pd.DataFrame(rows)
cls_map = cls.set_index("transcript_id")["class_code"]
tf["class_code"] = tf["transcript_id"].map(cls_map).fillna("NA")
tf["is_novel"] = tf["transcript_id"].isin(novel_ids)
comp_t = set(pairs["transcript_id"]); single_t = set(nt["transcript_id"])
tf["NCP_compatible"] = tf["transcript_id"].isin(comp_t)
tf["in_single_set"] = tf["transcript_id"].isin(single_t)
tf.to_csv(f"{D}/transcript_features.tsv", sep="\t", index=False)
print(f"transcript_features: {len(tf):,} (novel {tf['is_novel'].sum():,})")


tot_iso = cls.groupby("gene_id").size().rename("total_assembled_isoforms")
comp_by_locus = pairs.groupby("locus").agg(NCP_compatible_isoforms=("transcript_id","nunique"),
                                           distinct_NCPs=("peptide","nunique"))
lf = comp_by_locus.join(tot_iso, how="left").reset_index().rename(columns={"locus":"gene_id"})
lf["compatible_ratio"] = (lf["NCP_compatible_isoforms"]/lf["total_assembled_isoforms"]).round(4)
lf["total_assembled_isoforms"] = lf["total_assembled_isoforms"].fillna(0).astype(int)
cc_map = cls.groupby("gene_id")["class_code"].apply(lambda s: ",".join(sorted(set(s))))
lf["class_codes_all_isoforms"] = lf["gene_id"].map(cc_map)
lf.to_csv(f"{D}/locus_isoform_features.tsv", sep="\t", index=False)
print(f"locus features: {len(lf):,} loci; median total isoforms {lf['total_assembled_isoforms'].median():.0f}, "
      f"median compatible ratio {lf['compatible_ratio'].median()}")


bg = tf[tf["is_novel"] & ~tf["NCP_compatible"]]
sets = {"T-bg_novel_only": bg["class_code"], "T-match": tf.loc[tf["NCP_compatible"],"class_code"],
        "T-focus(single)": tf.loc[tf["in_single_set"],"class_code"]}
comp_tab = pd.DataFrame({k: v.value_counts() for k, v in sets.items()}).fillna(0).astype(int)
comp_tab.index.name = "class_code"
comp_tab.to_csv(f"{D}/fig3A_class_code_composition.tsv", sep="\t")
or_rows = []
for code in comp_tab.index:
    a = comp_tab.loc[code, "T-match"]; b = comp_tab["T-match"].sum() - a
    c = comp_tab.loc[code, "T-bg_novel_only"]; d = comp_tab["T-bg_novel_only"].sum() - c
    table = [[a, b], [c, d]]
    if min(a+b, c+d) == 0: continue
    orr, p = stats.fisher_exact(table, alternative="two-sided")
    or_rows.append({"class_code": code, "n_compatible": a, "n_background": c,
                    "OR": round(orr, 3), "p": p})
orf = pd.DataFrame(or_rows)
orf["padj_BH"] = np.minimum(1, orf["p"] * len(orf) / orf["p"].rank(method="min").clip(lower=1))
orf = orf.sort_values("p")
orf.to_csv(f"{D}/fig3A_class_code_enrichment.tsv", sep="\t", index=False)
print("\nclass composition:\n", comp_tab.to_string())
print("\nOR top:\n", orf.head(10).to_string(index=False))


f1 = pd.DataFrame({"k": list(psum["n_transcripts"].value_counts().sort_index().index),
                   "n_NCPs": list(psum["n_transcripts"].value_counts().sort_index().values)})
f1["bin"] = f1["k"].map(lambda k: str(k) if k < 10 else ">=10")
f1b = f1.groupby("bin", sort=False)["n_NCPs"].sum().reset_index()
f1b.to_csv(f"{D}/fig3C_transcripts_per_NCP.tsv", sep="\t", index=False)
n_per_tx = pairs.groupby("transcript_id")["peptide"].nunique().value_counts().sort_index()
f2 = pd.DataFrame({"k": n_per_tx.index, "n_transcripts": n_per_tx.values})
f2.to_csv(f"{D}/fig3C_NCPs_per_transcript.tsv", sep="\t", index=False)
print(f"\nComplexity: median transcripts per NCP: {psum['n_transcripts'].median()}, "
      f"median NCPs per transcript: {pairs.groupby('transcript_id')['peptide'].nunique().median():.0f}")


rel_rows = []
for r in pairs.itertuples():
    lst, strand = exon_list_of[r.transcript_id]
    sl = spliced_len_of[r.transcript_id]
    ps, pe = int(r.start), int(r.end)
    up, seg = 0, 0
    exon_idx = None
    for i, (s, e) in enumerate(lst):
        if s <= ps and pe <= e:
            exon_idx = i
            up += (ps - s)
            break
        up += (e - s + 1)
    if exon_idx is None:   
        rel_rows.append((r.peptide, r.transcript_id, r.class_code, r.strand, np.nan, np.nan, sl, np.nan, None, "NA"))
        continue
    seg = pe - ps + 1
    if strand == "+":
        cdna_s, cdna_e = up + 1, up + seg
    else:
        cdna_s, cdna_e = sl - up - seg + 1, sl - up
    rel_mid = ((cdna_s + cdna_e) / 2) / sl
    n_ex = len(lst)

    if strand == "-":
        exon_ord = n_ex - exon_idx          
    else:
        exon_ord = exon_idx + 1
    cat = ("single_exon" if n_ex == 1 else
           "first_exon" if exon_ord == 1 else
           "last_exon" if exon_ord == n_ex else "internal_exon")
    rel_rows.append((r.peptide, r.transcript_id, r.class_code, r.strand, cdna_s, cdna_e, sl, rel_mid, exon_ord, cat))
rel = pd.DataFrame(rel_rows, columns=["peptide","transcript_id","class_code","strand",
                                      "cdna_start","cdna_end","spliced_length_nt",
                                      "relative_midpoint","exon_index_1based","exon_category"])
ok = (rel["cdna_end"] - rel["cdna_start"] + 1 == pairs["nt_len"].values)
print(f"\nRelative position: {len(rel):,} pairs; fraction with cDNA length = 3L: {ok.mean():.4f}; relative midpoint range "
      f"{rel['relative_midpoint'].min():.3f}-{rel['relative_midpoint'].max():.3f}")
rel.to_csv(f"{D}/peptide_transcript_relative_positions.tsv", sep="\t", index=False)
rel["exon_category"].value_counts().to_csv(f"{D}/fig3E_exon_category_counts.tsv", sep="\t")


with open(f"{A}/stats/R3_summary.md", "w", encoding="utf-8") as fh:
    fh.write(f"""# R3 Key results
- Transcript features: {len(tf):,} (novel background: {tf['is_novel'].sum():,}; compatible: {tf['NCP_compatible'].sum():,}; single-match subset: {tf['in_single_set'].sum():,})
- Compatible transcripts, medians: spliced length {tf.loc[tf.NCP_compatible,'spliced_length_nt'].median():.0f} nt; exon count {tf.loc[tf.NCP_compatible,'exon_count'].median():.0f}. Novel background: {bg['spliced_length_nt'].median():.0f} nt / {bg['exon_count'].median():.0f} exons
- Loci: {len(lf):,} compatible loci; median total assembled isoforms: {lf['total_assembled_isoforms'].median():.0f}; median compatible ratio: {lf['compatible_ratio'].median()}
- Class composition (Fig3A), most common compatible classes: {comp_tab['T-match'].sort_values(ascending=False).head(4).to_dict()}
- Odds-ratio enrichment after BH correction: {orf[orf['p']<0.05/len(orf)][['class_code','OR']].to_dict('records') if len(orf) else 'See table'}
- Complexity: median compatible transcripts per NCP: {psum['n_transcripts'].median():.0f} (max {psum['n_transcripts'].max()}); median NCPs per transcript: {pairs.groupby('transcript_id')['peptide'].nunique().median():.0f}
- Relative position: {len(rel):,} matched pairs projected; correct cDNA length (3x peptide length): {ok.mean():.1%}; exon category: {rel['exon_category'].value_counts().to_dict()}
""")
print("DONE R3")
