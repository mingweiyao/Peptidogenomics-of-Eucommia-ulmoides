from pathlib import Path
import pandas as pd, numpy as np, os

B = str(Path(__file__).resolve().parents[2])
A = f"{B}/output/04_article"
NCP = pd.read_excel(f"{B}/output/01location/NCP_peptides.xlsx", sheet_name="NCP")
CP = pd.read_excel(f"{B}/output/01location/CP_peptides.xlsx", sheet_name="CP")
psum = pd.read_csv(f"{B}/output/02filter_trans/peptide_summary.tsv", sep="\t")
pairs = pd.read_csv(f"{B}/output/02filter_trans/peptide_transcript_pairs.tsv", sep="\t")
unmat = pd.read_csv(f"{B}/output/02filter_trans/unmatched_peptides.tsv", sep="\t")
vtpm = pd.read_csv(f"{B}/output/02filter_trans/verified_transcript_TPM.tsv", sep="\t")
lexpr = pd.read_csv(f"{B}/output/02filter_trans/locus_expression.tsv", sep="\t")
nt = pd.read_csv(f"{B}/output/03_ms_trans_analysis/00_input/ncp_to_transcript.tsv", sep="\t")
cc = pd.read_csv(f"{B}/output/03_ms_trans_analysis/01_deseq2/cross_candidates.tsv", sep="\t")
novel_ids = [l.strip() for l in open(f"{B}/output/02filter_trans/07new_transcript_quant/novel_transcripts.txt") if l.strip()]

rows = []
def add(set_id, unit, n, source, note=""):
    rows.append({"set": set_id, "unit": unit, "n": n, "source": source, "note": note})


ncp_pos = NCP.drop_duplicates(subset=["sequence","chrom","start","end","strand"])
add("P-all", "NCP unique peptide sequence = unique genomic position", len(ncp_pos),
    "01location/NCP_peptides.xlsx[NCP] dedup",
    f"unique seq {NCP['sequence'].nunique():,}; 1:1 check {'PASS' if ncp_pos['sequence'].nunique()==len(ncp_pos) else 'FAIL'}")
add("P-match", "matched NCP unique sequence", len(psum), "02filter_trans/peptide_summary.tsv")
add("P-unmatched", "unmatched unique sequence", unmat["peptide"].nunique(),
    "02filter_trans/unmatched_peptides.tsv",
    f"rows {len(unmat):,}; reasons {unmat['reason'].value_counts().to_dict()}")
add("P-T", "NCP-transcript compatible pairs", len(pairs), "02filter_trans/peptide_transcript_pairs.tsv")
add("T-match", "compatible unique transcripts", pairs["transcript_id"].nunique(), "pairs / verified_transcript_TPM")
add("L-match", "compatible unique loci", pairs["locus"].nunique(), "pairs / locus_expression.tsv")
add("P-single", "NCPs compatible with exactly 1 transcript", (psum["n_transcripts"]==1).sum(),
    "peptide_summary.n_transcripts==1", "Peptide set for expression analysis")
T_focus = set(nt["transcript_id"])
add("T-focus", "transcripts of P-single", len(T_focus), "03_ms_trans_analysis/00_input/ncp_to_transcript.tsv",
    f"loci {nt['locus'].nunique():,}")
add("T-bg", "assembled novel transcripts (background)", len(novel_ids),
    "07new_transcript_quant/novel_transcripts.txt", "headerless, 83,390")
add("T-full", "full merged library (ref-recovered + novel)", 108524, "07new_transcript_quant/classification.tsv")
add("DE-comparisons", "DESeq2 contrasts", 20, "01_deseq2/deseq2_summary.tsv")
add("DEG-any", "targets DEG in >=1 comparison", 1225, "README_analysis.md / deseq2_summary")
add("candidates", "integrated candidates", len(cc), "01_deseq2/cross_candidates.tsv",
    f"tiers {cc['tier'].value_counts().to_dict()}")

uni = pd.DataFrame(rows)
uni.to_csv(f"{A}/00_universe/analysis_universe.tsv", sep="\t", index=False)
print(uni.to_string(index=False))


MAT_ORDER = ["Root","Xylem","Seed","Leaf"]
cp_m = CP["material"].value_counts().reindex(MAT_ORDER).fillna(0).astype(int)
ncp_m = NCP["material"].value_counts().reindex(MAT_ORDER).fillna(0).astype(int)
r1 = pd.DataFrame({"tissue": MAT_ORDER, "CP_rows": cp_m.values, "NCP_rows": ncp_m.values})
r1.loc[len(r1)] = ["Total", len(CP), len(NCP)]
# B: CP vs NCP length per tissue (Wilcoxon rank-sum)
from scipy.stats import mannwhitneyu
rows_b = []
for m in MAT_ORDER:
    a = CP.loc[CP["material"]==m, "length"].astype(float)
    b = NCP.loc[NCP["material"]==m, "length"].astype(float)
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    rows_b.append({"tissue": m, "CP_median_len": a.median(), "NCP_median_len": b.median(),
                   "CP_mean_len": round(a.mean(),2), "NCP_mean_len": round(b.mean(),2),
                   "MW_U": u, "p_value": p})
# C/D: mean MW per tissue + overall
rows_c = []
for m in MAT_ORDER:
    rows_c.append({"tissue": m, "CP_mean_MW": round(CP.loc[CP['material']==m,'molecular_weight'].astype(float).mean(),2),
                   "NCP_mean_MW": round(NCP.loc[NCP['material']==m,'molecular_weight'].astype(float).mean(),2)})
rows_c.append({"tissue":"All","CP_mean_MW":round(CP['molecular_weight'].astype(float).mean(),2),
               "NCP_mean_MW":round(NCP['molecular_weight'].astype(float).mean(),2)})
r1b = pd.DataFrame(rows_b); r1c = pd.DataFrame(rows_c)
r1b.to_csv(f"{A}/R1_identity/fig1B_length_wilcoxon.tsv", sep="\t", index=False)
r1c.to_csv(f"{A}/R1_identity/fig1CD_mean_MW.tsv", sep="\t", index=False)
r1.to_csv(f"{A}/R1_identity/fig1A_counts.tsv", sep="\t", index=False)
print("\nR1 Fig1A:\n", r1.to_string(index=False))
print("\nR1 Fig1B (Wilcoxon):\n", r1b.to_string(index=False))
print("\nR1 Fig1CD:\n", r1c.to_string(index=False))
print("\nDONE -> 00_universe/analysis_universe.tsv + R1_identity/")
