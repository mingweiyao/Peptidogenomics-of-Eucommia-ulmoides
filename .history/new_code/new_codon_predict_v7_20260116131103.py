import os, re
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np

# =========================================================
# CONFIG
# =========================================================
OUT_DIR = r"F:\work_mechanism\new_file\02figure\figure4\codon\codon_prediction\codon_prediction_v7"
INPUT_DIR = "F:/work_mechanism/new_file/02figure/Eu_genome_modified"
CDS_FA = os.path.join(INPUT_DIR, "Eu_CDS.fasta")
GENOME_FA = os.path.join(INPUT_DIR, "Eu_genome.fasta")
GFF3_FA = os.path.join(INPUT_DIR, "GWHBISF00000000.gff")
CANDIDATES_XLSX = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\hit_transcript_predict_orf.csv"
TRANSCRIPT_FA = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\hit_transcripts_extra.fa"
CODON_EFFICIENCY = r"F:\work_mechanism\new_file\02figure\figure4\codon\codon_efficiency.xlsx"
OUTPUT_TP_TN = os.path.join(OUT_DIR, "TP_TN_with_motif_columns.xlsx")
OUT_XLSX = os.path.join(OUT_DIR, "candidates_scored.xlsx")
SHEET_NAME = "hit_transcript_predict_orf"

UP_LEN = 200
UP_START = 200
UP_END = 0

CODON = ['ATG', 'CTG', 'ACG', 'GTG', 'TTG']
MIN_INTERNAL_NT = 30
KMER = (3, 5)
DNA_RE_200 = re.compile(r"^[ACGT]{200}$")
DNA_RE = re.compile(r"^[ACGT]+$")
MOTIFS_DNA = ["TTC", "TCT", "CTC", "TCTTC", "TCTCT"]
# =========================================================
# CDS parsing
# =========================================================
def parse_tid_details(description, cds_seq, genome_dict):
    chrom, exon_part, strand = description.split("\t")[3].split("=")[1].split(": ")
    exon_coords = [(int(s), int(e)) for s, e in (b.split("-") for b in exon_part.split(","))]
    if chrom not in genome_dict: return None
    g = genome_dict[chrom]
    cds_seq = str(cds_seq).upper()
    if cds_seq[:3] != "ATG": return None
    if strand == "+":
        start = min(s for s, _ in exon_coords)
        if start - (UP_START + 1) < 0: return None
        up_seq = g[start - (UP_START + 1): start - 1] + cds_seq
    elif strand == "-":
        end = max(e for _, e in exon_coords)
        if end + UP_START > len(g): return None
        up_seq = str(Seq(g[end: end + UP_START]).reverse_complement()) + cds_seq
    else: return None
    return up_seq.upper()
def load_cds_fasta(cds_file, genome_dict):
    cds_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        tid = record.id
        merged = parse_tid_details(record.description, str(record.seq), genome_dict)
        if merged is None: continue
        cds_dict[tid] = merged
    return cds_dict
# =========================================================
# Building TP/TN
# =========================================================
def sample_tp_tn_all_internal_atg(cds_dict, min_internal_nt=MIN_INTERNAL_NT):
    rows_tp = []
    rows_tn = []
    for tid, seq in cds_dict.items():
        seq = seq.upper()
        if seq[UP_START:UP_START+3] != "ATG": continue
        kozak_tp = seq[UP_START-4: UP_START+4]
        rows_tp.append(('ATG', tid, UP_START, seq, 1, kozak_tp))
        for i in range(UP_START + 3, len(seq) - 2, 3):
            if (i - UP_START) < min_internal_nt: continue
            codon = seq[i:i+3]
            kozak_tn = seq[i-4:i+4]
            if codon in CODON: rows_tn.append((codon, tid, i, seq, 0, kozak_tn))
    return pd.DataFrame(rows_tp + rows_tn, columns=["codon", "tid", "tis_pos", "seq", "label", "kozak_seq"])
# =========================================================
# k-mer LM
# =========================================================
def count_motifs(up, motifs=MOTIFS_DNA):
    up = up.upper()
    out = {}
    for m in motifs:
        m = m.upper()
        Lm = len(m)
        cnt = 0
        if Lm == 5:
            up_5 = up[100:]
            for i in range(len(up_5) - Lm + 1):
                if up_5[i:i+Lm] == m:
                    cnt += 1
        else:
            for i in range(len(up) - Lm + 1):
                if up[i:i+Lm] == m:
                    cnt += 1
        out[m] = cnt
    return out
def cu_fraction(up):
    return (sum(1 for ch in up if ch in "CT") / len(up)) if up else float("nan")
def upstream_window(seq, tis_pos0):
    s = tis_pos0 - UP_START
    e = tis_pos0 - UP_END
    if s < 0 or e > len(seq) or e <= s: return None
    up = seq[s:e].upper()
    return up if len(up) == UP_LEN else None
def summarize_tp_tn(df, codon_efficiency):
    rows = []
    for _, r in df.iterrows():
        tis_pos = int(r["tis_pos"])
        kozak_seq = r['kozak_seq']
        up = upstream_window(r["seq"], tis_pos)
        cu = float("nan")
        motif_counts = {m: float("nan") for m in MOTIFS_DNA}
        if up is not None:
            up = up.upper()
            if DNA_RE_200.fullmatch(up):
                cu = cu_fraction(up)
                motif_counts = count_motifs(up)
        if DNA_RE.fullmatch(kozak_seq):
            eff = codon_efficiency.get(kozak_seq)
            codon_score = (eff if (eff is not None and not pd.isna(eff)) else np.nan)
        else:
            codon_score = np.nan
        row = {
            "codon": r["codon"],
            "tid": r["tid"],
            "tis_pos": tis_pos,
            "label": r["label"],
            "cu_fraction": cu,
            **motif_counts,
            "codon_score": codon_score
        }
        rows.append(row)
    df2 = pd.DataFrame(rows)
    return df2
# =========================================================
# CODON efficiency parsing
# =========================================================
def parse_codon_file(CODON_EFFICIENCY):
    df = pd.read_excel(CODON_EFFICIENCY)
    codon_efficiency = {}
    for _, r in df.iterrows():
        codon_efficiency[str(r['TIS Sequence'])] = r['TIS Efficiency']
    return codon_efficiency
# =========================================================
# Candidate scoring helpers
# =========================================================
def make_candidate_key(row):
    return f'{row["chrom"]}|{row["strand"]}|{int(row["phy_start"])}|{int(row["phy_end"])}|{str(row["codon"]).upper()}'
def score_candidates_excel():
    os.makedirs(OUT_DIR, exist_ok=True)
    # -----0) Analysis rules of TP/TN in Eu_genome -----
    genome_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(GENOME_FA, "fasta")}
    cds_dict = load_cds_fasta(CDS_FA, genome_dict)
    codon_efficiency = parse_codon_file(CODON_EFFICIENCY)
    df = sample_tp_tn_all_internal_atg(cds_dict)
    df2 = summarize_tp_tn(df, codon_efficiency)
    df2.to_excel(OUTPUT_TP_TN, index=False)
    # -----1) Analysis for candidate seq -----
    transcript_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(TRANSCRIPT_FA, "fasta")}
    df_sp = pd.read_csv(CANDIDATES_XLSX)
    required_cols = ["trans_id", "strand", "phy_start", "phy_end", "chrom", "codon"]
    missing = [c for c in required_cols if c not in df_sp.columns]
    if missing:
        raise ValueError(f"Candidate sheet 缺少必要列: {missing}")
    candidate_keys = []
    cu_fracs = []
    motif_cols = {m: [] for m in MOTIFS_DNA}
    codon_score = []
    for _, r in df_sp.iterrows():
        tran_id = r['trans_id']
        strand = str(r.get("strand", ""))
        start = int(r.get("phy_start", -1))
        end = int(r.get("phy_end", -1))
        trans_start = start + UP_LEN
        trans_end = end + UP_LEN
        candidate_keys.append(make_candidate_key(r))
        cu = np.nan
        mc = {m: np.nan for m in MOTIFS_DNA}
        g = transcript_dict[tran_id]
        if strand == "+":
            left = trans_start - (UP_START + 1)
            right = trans_start - 1
            up = g[left:right] if (left >= 0 and right <= len(g)) else None
        else:
            left = trans_end
            right = trans_end + UP_START
            up = str(Seq(g[left:right]).reverse_complement()) if (left >= 0 and right <= len(g)) else None
        if up is None or len(up) != UP_LEN or (not DNA_RE_200.fullmatch(up)):
            cu_fracs.append(np.nan); 
            continue
        cu = cu_fraction(up)
        mc = count_motifs(up)
        cu_fracs.append(cu)
        for m in MOTIFS_DNA:
            motif_cols[m].append(mc.get(m, np.nan))
        # codon and background
        w0 = str(r.get("kozak_seq", "")).upper()[2:10]
        if DNA_RE.fullmatch(w0):
            eff = codon_efficiency.get(w0)
            codon_score.append(eff if (eff is not None and not pd.isna(eff)) else np.nan)
        else:
            codon_score.append(np.nan)
    df_sp["candidate_key"] = candidate_keys
    df_sp["cu_fraction"] = cu_fracs
    for m in MOTIFS_DNA:
        df_sp[m] = motif_cols[m]
    df_sp["codon_efficiency"] = codon_score
    df_sp["codon_rank"] = (
        df_sp
        .groupby("trans_id")["codon_efficiency"]
        .rank(method="min", ascending=False)
        .astype("Int64")
    )
    df_sp.to_excel(OUT_XLSX, index=False)
if __name__ == "__main__":
    score_candidates_excel()