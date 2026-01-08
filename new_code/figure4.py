# # 转录组数据基础的NCP分布veen图
# import pandas as pd
# gene_expression_df = pd.read_csv(r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\03ouput\finally_filtered_genes.csv", index_col=0)
# tissue_mapping_df = pd.read_excel(r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\Total_rna_seq.xlsx", sheet_name="group")
# tissue_group_to_samples = (tissue_mapping_df.groupby(['Tissues', 'Group'])['Sample'].apply(list).to_dict())
# peptide_ids_by_tissue = {}
# for (tissue, group), samples in tissue_group_to_samples.items():
#     samples = [s for s in samples if s in gene_expression_df.columns]
#     group_expr = gene_expression_df.loc[:, samples]
#     group_mask = (group_expr >= 1).all(axis=1)
#     if tissue not in peptide_ids_by_tissue:
#         peptide_ids_by_tissue[tissue] = set()
#     peptide_ids_by_tissue[tissue].update(group_expr.index[group_mask].tolist())
# peptide_ids_by_tissue = {tissue: list(ids) for tissue, ids in peptide_ids_by_tissue.items()}
# max_length = max((len(ids) for ids in peptide_ids_by_tissue.values()), default=0)
# peptide_ids_df = pd.DataFrame({
#     tissue: ids + [None] * (max_length - len(ids))
#     for tissue, ids in peptide_ids_by_tissue.items()
# })
# output_excel_path = r"D:\Desktop\peptidemicro\00file\01figure\figure4\rnaseq_veen.xlsx"
# peptide_ids_df.to_excel(output_excel_path, index=False)

# # 数量-肽统计
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# from brokenaxes import brokenaxes
# import warnings
# warnings.filterwarnings('ignore')
# def plot_peptide_distribution_bar_broken_y(expression_matrix_path, output_pdf_path, output_csv_path, threshold=1, ylims=((0, 400), (600, 800), (2800, 3000)), bar_color="#4C72B0", edge_color="white"):
#     if expression_matrix_path.endswith('.csv'):
#         df = pd.read_csv(expression_matrix_path, index_col=0)
#     elif expression_matrix_path.endswith(('.xls', '.xlsx')):
#         df = pd.read_excel(expression_matrix_path, index_col=0)
#     else:
#         raise ValueError("仅支持 csv / xlsx")
#     print(f"读取矩阵：{df.shape[0]} 个ID × {df.shape[1]} 个列")
#     # ---------- 2. 二值化 & 统计 ----------
#     binary_df = (df >= threshold).astype(int)
#     counts_per_id = binary_df.sum(axis=1)
#     count_dist = counts_per_id.value_counts().sort_index()
#     total_ids = len(counts_per_id)
#     plot_data = pd.DataFrame({
#         "Number_of_Columns": count_dist.index.astype(int),
#         "Number_of_IDs": count_dist.values.astype(int),
#         "Percentage": (count_dist.values / total_ids * 100).round(2)
#     })
#     # ---------- 3) 输出 CSV ----------
#     plot_data.to_csv(output_csv_path, index=False)
#     print(f"✅ 统计结果已保存为：{output_csv_path}")
#     # ---------- 4) 画柱状图（broken y-axis） ----------
#     fig = plt.figure(figsize=(10, 8))
#     bax = brokenaxes(ylims=ylims, hspace=0.04, fig=fig)
#     x = plot_data["Number_of_Columns"].values
#     y = plot_data["Number_of_IDs"].values
#     bax.bar(x, y, width=0.8, color=bar_color, edgecolor=edge_color, linewidth=1.0)
#     bax.set_xlabel(f"Number of columns with CPM ≥ {threshold}")
#     bax.set_ylabel("Count of IDs")
#     bax.set_title("Peptide/ID Expression Distribution (Y-axis truncated)")
#     fig.savefig(output_pdf_path, bbox_inches="tight")
#     plt.close(fig)
#     print(f"✅ 截断柱状图已保存为：{output_pdf_path}")
#     return plot_data
# if __name__ == "__main__":
#     expression_file = r"D:\Desktop\peptidemicro\00file\01figure\finally_filtered_genes.csv"
#     output_pdf = r"D:\Desktop\peptidemicro\00file\01figure\figure4\peptide_distribution_bar_truncated.pdf"
#     output_csv = r"D:\Desktop\peptidemicro\00file\01figure\figure4\peptide_distribution_stat.csv"
#     plot_peptide_distribution_bar_broken_y(
#         expression_file,
#         output_pdf,
#         output_csv,
#         threshold=1,
#         ylims=((0, 400), (600, 800), (2800, 3000))
#     )


# # 起始密码子预测
# import pandas as pd
# from Bio import SeqIO
# from Bio.Seq import Seq
# from tqdm import tqdm
# from multiprocessing import Pool
# import warnings
# warnings.filterwarnings("ignore")
# # -----------------------------
# # 1) 常量与全局变量
# # -----------------------------
# START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
# STOP_CODONS  = {"TAA", "TAG", "TGA"}
# # 负链：沿用你原脚本的“在正向基因组序列上匹配的三联体集合”
# MINUS_START_CODONS = {"CAT", "CAG", "CAC", "CAA", "CGT"}
# MINUS_STOP_CODONS  = {"TTA", "CTA", "TCA"}
# CODON_PRIOR = {
#     "ATG": 0.0,      # 95%以上，基准
#     "CTG": -0.5,     # ~2%，相对常见
#     "GTG": -1.0,     # ~1%
#     "TTG": -1.5,     # <1%
#     "ACG": -3.0,     # 极罕见
#     "ATA": -4.0,     # 几乎不存在
#     "ATT": -4.0,     # 几乎不存在
#     "ATC": -4.0,     # 几乎不存在
#     # 其他罕见密码子：-4.0
# }
# MINUS_CODON_PRIOR = {
#     "CAT": 0.0,    # 对应ATG
#     "CAG": -0.5,   # 对应CTG
#     "CAC": -1.0,   # 对应GTG
#     "CAA": -1.5,   # 对应TTG
#     "CGT": -3.0,   # 对应ACG
#     "TAT": -4.0,   # 对应ATA
#     "AAT": -4.0,   # 对应ATT
#     "GAT": -4.0    # 对应ATC
# }
# _genome_dict = None
# _max_scan_nt = None
# def init_worker(genome_file, max_scan_nt):
#     global _genome_dict, _max_scan_nt
#     _max_scan_nt = max_scan_nt
#     _genome_dict = {}
#     for rec in SeqIO.parse(genome_file, "fasta"):
#         _genome_dict[rec.id] = rec.seq
# # -----------------------------
# # 2) 读取 peptide 表、按 accession 合并区间、算 coverage
# # -----------------------------
# def merge_segments(segments):
#     min_start = min(segments, key=lambda x: x[0])[0]
#     max_end = max(segments, key=lambda x: x[1])[1]
#     return min_start, max_end
# def filter_peptide_seq_cal_cov(peptide_file, database_file):
#     database_dict = {}
#     for rec in SeqIO.parse(database_file, "fasta"):
#         database_dict[rec.id] = rec.seq
#     df_NCP = pd.read_excel(peptide_file, sheet_name="unique")
#     accession_segments = {}
#     for _, row in df_NCP.iterrows():
#         accession = row['accessions']
#         start = int(row['start'])
#         end = int(row['end'])
#         chrom = row['chrom']
#         strand = row['strand']
#         accession_segments.setdefault(accession, []).append((start, end, chrom, strand))
#     stats = []
#     for accession, segments in accession_segments.items():
#         if accession not in database_dict:
#             continue
#         seq = database_dict[accession]
#         length_aa = len(seq)
#         min_start, max_end = merge_segments(segments)
#         cov_length_aa = (max_end - min_start + 1) / 3.0
#         coverage = cov_length_aa / length_aa if length_aa > 0 else None
#         chrom = segments[0][2]
#         strand = segments[0][3]
#         peptide_count = len(segments)
#         stats.append({
#             'accession': accession,
#             'chrom': chrom,
#             'strand': strand,
#             'min_start': min_start,
#             'max_end': max_end,
#             'coverage': coverage,
#             'peptide_count': peptide_count,
#             'sequence_length_aa': length_aa
#         })
#     return stats
# # -----------------------------
# # 3) Kozak 评分
# # -----------------------------
# def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank=6):
#     genome_seq = Seq(str(genome_seq).upper())
#     if strand == '+':
#         ctx = genome_seq[phy_start - 1 - flank: phy_start + flank + 2]
#     else:
#         ctx = genome_seq[phy_end - 3 - flank: phy_end + flank].reverse_complement()
#     codon_pos = flank
#     if len(ctx) < codon_pos + 4:
#         return {"core": None, "extended": None, "total": None, "context": None}
#     core = 0.0
#     if codon_pos - 3 >= 0 and ctx[codon_pos - 3] in ('A', 'G'):
#         if ctx[codon_pos - 3] == 'A':
#             core += 2.5
#         else:
#             core += 2.0
#     if codon_pos + 3 < len(ctx) and ctx[codon_pos + 3] == 'G':
#         core += 2.0
#     return {"core": core, "extended": None, "total": core, "context": str(ctx)}
# # -----------------------------
# # 4) 枚举候选 ORF（+）
# # -----------------------------
# def _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         s = min_start - 3 * i
#         if s - 1 < 0 or s + 2 > L:
#             continue
#         triplet = seq_str[s - 1:s + 2]
#         if triplet in START_CODONS:
#             yield s, triplet, CODON_PRIOR.get(triplet, -10.0)
# def _find_stop_plus(seq_str, max_end, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for j in range(max_steps):
#         e = max_end + 3 * j
#         if e + 3 > L:
#             break
#         triplet = seq_str[e:e + 3]
#         if triplet in STOP_CODONS:
#             return e
#     return None
# def enumerate_orf_candidates_plus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
#     seq_str = str(Seq(str(genome_seq)).upper())
#     candidates = []
#     for phy_start, triplet, prior_value in _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
#         phy_end = _find_stop_plus(seq_str, max_end, max_scan_nt)
#         if phy_end is None:
#             continue
#         if not (phy_start <= min_start and phy_end >= max_end):
#             continue
#         kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='+', flank=flank)
#         if kozak["total"] is None:
#             continue
#         total_score = prior_value * 0.6 + kozak["total"] * 0.4
#         candidates.append({
#             "phy_start": phy_start,
#             "phy_end": phy_end,
#             "prior": triplet,
#             "prior_value": prior_value,
#             "kozak_total": kozak["total"],
#             "kozak_seq": kozak["context"],
#             "total_score": total_score,
#             "start_to_peptide_nt": int(min_start - phy_start)
#         })
#     candidates.sort(key=lambda d: d["total_score"], reverse=True)
#     return candidates
# # -----------------------------
# # 5) 枚举候选 ORF（-）——沿用你原坐标/切片定义
# # -----------------------------
# def _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         left = max_end + 3 * (i - 1)
#         right = max_end + 3 * i
#         if right <= 0 or left >= L:
#             continue
#         triplet = seq_str[left:right]
#         if len(triplet) != 3:
#             continue
#         if triplet in MINUS_START_CODONS:
#             phy_end = max_end + 3 * i
#             prior_value = MINUS_CODON_PRIOR.get(triplet, -10.0)
#             yield phy_end, triplet, prior_value
# def _find_stop_minus(seq_str, min_start, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         left = (min_start - 1) - 3 * (i + 1)
#         right = (min_start - 1) - 3 * i
#         if right <= 0 or left >= L:
#             continue
#         triplet = seq_str[left:right]
#         if len(triplet) != 3:
#             continue
#         phy_start = min_start - 3 * i
#         if triplet in MINUS_STOP_CODONS:
#             return phy_start
#     return None
# def enumerate_orf_candidates_minus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
#     seq_str = str(Seq(str(genome_seq)).upper())
#     candidates = []
#     for phy_end, triplet_raw, prior_value in _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
#         phy_start = _find_stop_minus(seq_str, min_start, max_scan_nt)
#         if phy_start is None:
#             continue
#         if not (phy_start <= min_start and phy_end >= max_end):
#             continue
#         kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='-', flank=flank)
#         if kozak["total"] is None:
#             continue
#         total_score = prior_value * 0.6 + kozak["total"] * 0.4
#         triplet_rc = str(Seq(triplet_raw).reverse_complement())
#         candidates.append({
#             "phy_start": phy_start,
#             "phy_end": phy_end,
#             "prior": triplet_rc,
#             "prior_value": prior_value,
#             "kozak_total": kozak["total"],
#             "kozak_seq": kozak["context"],
#             "total_score": total_score,
#             "start_to_peptide_nt": int(phy_end - max_end)
#         })
#     candidates.sort(key=lambda d: d["total_score"], reverse=True)
#     return candidates
# # -----------------------------
# # 6) worker：返回 best + candidates
# # -----------------------------
# def run_scan_and_output_for_item(item):
#     global _genome_dict, _max_scan_nt
#     chrom = item['chrom']
#     strand = item['strand']
#     min_start = int(item['min_start'])
#     max_end = int(item['max_end'])
#     if chrom not in _genome_dict:
#         item['note'] = 'chrom_not_found'
#         item['best'] = None
#         item['candidates'] = []
#         return item
#     gseq = _genome_dict[chrom]
#     if strand == '+':
#         candidates = enumerate_orf_candidates_plus(gseq, min_start, max_end, _max_scan_nt, flank=6)
#     else:
#         candidates = enumerate_orf_candidates_minus(gseq, min_start, max_end, _max_scan_nt, flank=6)
#     item['best'] = candidates[0] if candidates else None
#     item['candidates'] = candidates
#     return item
# def run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=20):
#     with Pool(processes=nproc, initializer=init_worker, initargs=(genome_file, max_scan_nt)) as pool:
#         stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
#     return stats_update
# # -----------------------------
# # 7) main：输出两张表
# # -----------------------------
# def main():
#     peptide_file  = "finally_expressed_sp_info.xlsx"
#     database_file = "Eu_peptide_database_customized_5.fa"
#     genome_file   = "Eu_genome.fasta"
#     out_best = "output_best.csv"
#     out_cand = "output_candidates.csv"
#     # 1) 汇总 accession 覆盖信息
#     stats = filter_peptide_seq_cal_cov(peptide_file, database_file)
#     # 2) 多进程扫描并枚举候选
#     stats_update = run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=100)
#     # 3) 拆成 best 表与 candidates 表
#     best_rows = []
#     cand_rows = []
#     for it in stats_update:
#         base = {k: it[k] for k in it.keys() if k not in ("best", "candidates")}
#         # best：每个 accession 一行
#         if it.get("best") is not None:
#             best_rows.append({**base, **it["best"]})
#         else:
#             best_rows.append({
#                 **base,
#                 "phy_start": None, "phy_end": None, "prior": None,
#                 "prior_value": None, "kozak_total": None, "kozak_seq": None,
#                 "total_score": None, "start_to_peptide_nt": None
#             })
#         # candidates：每个 accession 多行（rank 从 1 开始）
#         for rank, cand in enumerate(it.get("candidates", []), start=1):
#             cand_rows.append({**base, "rank": rank, **cand})
#     df_best = pd.DataFrame(best_rows)
#     df_cand = pd.DataFrame(cand_rows)
#     df_best.to_csv(out_best, index=False)
#     df_cand.to_csv(out_cand, index=False)
#     print(f"✅ 已输出 best：{out_best}（每个 accession 最优 ORF）")
#     print(f"✅ 已输出 candidates：{out_cand}（每个 accession 所有候选 ORF，含 rank）")
# if __name__ == "__main__":
#     main()

# # 提取基因的Kozak序列
# import pandas as pd
# from Bio import SeqIO
# import re
# from tqdm import tqdm
# position_regex = re.compile(r'Position=([GWHBISF\d]+): (\d+).*?(\d+):\s([+-])')
# def extract_translation_start_positions(protein_file):
#     protein_start_positions = {}
#     for rec in SeqIO.parse(protein_file, "fasta"):
#         description = rec.description
#         match = position_regex.search(description)
#         if match:
#             chrom = match.group(1) 
#             start_position = match.group(2)
#             end_position = match.group(3)
#             strand = match.group(4)
#             protein_start_positions[rec.id] = (chrom, int(start_position), int(end_position), strand)
#     return protein_start_positions
# def extract_kozak_sequence_from_genome(genome_dict, chrom, start, end, strand, flank=6):
#     if chrom in genome_dict:
#         seq = genome_dict[chrom]
#         if strand == '-':
#             kozak_seq = seq[end-9:end+6].reverse_complement()
#         else:
#             kozak_seq = seq[start-7:start+8]
#     return kozak_seq
# def extract_kozak_sequences(protein_file, genome_file):
#     protein_positions = extract_translation_start_positions(protein_file)
#     kozak_sequences = {}
#     genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(genome_file, "fasta")}
#     for protein_id, (gene_id, start, end, strand) in tqdm(protein_positions.items(), desc="Extracting Kozak sequences", unit="protein"):
#         chrom = gene_id
#         kozak_seq = extract_kozak_sequence_from_genome(genome_dict, chrom, start, end, strand)
#         if kozak_seq:
#             kozak_sequences[protein_id] = kozak_seq
#         else:
#             print(f"警告: 未找到基因组序列 {protein_id} 对应的 Kozak 序列")
#     return kozak_sequences
# def main():
#     protein_file = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_peptide_database.fa"
#     genome_file = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_genome.fasta"
#     output_file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\gene_kozak_sequences.xlsx"
#     kozak_sequences = extract_kozak_sequences(protein_file, genome_file)
#     df_kozak = pd.DataFrame(list(kozak_sequences.items()), columns=['Protein ID', 'Kozak Sequence'])
#     df_kozak.to_excel(output_file, index=False)
# if __name__ == "__main__":
#     main()

# # Kozak序列位置分析
# import pandas as pd
# file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\output_best.xlsx"
# df = pd.read_excel(file, sheet_name="GTG")
# seqs = df['kozak_seq'].dropna().astype(str)
# seqs = seqs[seqs.str.len() == 15]
# pos_counts = {i: {} for i in range(1, 16)}
# for s in seqs:
#     for i, ch in enumerate(s, start=1):
#         pos_counts[i][ch] = pos_counts[i].get(ch, 0) + 1
# all_chars = sorted({ch for counts in pos_counts.values() for ch in counts.keys()})
# result_df = pd.DataFrame(
#     index=all_chars,
#     columns=[f"X{i}" for i in range(1, 16)]
# )
# for pos in range(1, 16):
#     for ch, cnt in pos_counts[pos].items():
#         result_df.loc[ch, f"X{pos}"] = cnt
# result_df = result_df.fillna(0).astype(int)
# result_df.to_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\sp_GTG_kozak_seq_statistic.xlsx")











import os, re, json, random, math
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

# =========================================================
# CONFIG
# =========================================================
CDS_FA = "data/CDS.fa"
GENOME_FA = "data/genome.fa"
CANDIDATES_XLSX = "data/candidates.xlsx"
GFF3_FA = "data/genome.gff3"
OUT_XLSX = "out/candidates_scored.xlsx"
SHEET_NAME = ""

TP_ONLY_KMER_SHUFFLES = 200
TP_ONLY_KMER_MIN_Z = 2.5
TP_ONLY_KMER_ARTIFACT = "out/tp_only_kmer_shuffle_weights.json"

CODON_BONUS_SCHEMES = {
    "weak":   {"ATG": 0.0, "CTG": -0.1, "GTG": -0.2, "TTG": -0.2, "ACG": -0.3},
    "medium": {"ATG": 0.0, "CTG": -0.5, "GTG": -1.0, "TTG": -1.0, "ACG": -2.0},
    "strong": {"ATG": 0.0, "CTG": -1.0, "GTG": -2.0, "TTG": -2.0, "ACG": -3.0},
}
DEFAULT_OTHER_CODON_BONUS = -3.0

# =========================================================
# CDS parsing
# =========================================================
def parse_gff_file(gff3_file):
    has_5utr = set()
    with open(gff3_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            if parts[2] != "five_prime_UTR": continue
            m = re.search(r"Parent_Accession=([^;]+)", parts[8])
            if m: has_5utr.add(m.group(1))
    return has_5utr
def parse_tid_details(description, cds_seq, genome_dict):
    chrom, exon_part, strand = description.split("\t")[3].split("=")[1].split(": ")
    exon_coords = [(int(s), int(e)) for s, e in (b.split("-") for b in exon_part.split(","))]
    g = genome_dict[chrom]
    if strand == "+":
        start = min(s for s, _ in exon_coords)
        if start - 101 < 0: return None
        up_seq = g[start - 101: start - 1] + cds_seq
    elif strand == "-":
        end = max(e for _, e in exon_coords)
        if end + 100 > len(g): return None
        up_seq = str(Seq(g[end: end + 100]).reverse_complement()) + cds_seq
    return up_seq
def load_cds_fasta(cds_file, gff3_file, genome_dict):
    has_5utr = parse_gff_file(gff3_file)
    cds_dict = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        tid = record.id
        if tid not in has_5utr: continue
        merged = parse_tid_details(record.description, str(record.seq), genome_dict)
        cds_dict[tid] = merged
    return cds_dict    
# =========================================================
# Training data: TP/TN sampling (only used for LR training on context, not for motif discovery)
# =========================================================
def upstream_window(seq, tis_pos0):
    s = tis_pos0 - 100; e = tis_pos0 - 0
    if s < 0 or e > len(seq) or e <= s: return None
    up = seq[s:e]
    return up if len(up) == 100 else None
def sample_tp_tn(cds_dict, tn_per_tx=5, min_internal_nt=30, seed=13):
    import random
    random.seed(seed)
    rows = []
    for tid, seq in cds_dict.items():
        if seq[100:103] != "ATG": continue
        rows.append((tid, 100, seq, 1))
        cand = []
        for i in range(103, len(seq) - 2, 3):
            if i - 100 < min_internal_nt: continue
            if seq[i:i+3] == "ATG": cand.append(i)
        if cand:
            picks = random.sample(cand, k=min(tn_per_tx, len(cand)))
            for i in picks:
                rows.append((tid, i, seq, 0))
    return pd.DataFrame(rows, columns=["tid", "tis_pos", "seq", "label"])
def compute_bg(cds_dict):
    c = {b: 0 for b in "ACGT"}
    for seq in cds_dict.values():
        for ch in seq[100:]:
            if ch in c: c[ch] += 1
    tot = sum(c.values())
    return {b: (c[b] / tot if tot else 0.25) for b in c}
def kozak_window(seq, tis_pos):
    s = tis_pos - 6
    e = tis_pos + 7
    if s < 0 or e > len(seq):
        return None
    w = seq[s:e]
    return w if len(w) == 13 else None
def build_pwm(tp_windows, pseudocount=1.0):
    if not tp_windows:
        raise ValueError("TP kozak windows 为空，无法建 PWM")
    L = len(tp_windows[0])
    counts = {b: np.full(L, pseudocount, dtype=float) for b in "ACGT"}
    for w in tp_windows:
        if len(w) != L: continue
        for i, ch in enumerate(w):
            if ch in counts:
                counts[ch][i] += 1.0
    totals = sum(counts[b] for b in "ACGT")
    return {b: (counts[b] / totals).tolist() for b in "ACGT"}
def kozak_logodds(w, pwm, bg, eps=1e-12):
    s = 0.0
    for i, ch in enumerate(w):
        if ch in "ACGT":
            ratio = pwm[ch][i] / max(bg[ch], eps)
            s += math.log2(max(ratio, eps))
    return s
def cu_fraction(up):
    return (sum(1 for ch in up if ch in "CT") / len(up)) if up else float("nan")
def tp_only_kmer_score(up, weights_z, k):
    if up is None or len(up) < k:
        return float("nan")
    up = up.upper()
    s = 0.0
    L = len(up)
    for i in range(L - k + 1):
        km = up[i:i+k]
        w = weights_z.get(km)
        if w is not None:
            s += w
    return float(s)
def train_pwm_and_weights(cds_dict, w3, w5):
    bg = compute_bg(cds_dict)
    df = sample_tp_tn(cds_dict)
    tp_ws = []
    for _, r in df[df["label"] == 1].iterrows():
        w = kozak_window(r["seq"], int(r["tis_pos"]))
        if w is not None:
            tp_ws.append(w)
    pwm = build_pwm(tp_ws)
    # LR features
    X = []; y = []
    for _, r in df.iterrows():
        seq = r["seq"]; tis = int(r["tis_pos"])
        w = kozak_window(seq, tis)
        up = upstream_window(seq, tis)
        if w is None or up is None: continue
        kz = kozak_logodds(w, pwm, bg)
        cu = cu_fraction(up)
        s3 = tp_only_kmer_score(up, w3, 3)
        s5 = tp_only_kmer_score(up, w5, 5)
        X.append([kz, cu, s3, s5])
        y.append(int(r["label"]))
    X = np.asarray(X, float)
    y = np.asarray(y, int)
    if X.shape[0] < 50: raise ValueError(f"训练样本过少：{X.shape[0]}")
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    lr = LogisticRegression(penalty="l2", solver="liblinear", max_iter=2000, class_weight="balanced")
    lr.fit(Z, y)
    weights = {
        "w1": float(lr.coef_[0][0]),  # kozak
        "w2": float(lr.coef_[0][1]),  # cu_fraction
        "w3": float(lr.coef_[0][2]),  # tp_only_3mer_score
        "w4": float(lr.coef_[0][3]),  # tp_only_5mer_score
        "beta0": float(lr.intercept_[0]),
    }
    std = {
        "mu": mu.tolist(),
        "sd": sd.tolist(),
        "feat": ["kozak", "cu_fraction", "tp_only_3mer_score", "tp_only_5mer_score"],
    }
    return pwm, bg, weights, std
# =========================================================
# Normalized background influence of TP
# =========================================================
def _eulerian_trail_from_adj(adj, start):
    local = {u: adj[u][:] for u in adj}
    for u in local: random.shuffle(local[u])
    stack = [start]
    path = []
    while stack:
        v = stack[-1]
        if local[v]: stack.append(local[v].pop())
        else: path.append(stack.pop())
    path.reverse()
    return path
def _build_dinucl_graph(seq):
    seq = seq.upper()
    adj = {b: [] for b in "ACGT"}
    for a, b in zip(seq[:-1], seq[1:]):
        if a in adj and b in "ACGT": adj[a].append(b)
        else: return None
    return adj
def dinuc_shuffle(seq, max_tries=20):
    seq = seq.upper()
    if len(seq) < 2 or any(ch not in "ACGT" for ch in seq): return None
    adj = _build_dinucl_graph(seq)
    if adj is None: return None
    start = seq[0]
    for _ in range(max_tries):
        path = _eulerian_trail_from_adj(adj, start)
        if len(path) == len(seq):
            return "".join(path)
    return None    
def count_kmers(seq, k):
    c = Counter()
    seq = seq.upper()
    L = len(seq)
    pat = re.compile(rf"[ACGT]{{{k}}}")
    for i in range(L - k + 1):
        s = seq[i:i+k]
        if pat.fullmatch(s):
            c[s] += 1
    return c
def learn_tp_only_kmer_weights(tp_ups, k, n_shuffle=200, min_z=2.5, seed=7, top_keep=None):
    random.seed(seed)
    obs = Counter()
    for up in tp_ups: obs.update(count_kmers(up, k))
    kmers = [''.join(p) for p in __import__("itertools").product("ACGT", repeat=k)]
    k2i = {km:i for i,km in enumerate(kmers)}
    rep_sums = np.zeros((n_shuffle, len(kmers)), dtype=np.float64)
    for r in range(n_shuffle):
        c = np.zeros(len(kmers), dtype=np.float64)
        for up in tp_ups:
            sh = dinuc_shuffle(up)
            if sh is None: continue
            L = len(sh)
            for i in range(L - k + 1):
                km = sh[i:i+k]
                j = k2i.get(km)
                if j is not None: c[j] += 1.0
        rep_sums[r] = c
    exp = rep_sums.mean(axis=0)
    sd = rep_sums.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    obs_arr = np.zeros(len(kmers), dtype=np.float64)
    for km, i in k2i.items():
        obs_arr[i] = float(obs.get(km, 0))
    z = (obs_arr - exp) / sd
    weights_z = {km: float(z[k2i[km]]) for km in kmers if z[k2i[km]] >= min_z}
    if top_keep is not None and len(weights_z) > top_keep:
        weights_z = dict(sorted(weights_z.items(), key=lambda x: x[1], reverse=True)[:top_keep])
    full = {
        "params": {"k": k, "n_shuffle": n_shuffle, "min_z": min_z, "seed": seed,
                   "window": f"[-{100}, -{0})", "up_len": 100, "top_keep": top_keep},
        "weights_z": weights_z,
    }
    return full
def ensure_tp_only_artifact(cds_dict):
    os.makedirs(os.path.dirname(TP_ONLY_KMER_ARTIFACT), exist_ok=True)
    if os.path.exists(TP_ONLY_KMER_ARTIFACT):
        with open(TP_ONLY_KMER_ARTIFACT, "r", encoding="utf-8") as f:
            return json.load(f)
    tp_ups = []
    for _, seq in cds_dict.items():
        if seq[100:103] != "ATG": continue
        up = upstream_window(seq, 100)
        if up is not None and re.fullmatch(r"[ACGT]{100}", up): tp_ups.append(up)
    if len(tp_ups) < 50: raise ValueError(f"TP upstream windows 太少（{len(tp_ups)}），无法建 TP-only kmer 背景模型")
    art3 = learn_tp_only_kmer_weights(tp_ups, k=3, n_shuffle=TP_ONLY_KMER_SHUFFLES, min_z=TP_ONLY_KMER_MIN_Z, seed=7, top_keep=None)
    art5 = learn_tp_only_kmer_weights(tp_ups, k=5, n_shuffle=TP_ONLY_KMER_SHUFFLES, min_z=TP_ONLY_KMER_MIN_Z, seed=11, top_keep=500)
    obj = {
        "params": {"window": f"[-{100}, -{0})", "up_len": 100},
        "3mer": art3,
        "5mer": art5,
    }
    with open(TP_ONLY_KMER_ARTIFACT, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
    return obj                  
# =========================================================
# Candidate scoring helpers
# =========================================================
def codon_bonus_by_scheme(codon, scheme_dict):
    c = str(codon).upper()
    return float(scheme_dict.get(c, DEFAULT_OTHER_CODON_BONUS))
def make_candidate_key(row):
    return f'{row["chr"]}|{row["strand"]}|{int(row["phy_start"])}|{int(row["phy_end"])}|{str(row["codon"]).upper()}'
# =========================================================
# Main scoring pipeline
# =========================================================
def score_candidates_excel():
    os.makedirs(os.path.dirname(OUT_XLSX), exist_ok=True)
    genome_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(GENOME_FA, "fasta")}
    cds_dict = load_cds_fasta(CDS_FA, GFF3_FA, genome_dict)
    art = ensure_tp_only_artifact(cds_dict)
    w3 = art["3mer"]["weights_z"]
    w5 = art["5mer"]["weights_z"]
    pwm, bg, weights, std = train_pwm_and_weights(cds_dict, w3, w5)
    df_sp = pd.read_excel(CANDIDATES_XLSX, sheet_name=SHEET_NAME)
    mu = np.array(std["mu"], float)
    sd = np.array(std["sd"], float)
    candidate_keys = []
    kozak_scores = []
    cu_fracs = []
    tp3_scores = []
    tp5_scores = []
    tis_scores = {name: [] for name in CODON_BONUS_SCHEMES}   
    for _, r in df_sp.iterrows():
        chrom = str(r["chr"])
        strand = str(r["strand"])
        start = int(r["phy_start"])
        end = int(r["phy_end"])
        codon = str(r.get("codon", "")).upper()
        candidate_keys.append(make_candidate_key(r))
        w = str(r.get("kozak", "")).upper()[:13]
        g = genome_dict[chrom]
        if strand == "+":
            left = start - 101; right = start - 1
            up = g[left:right] if (left >= 0 and right <= len(g)) else None
        else:
            left = end; right = end + 100
            up = str(Seq(g[left:right]).reverse_complement()) if (left >= 0 and right <= len(g)) else None
        kz = kozak_logodds(w, pwm, bg)
        cu = cu_fraction(up)
        s3 = tp_only_kmer_score(up, w3, 3)
        s5 = tp_only_kmer_score(up, w5, 5)
        z = (np.array([kz, cu, s3, s5], float) - mu) / sd
        # LR logit score
        score_base = (
            weights["beta0"]
            + weights["w1"] * z[0]
            + weights["w2"] * z[1]
            + weights["w3"] * z[2]
            + weights["w4"] * z[3]
        )
        kozak_scores.append(kz)
        cu_fracs.append(cu)
        tp3_scores.append(s3)
        tp5_scores.append(s5)
        for scheme_name, scheme_dict in CODON_BONUS_SCHEMES.items():
            tis_scores[scheme_name].append(score_base + codon_bonus_by_scheme(codon, scheme_dict))
    df_sp["candidate_key"] = candidate_keys
    df_sp["kozak_score"] = kozak_scores
    df_sp["cu_fraction"] = cu_fracs
    df_sp["tp_only_3mer_score"] = tp3_scores
    df_sp["tp_only_5mer_score"] = tp5_scores
    for scheme_name in CODON_BONUS_SCHEMES:
        score_col = f"tis_score_{scheme_name}"
        rank_col = f"rank_{scheme_name}"
        df_sp[score_col] = tis_scores[scheme_name]
        df_sp[rank_col] = df_sp.groupby("peptide_id")[score_col].rank(ascending=False, method="first")
        top3 = df_sp[df_sp[rank_col] <= 3].sort_values(["peptide_id", rank_col])
        top3.to_excel(OUT_XLSX.replace(".xlsx", f"_top3_{scheme_name}.xlsx"), index=False)
    df_sp.to_excel(OUT_XLSX, index=False)
    top1 = {}
    for scheme_name in CODON_BONUS_SCHEMES:
        score_col = f"tis_score_{scheme_name}"
        idx = df_sp.groupby("peptide_id")[score_col].idxmax()
        idx = idx.dropna().astype(int)
        tmp = df_sp.loc[idx, ["peptide_id", "candidate_key"]].set_index("peptide_id")["candidate_key"]
        top1[scheme_name] = tmp
    summary = pd.DataFrame({
        "top1_weak": top1.get("weak"),
        "top1_medium": top1.get("medium"),
        "top1_strong": top1.get("strong"),
    })
    summary["stable_top1"] = (summary["top1_weak"] == summary["top1_medium"]) & (summary["top1_medium"] == summary["top1_strong"])
    summary = summary.reset_index()
    summary.to_excel(OUT_XLSX.replace(".xlsx", "_sensitivity_summary.xlsx"), index=False)
    out_dir = os.path.dirname(OUT_XLSX)
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "pwm.json"), "w", encoding="utf-8") as f:
        json.dump(pwm, f, indent=2)
    with open(os.path.join(out_dir, "weights.json"), "w", encoding="utf-8") as f:
        json.dump(weights, f, indent=2)
    with open(os.path.join(out_dir, "standardizer.json"), "w", encoding="utf-8") as f:
        json.dump(std, f, indent=2)
    with open(os.path.join(out_dir, "tp_only_3mer5mer_shuffle_weights.json"), "w", encoding="utf-8") as f:
        json.dump(art, f, indent=2)                         
if __name__ == "__main__":
    score_candidates_excel()