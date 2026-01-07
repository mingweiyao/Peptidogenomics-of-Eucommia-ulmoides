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



















import os, json, re, math
import numpy as np
import pandas as pd
from pyfaidx import Fasta
from sklearn.linear_model import LogisticRegression
from Bio import SeqIO

# ===========
# CONFIG
# ===========
CDS_FA = "data/CDS.fa"
GENOME_FA = "data/genome.fa"
CANDIDATES_XLSX = "data/candidates.xlsx"
OUT_XLSX = "out/candidates_scored.xlsx"
GFF3_FA = "data/genome.gff3"

# windows (DNA, genome-based)
KOZAK_UP = 6
KOZAK_DOWN = 4     # downstream after ATG (nt)
UP_START = 120
UP_END = 20

MOTIFS = ["TCTTC", "TCTCT"]  # UCUUC / UCUCU
MOTIF_MODE = "binary"        # or "count"

STOP = {"TAA","TAG","TGA"}
DNA_COMP = str.maketrans("ACGTacgtNn","TGCAtgcaNn")

def revcomp(s): return s.translate(DNA_COMP)[::-1]
def safe_log2(x, eps=1e-12): return math.log2(max(x, eps))

def parse_gff3_for_5utr(gff3_path):
    has_5utr = set()
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attributes = parts[8]
            if feature_type == "five_prime_UTR":
                match = re.search(r'Parent_Accession=([^;]+)', attributes)
                tid = match.group(1)
                has_5utr.add(tid)
    return has_5utr
def parse_tid_details(description):
    m = re.search(r"Position=([^ ]+)", description)
    pos_str = m.group(1)
    _, exon_part, strand = pos_str.split(":", 2)
    exon_coords = []
    for block in exon_part.split(","):
        s, e = block.split("-")
        exon_coords.append((int(s), int(e)))
    if strand == "+":
        A_genome = min(s for s, _ in exon_coords)
        start, end = A_genome, A_genome + 2
    else:
        A_genome = max(e for _, e in exon_coords)
        start, end = A_genome - 2, A_genome
    # 注意，这里的start/end仅代表起始密码子的基因组位置，而不代表真实的翻译方向
    return start, end, strand
def load_cds_fasta(path, gff3_path):
    has_5utr = parse_gff3_for_5utr(gff3_path)
    d = {}
    for record in SeqIO.parse(path, "fasta"):
        tid = record.id
        if tid not in has_5utr:
            continue
        codon_start, codon_end, strand = parse_tid_details(record.description)
        seq = str(record.seq).upper()
        d[tid] = {"seq": seq, "codon_start": codon_start, "codon_end": codon_end, "strand": strand}
    return d



def upstream_window(seq, tis_pos0):
    s = tis_pos0 - UP_START
    e = tis_pos0 - UP_END
    if s < 0 or e > len(seq) or e <= s: return None
    return seq[s:e]

def build_pwm(tp_windows, pseudocount=1.0):
    L = len(tp_windows[0])
    counts = {b: np.full(L, pseudocount) for b in "ACGT"}
    for w in tp_windows:
        for i,ch in enumerate(w):
            if ch in counts: counts[ch][i]+=1
    totals = sum(counts[b] for b in "ACGT")
    return {b:(counts[b]/totals).tolist() for b in "ACGT"}

def kozak_logodds(w, pwm, bg):
    s=0.0
    for i,ch in enumerate(w):
        if ch in "ACGT":
            s += safe_log2(pwm[ch][i] / max(bg[ch], 1e-12))
    return s

def cu_fraction(up):
    return (sum(1 for ch in up if ch in "CT")/len(up)) if up else float("nan")

def cu_motif_score(up):
    if MOTIF_MODE=="binary":
        return float(sum(1 for m in MOTIFS if m in up))
    tot=0
    for m in MOTIFS:
        tot += len(re.findall(f"(?={m})", up))
    return float(tot)


def compute_bg(cds_dict):
    c = {b:0 for b in "ACGT"}; tot=0
    for rec in cds_dict.values():
        for ch in rec["seq"]:
            if ch in c:
                c[ch] += 1
                tot += 1
    return {b:(c[b]/tot if tot else 0.25) for b in "ACGT"}
def sample_tp_tn(cds_dict, tn_per_tx=5, min_internal_nt=30, seed=13):
    import random
    random.seed(seed)
    rows=[]
    for tid, values in cds_dict.items():
        s = values[0]
        if len(s) < 200: # 谨慎判断是否可行
            continue
        if s[:3]=="ATG":
            rows.append((tid,0,s,1))
        cand=[]
        for i in range(3, len(s)-2, 3):
            if i < min_internal_nt: 
                continue
            if s[i:i+3]=="ATG":
                cand.append(i)
        if cand:
            picks=random.sample(cand, k=min(tn_per_tx, len(cand)))
            for i in picks:
                rows.append((tid,i,s,0))
    return pd.DataFrame(rows, columns=["tid","tis_pos","seq","label"])
def kozak_window(seq, tis_pos0):
    s = tis_pos0 - KOZAK_UP
    e = tis_pos0 + 3 + KOZAK_DOWN
    if s < 0 or e > len(seq): return None
    return seq[s:e]
def train_pwm_and_weights(cds_dict, genome):
    bg = compute_bg(cds_dict)
    df = sample_tp_tn(cds_dict)

    tp_ws=[]
    for _,r in df[df["label"]==1].iterrows():
        w = kozak_window(r["seq"], int(r["tis_pos"]))
        if w is not None:
            tp_ws.append(w)
    pwm = build_pwm(tp_ws)

    # features for LR
    X=[]; y=[]
    for _,r in df.iterrows():
        seq=r["seq"]; tis=int(r["tis_pos"])
        w=kozak_window(seq,tis)
        up=upstream_window(seq,tis)
        if w is None or up is None: 
            continue
        X.append([kozak_logodds(w,pwm,bg), cu_fraction(up), cu_motif_score(up)])
        y.append(int(r["label"]))
    X=np.asarray(X,float); y=np.asarray(y,int)

    mu=X.mean(axis=0); sd=X.std(axis=0); sd[sd==0]=1.0
    Z=(X-mu)/sd

    lr=LogisticRegression(penalty="l2", solver="liblinear", max_iter=2000)
    lr.fit(Z,y)

    weights={
        "w1": float(lr.coef_[0][0]),
        "w2": float(lr.coef_[0][1]),
        "w3": float(lr.coef_[0][2]),
        "beta0": float(lr.intercept_[0]),
    }
    std={"mu": mu.tolist(), "sd": sd.tolist(), "feat":["kozak","cu_fraction","cu_motif"]}
    return pwm, bg, weights, std

def codon_bonus(codon):
    codon=codon.upper()
    if codon=="ATG": return 0.0
    if codon=="CTG": return -0.1
    if codon in ("ACG","GTG"): return -0.2
    return -0.3

def fetch_genome_seq(genome, chrom, start1, end1, strand):
    # 1-based inclusive -> python slice (0-based, end exclusive)
    seq = str(genome[chrom][start1-1:end1]).upper()
    if strand=="-":
        seq = revcomp(seq)
    return seq

def score_candidates_excel():
    os.makedirs(os.path.dirname(OUT_XLSX), exist_ok=True)
    genome = Fasta(GENOME_FA, as_raw=True, sequence_always_upper=True)
    cds = load_cds_fasta(CDS_FA, GFF3_FA)
    pwm, bg, weights, std = train_pwm_and_weights(cds, genome)
    df = pd.read_excel(CANDIDATES_XLSX)

    # 你需要确保这些列存在：peptide_id, chr, strand, start, end
    need={"peptide_id","chr","strand","start","end"}
    miss=need - set(df.columns)
    if miss:
        raise ValueError(f"candidates.xlsx 缺列: {sorted(miss)}")

    # 取出 codon（如果 excel 没给，就从 genome 取 start..end 的3nt）
    if "codon" not in df.columns:
        df["codon"] = df.apply(lambda r: fetch_genome_seq(genome, r["chr"], int(r["start"]), int(r["end"]), r["strand"]), axis=1)

    mu=np.array(std["mu"],float); sd=np.array(std["sd"],float)

    # 对每个候选位点：从 genome 取窗口（注意：这是 genome-based 窗口）
    kozak_scores=[]; cu_fracs=[]; cu_motifs=[]; tis_scores=[]
    for _,r in df.iterrows():
        chrom=str(r["chr"]); strand=str(r["strand"])
        start1=int(r["start"]); end1=int(r["end"])
        # 取一个更大的片段，以便切 Kozak/upstream
        # 我们用 start1 作为 codon 的第1个碱基位置（A/T/C/G），对 '+' strand A 位点是 start1 (0-based => start1-1)
        # 对 '-' strand 我们已经 revcomp 了，所以仍把 A 位点当成该3nt的第1位
        big_left = start1 - UP_START - 10
        big_right = end1 + KOZAK_DOWN + 10
        if big_left < 1:
            kozak_scores.append(np.nan); cu_fracs.append(np.nan); cu_motifs.append(np.nan); tis_scores.append(np.nan)
            continue
        big_seq = fetch_genome_seq(genome, chrom, big_left, big_right, strand)

        # 在 big_seq 里，A 位点坐标（0-based）
        # 对 '+'：start1 对应 big_seq 的位置 = (start1 - big_left)
        # 对 '-'：fetch_genome_seq 已 revcomp，所以同样成立
        tis_pos0 = start1 - big_left  # position of first base of codon in big_seq
        # 这里我们把 tis_pos0 作为“start codon 的第1个碱基位置”，Kozak_window 实现是围绕 A 位点
        # 更严格可把 A 位置固定为 codon 第1位（对 ATG 这就是 A）
        w = kozak_window(big_seq, tis_pos0)
        up = upstream_window(big_seq, tis_pos0)
        if w is None or up is None:
            kozak_scores.append(np.nan); cu_fracs.append(np.nan); cu_motifs.append(np.nan); tis_scores.append(np.nan)
            continue

        kz = kozak_logodds(w, pwm, bg)
        cu = cu_fraction(up)
        cm = cu_motif_score(up)
        z = (np.array([kz,cu,cm],float) - mu) / sd
        score = weights["w1"]*z[0] + weights["w2"]*z[1] + weights["w3"]*z[2] + codon_bonus(r["codon"])

        kozak_scores.append(kz); cu_fracs.append(cu); cu_motifs.append(cm); tis_scores.append(score)

    df["kozak_score"]=kozak_scores
    df["cu_fraction"]=cu_fracs
    df["cu_motif_score"]=cu_motifs
    df["tis_score"]=tis_scores

    # 组内排序：同一个 peptide_id 选 top1/top3
    df["rank_in_peptide"] = df.groupby("peptide_id")["tis_score"].rank(ascending=False, method="first")

    df.to_excel(OUT_XLSX, index=False)

    # 另外输出 top3
    top3 = df[df["rank_in_peptide"]<=3].sort_values(["peptide_id","rank_in_peptide"])
    top3.to_excel(OUT_XLSX.replace(".xlsx","_top3.xlsx"), index=False)

    # 保存模型产物
    with open(os.path.join(os.path.dirname(OUT_XLSX), "pwm.json"), "w") as f:
        json.dump(pwm,f,indent=2)
    with open(os.path.join(os.path.dirname(OUT_XLSX), "weights.json"), "w") as f:
        json.dump(weights,f,indent=2)
    with open(os.path.join(os.path.dirname(OUT_XLSX), "standardizer.json"), "w") as f:
        json.dump(std,f,indent=2)

if __name__ == "__main__":
    score_candidates_excel()
