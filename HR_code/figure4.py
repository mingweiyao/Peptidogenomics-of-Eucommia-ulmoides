# # figure4b：不同组织中肽表达数量
# import pandas as pd
# gene_expression_df = pd.read_csv(r"F:\Eu_peptido\new_prepare\sp_expression.csv", index_col=0)
# tissue_mapping_df = pd.read_excel(r"F:\Eu_peptido\00file\00raw\rnaseq\Total_rna_seq.xlsx", sheet_name="Sheet4")
# tissue_to_samples = tissue_mapping_df.groupby('Tissues')['Sample'].apply(list).to_dict()
# peptide_ids_by_tissue = {}
# for tissue, samples in tissue_to_samples.items():
#     tissue_expr = gene_expression_df[samples]
#     expressed_peptides = tissue_expr[(tissue_expr > 10).any(axis=1)].index.tolist()
#     peptide_ids_by_tissue[tissue] = expressed_peptides
# max_length = max(len(ids) for ids in peptide_ids_by_tissue.values())
# peptide_ids_df = pd.DataFrame({
#     tissue: ids + [None]*(max_length - len(ids)) 
#     for tissue, ids in peptide_ids_by_tissue.items()
# })
# output_csv_path = r"F:\Eu_peptido\new_prepare\figure4_peptide_expression_ids.csv"
# peptide_ids_df.to_csv(output_csv_path, index=False)
# print("✅ 结果已保存至:", output_csv_path)
# print("\n示例数据：")
# print(peptide_ids_df.head())

# # figure1:数量-肽统计
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.backends.backend_pdf import PdfPages
# import warnings
# warnings.filterwarnings('ignore')
# def plot_peptide_distribution(expression_matrix_path, output_pdf_path, threshold=10):
#     sns.set_theme(style="whitegrid")
#     plt.rcParams.update({
#         'font.size': 12,
#         'font.family': 'DejaVu Sans',
#         'axes.titlepad': 20,
#         'axes.labelpad': 10,
#         'figure.dpi': 300
#     })
#     try:
#         if expression_matrix_path.endswith('.csv'):
#             df = pd.read_csv(expression_matrix_path, index_col=0)
#         elif expression_matrix_path.endswith('.txt'):
#             df = pd.read_csv(expression_matrix_path, sep='\t', index_col=0)
#         elif expression_matrix_path.endswith(('.xls', '.xlsx')):
#             df = pd.read_excel(expression_matrix_path, index_col=0)
#         else:
#             raise ValueError("Unsupported file format. Please use .csv, .txt, .xls, or .xlsx")
#         print(f"成功读取数据: {df.shape[0]}个肽段 × {df.shape[1]}个样本")        
#     except Exception as e:
#         print(f"文件读取错误: {e}")
#         return None
#     binary_df = (df >= threshold).astype(int)
#     peptide_counts = binary_df.sum(axis=1)
#     total_peptides = len(peptide_counts)
#     count_dist = peptide_counts.value_counts().sort_index()
#     plot_data = pd.DataFrame({
#         'Number_of_Samples': count_dist.index,
#         'Number_of_Peptides': count_dist.values,
#         'Percentage': (count_dist.values / total_peptides * 100).round(2)
#     })
#     colors = ["#1f78b4", "#66c2a5", "#fc8d62", "#e31a1c"]
#     cmap = LinearSegmentedColormap.from_list("custom_diverging", colors)
#     norm = plt.Normalize(plot_data['Number_of_Samples'].min(), 
#                          plot_data['Number_of_Samples'].max())
#     with PdfPages(output_pdf_path) as pdf:
#         fig1, ax1 = plt.subplots(figsize=(12, 7))
#         bars = ax1.bar(plot_data['Number_of_Samples'], 
#                        plot_data['Number_of_Peptides'], 
#                        width=0.8,
#                        edgecolor='white', 
#                        linewidth=1.0,
#                        color=[cmap(norm(val)) for val in plot_data['Number_of_Samples']],
#                        alpha=0.8)
#         for i, (x, y) in enumerate(zip(plot_data['Number_of_Samples'], 
#                                        plot_data['Number_of_Peptides'])):
#             if y > max(plot_data['Number_of_Peptides']) * 0.05:
#                 ax1.text(x, y + max(plot_data['Number_of_Peptides']) * 0.01, 
#                         f'{y:,}', ha='center', va='bottom', fontsize=9)
#         x = plot_data['Number_of_Samples'].values
#         y = plot_data['Number_of_Peptides'].values
#         if len(x) > 2:
#             try:
#                 coeffs = np.polyfit(x, y, 2)
#                 poly = np.poly1d(coeffs)
#                 trend_x = np.linspace(min(x), max(x), 100)
#                 trend_y = poly(trend_x)
#                 trend_y[trend_y < 0] = 0
#                 ax1.plot(trend_x, trend_y, 
#                          color='black', 
#                          linestyle='--', 
#                          linewidth=2,
#                          label='Trend line (quadratic fit)')
#                 ax1.legend(loc='best')
#             except:
#                 print("趋势线拟合失败，可能数据点不足")
#         ax1.set_title('Distribution of Peptides by Number of Expressing Samples\n'
#                      f'(Expression Threshold: {threshold})', 
#                      fontweight='bold', fontsize=14)
#         ax1.set_xlabel('Number of Samples Where Expressed', fontweight='bold')
#         ax1.set_ylabel('Number of Peptides', fontweight='bold')
#         ax1.set_xticks(plot_data['Number_of_Samples'])
#         ax1.set_xticklabels(plot_data['Number_of_Samples'].astype(int))
#         sns.despine()
#         plt.tight_layout()
#         pdf.savefig(fig1, bbox_inches='tight')
#         plt.close(fig1)
#         fig2, ax2 = plt.subplots(figsize=(12, 7))
#         bars = ax2.bar(plot_data['Number_of_Samples'], 
#                        plot_data['Percentage'],
#                        width=0.8,
#                        edgecolor='white', 
#                        linewidth=1.0,
#                        color=[cmap(norm(val)) for val in plot_data['Number_of_Samples']],
#                        alpha=0.8)
#         for i, (x, y) in enumerate(zip(plot_data['Number_of_Samples'], 
#                                        plot_data['Percentage'])):
#             if y > 1.0:
#                 ax2.text(x, y + 1, 
#                         f'{y:.1f}%', ha='center', va='bottom', fontsize=9)
#         ax2.set_title('Percentage of Peptides by Number of Expressing Samples\n'
#                      f'(Expression Threshold: {threshold})', 
#                      fontweight='bold', fontsize=14)
#         ax2.set_xlabel('Number of Samples Where Expressed', fontweight='bold')
#         ax2.set_ylabel('Percentage of Total Peptides (%)', fontweight='bold')
#         ax2.set_xticks(plot_data['Number_of_Samples'])
#         ax2.set_xticklabels(plot_data['Number_of_Samples'].astype(int))
#         ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
#         ax2.grid(True, axis='y', alpha=0.3, linestyle='--')
#         sns.despine()
#         plt.tight_layout()
#         pdf.savefig(fig2, bbox_inches='tight')
#         plt.close(fig2)
#         fig3, ax3 = plt.subplots(figsize=(12, 7))
#         cumulative_percentage = np.cumsum(plot_data['Percentage'])
#         ax3.plot(plot_data['Number_of_Samples'], cumulative_percentage, 
#                 marker='o', markersize=6, linewidth=3, color='#2b8cbe')
#         ax3.fill_between(plot_data['Number_of_Samples'], cumulative_percentage, 
#                         alpha=0.2, color='#a6bddb')
#         for i, (x, y) in enumerate(zip(plot_data['Number_of_Samples'], cumulative_percentage)):
#             if i % max(1, len(plot_data) // 5) == 0:
#                 ax3.annotate(f'{y:.1f}%', xy=(x, y), xytext=(5, 5), 
#                            textcoords='offset points', fontweight='bold')        
#         ax3.set_title('Cumulative Distribution of Peptides\n'
#                      f'(Expression Threshold: {threshold})', 
#                      fontweight='bold', fontsize=14)
#         ax3.set_xlabel('Number of Samples Where Expressed', fontweight='bold')
#         ax3.set_ylabel('Cumulative Percentage (%)', fontweight='bold')
#         ax3.grid(True, alpha=0.3)
#         ax3.set_ylim(0, 105)        
#         sns.despine()
#         plt.tight_layout()
#         pdf.savefig(fig3, bbox_inches='tight')
#         plt.close(fig3)
#     print(f"\n{'='*60}")
#     print("肽段表达分布统计摘要")
#     print(f"{'='*60}")
#     print(f"总肽段数量: {total_peptides:,}")
#     print(f"表达阈值: {threshold}")
#     print(f"组织特异性肽段 (仅在1个样本中表达): {plot_data.iloc[0]['Number_of_Peptides']:,} "
#           f"({plot_data.iloc[0]['Percentage']:.1f}%)")
#     half_samples = binary_df.shape[1] // 2
#     widespread_peptides = sum(peptide_counts >= half_samples)
#     widespread_percentage = (widespread_peptides / total_peptides * 100)
#     print(f"广泛表达肽段 (在{half_samples}+个样本中表达): {widespread_peptides:,} "
#           f"({widespread_percentage:.1f}%)")
#     print(f"\nPDF文件已保存至: {output_pdf_path}")
#     return plot_data
# if __name__ == "__main__":
#     expression_file = r"F:\Eu_peptido\new_prepare\sp_expression.csv"
#     output_pdf = r"F:\Eu_peptido\new_prepare\peptide_expression_distribution.pdf"
#     results = plot_peptide_distribution(expression_file, output_pdf, threshold=10)
#     if results is not None:
#         print(f"\n{'='*60}")
#         print("详细分布统计结果:")
#         print(f"{'='*60}")
#         print(results.to_string(index=False, formatters={
#             'Number_of_Peptides': '{:,}'.format,
#             'Percentage': '{:.1f}%'.format
#         }))

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
