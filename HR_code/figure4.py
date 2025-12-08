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
# START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
# STOP_CODONS  = {"TAA", "TAG", "TGA"}
# MINUS_START_CODONS = {"CAT","CAG","CAC","CAA","CGT"}
# MINUS_STOP_CODONS = {"TTA","CTA","TCA"}
# CODON_PRIOR = {
#     "ATG": 0.0,

#     "CTG": -1.0,
#     "GTG": -1.0,
#     "TTG": -1.0,

#     "ACG": -1.5,
#     "AUA": -1.5,
#     "AUU": -1.5,
#     "AUC": -1.5,

#     "AAG": -2.0,
#     "AGG": -2.0,
#     "CGU": -2.0,
#     "CGC": -2.0,
#     "CGG": -2.0,
#     "CAG": -2.0
# }
# MINUS_CODON_PRIOR = {
#     "CAT": 0.0,

#     "CAG": -1.0,
#     "CAC": -1.0,
#     "CAA": -1.0,

#     "CGT": -1.5,
#     "TAT": -1.5,
#     "AAT": -1.5,
#     "GAT": -1.5,

#     "CTT": -2.0,
#     "CCT": -2.0,
#     "ACG": -2.0,
#     "GCG": -2.0,
#     "CCG": -2.0,
#     "CTG": -2.0
# }
# _genome_dict = None
# _max_scan_nt = None
# def init_worker(genome_file, max_scan_nt):
#     global _genome_dict, _max_scan_nt
#     _max_scan_nt = max_scan_nt
#     _genome_dict = {}
#     for rec in SeqIO.parse(genome_file, "fasta"):
#         _genome_dict[rec.id] = rec.seq
# def merge_segments(segments):
#     min_start = min(segments, key=lambda x: x[0])[0]
#     max_end = max(segments, key=lambda x: x[1])[1]  
#     return min_start, max_end
# def filter_peptide_seq_cal_cov(peptide_file, database_file):
#     database_dict = {}
#     for rec in SeqIO.parse(database_file, "fasta"):
#         database_dict[rec.id] = rec.seq
#     df_NCP = pd.read_excel(peptide_file)
#     df_NCP['proteins'] = df_NCP['proteins'].astype(int)
#     df_NCP_filter = df_NCP[(df_NCP['type'] != 'exon_diff') & (df_NCP['proteins'] == 1)]
#     accession_segments = {}
#     for _, row in df_NCP_filter.iterrows():
#         accession = row['accessions']
#         start = row['start']
#         end = row['end']
#         chrom = row['chrom']
#         strand = row['strand']
#         if accession not in accession_segments:
#             accession_segments[accession] = []
#         accession_segments[accession].append((start, end, chrom, strand))
#     stats = []
#     for accession, segments in accession_segments.items():
#         if accession in database_dict:
#             seq = database_dict[accession]
#             length = len(seq)
#             min_start, max_end = merge_segments(segments)
#             cov_length = (max_end-min_start+1) / 3
#             coverage = cov_length / length
#             chrom = segments[0][2]
#             strand = segments[0][3]
#             peptide_count = len(segments)
#             stats.append({
#                 'accession': accession,
#                 'chrom': chrom,
#                 'strand': strand,
#                 'min_start': min_start,
#                 'max_end': max_end,
#                 'coverage': coverage,
#                 'peptide_count': peptide_count,
#                 'sequence_length_aa': length
#             })
#     return stats
# def codon_test(seq, min_start, max_end, max_scan_nt):
#     seq = str(Seq(str(seq)).upper())
#     max_length = int(max_scan_nt / 3)
#     for i in range(max_length):
#         if seq[min_start-1-3*i:min_start+2-3*i] in START_CODONS:
#             phy_start = min_start-3*i
#             phy_value = CODON_PRIOR[seq[min_start-1-3*i:min_start+2-3*i]]
#             break
#     else:
#         phy_start = None
#         phy_value = None
#     for i in range(max_length):
#         if seq[max_end+i*3 : max_end+3*(i+1)] in STOP_CODONS:
#             phy_end = max_end + i*3
#             break
#     else:
#         phy_end = None
#     return phy_start, phy_value, phy_end
# def codon_test_minus(seq, min_start, max_end, max_scan_nt):
#     seq = str(Seq(str(seq)).upper())
#     max_length = int(max_scan_nt / 3)
#     for i in range(max_length):
#         if seq[min_start-1-3*(i+1):min_start-1-3*i] in MINUS_STOP_CODONS:
#             phy_start = min_start-3*i
#             break
#     else:
#         phy_start = None
#     for i in range(max_length):
#         if seq[max_end+3*(i-1):max_end+3*i] in MINUS_START_CODONS:
#             phy_end = max_end+3*i
#             phy_value = MINUS_CODON_PRIOR[seq[max_end+3*(i-1):max_end+3*i]]
#             break
#     else:
#         phy_end = None
#         phy_value = None
#     return phy_start, phy_value, phy_end
# def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank):
#     genome_seq = Seq(str(genome_seq).upper())
#     if strand == '+':
#         ctx = genome_seq[phy_start-1-flank:phy_start+flank+2]
#     else:
#         ctx = genome_seq[phy_end-3-flank:phy_end+flank].reverse_complement()
#     codon_pos = flank
#     core = 0.0
#     if ctx[codon_pos-3] in ('A', 'G'): core += 2.0
#     if ctx[codon_pos+3] == 'G':       core += 2.0
#     extended = 0.0
#     for rel, ref in { -6:'G', -5:'C', -4:'C', -2:'C', +4:'G' }.items():
#         if ctx[codon_pos+rel] == ref: extended += 0.5
#     return {"core": core, "extended": extended, "total": core+extended, "context": ctx}
# def run_scan_and_output_for_item(item):
#     global _genome_dict, _max_scan_nt
#     chrom = item['chrom']
#     strand = item['strand']
#     min_start = int(item['min_start'])
#     max_end = int(item['max_end'])
#     gseq = _genome_dict[chrom]

#     if strand == '+':
#         phy_start, phy_value, phy_end = codon_test(gseq, min_start, max_end, _max_scan_nt)
#         prior_triplet = str(gseq[phy_start-1:phy_start+2]) if phy_start else None
#     else:
#         phy_start, phy_value, phy_end = codon_test_minus(gseq, min_start, max_end, _max_scan_nt)
#         prior_triplet = str(gseq[phy_end-3:phy_end].reverse_complement()) if phy_end else None
#     item['phy_start'] = phy_start
#     item['phy_end'] = phy_end
#     item['prior'] = prior_triplet
#     item['total_score'] = None
#     if phy_start and phy_end:
#         kozak_results = kozak_score_cal(gseq, phy_start, phy_end, strand, flank=6)
#         item['total_score'] = phy_value + kozak_results['total']
#         item['kozak_seq'] = kozak_results['context']
#     return item
# def run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=100):
#     with Pool(processes=nproc, initializer=init_worker, initargs=(genome_file, max_scan_nt)) as pool:
#         stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
#     return stats_update
# def main():
#     peptide_file = "/media/wanglab/caca/Eu_peptido/20251113 horticulture research/initator_codon/sp_express_info.xlsx"
#     database_file = "/media/wanglab/caca/Eu_peptido/20251113 horticulture research/initator_codon/Eu_peptide_database_customized_5.fa"
#     genome_file = "/media/wanglab/caca/Eu_peptido/20251113 horticulture research/initator_codon/Eu_genome.fasta"
#     output_file = "/media/wanglab/caca/Eu_peptido/20251113 horticulture research/initator_codon/output.csv"    
#     stats = filter_peptide_seq_cal_cov(peptide_file, database_file)
#     stats_update = run_scan_and_output(stats, genome_file, max_scan_nt=300)    
#     df = pd.DataFrame(stats_update)
#     df.to_csv(output_file, index=False)
# if __name__ == "__main__":
#     main()

# # Kozak序列位置分析
# import pandas as pd
# file = r"F:\Eu_peptido\new_prepare\figure4d.xlsx"
# df = pd.read_excel(file, sheet_name="Sheet1")
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
# result_df.to_excel(r"F:\Eu_peptido\new_prepare\figure4e_seq_position_statistics.xlsx")

import pandas as pd
from collections import Counter
txt_file = "/Volumes/caca/test_fractionation/01figure/figure5/DEGs_ST_DT/dt_relation.txt"
excel_file = "/Volumes/caca/test_fractionation/01figure/sp_express_info.xlsx"
sequence_column = "sequence"
output_csv = "/Volumes/caca/test_fractionation/01figure/figure5/DEGs_ST_DT/dt_relation_amino_acid.csv"
try:
    with open(txt_file, 'r') as file:
        ids = {line.strip() for line in file if line.strip()}
    df = pd.read_excel(excel_file)
    if "ID" not in df.columns:
        raise ValueError("Excel文件中缺少'ID'列")
    filtered_df = df[df["ID"].isin(ids)]
    if sequence_column not in filtered_df.columns:
        raise ValueError(f"列 '{sequence_column}' 不存在于Excel文件中")
    first_amino_acids = []
    for seq in filtered_df[sequence_column]:
        if pd.notna(seq) and isinstance(seq, str) and len(seq) > 0:
            first_amino_acids.append(seq[0].upper())
    amino_acid_counts = Counter(first_amino_acids)
    sorted_counts = sorted(amino_acid_counts.items(), key=lambda x: x[0])
    result_df = pd.DataFrame(sorted_counts, columns=["Amino_Acid", "Count"])
    result_df.to_csv(output_csv, index=False)
    print(f"✅ 统计结果已保存到 {output_csv}")
    print("\n统计结果预览:")
    print(result_df.head())  # 输出前几行作为预览
except FileNotFoundError:
    print(f"错误: 文件 {excel_file} 或 {txt_file} 未找到")
except Exception as e:
    print(f"发生错误: {str(e)}")
