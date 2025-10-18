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
#     expression_file = "/Volumes/caca/test_fractionation/01figure/sp_expressed.csv"
#     output_pdf = "/Volumes/caca/test_fractionation/01figure/figure4/peptide_expression_distribution.pdf"
#     results = plot_peptide_distribution(expression_file, output_pdf, threshold=10)
#     if results is not None:
#         print(f"\n{'='*60}")
#         print("详细分布统计结果:")
#         print(f"{'='*60}")
#         print(results.to_string(index=False, formatters={
#             'Number_of_Peptides': '{:,}'.format,
#             'Percentage': '{:.1f}%'.format
#         }))

# # figure2：不同组织中肽表达数量
# import pandas as pd
# gene_expression_df = pd.read_csv("/Volumes/caca/test_fractionation/01figure/sp_expressed.csv", index_col=0)
# tissue_mapping_df = pd.read_excel('/Volumes/caca/test_fractionation/00raw/rnaseq/Total_rna_seq.xlsx', sheet_name="Sheet4")
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
# output_csv_path = '/Volumes/caca/test_fractionation/01figure/figure4/peptide_expression_ids.csv'
# peptide_ids_df.to_csv(output_csv_path, index=False)
# print("✅ 结果已保存至:", output_csv_path)
# print("\n示例数据：")
# print(peptide_ids_df.head())

# # figure63 统计起始氨基酸数量
# import pandas as pd
# from collections import Counter
# excel_file = "/Volumes/caca/test_fractionation/01figure/sp_express_info.xlsx"
# sequence_column = "sequence"
# output_csv = "/Volumes/caca/test_fractionation/01figure/figure4/amino_acid_start_counts.csv"
# try:
#     df = pd.read_excel(excel_file)
#     if sequence_column not in df.columns:
#         raise ValueError(f"列 '{sequence_column}' 不存在于Excel文件中")
#     first_amino_acids = []
#     for seq in df[sequence_column]:
#         if pd.notna(seq) and isinstance(seq, str) and len(seq) > 0:
#             first_amino_acids.append(seq[0].upper())
#     amino_acid_counts = Counter(first_amino_acids)
#     sorted_counts = sorted(amino_acid_counts.items(), key=lambda x: x[0])
#     result_df = pd.DataFrame(sorted_counts, columns=["Amino_Acid", "Count"])
#     result_df.to_csv(output_csv, index=False)
#     print(f"✅ 统计结果已保存到 {output_csv}")
#     print("\n统计结果预览:")
#     print(result_df.head())
# except FileNotFoundError:
#     print(f"错误: 文件 {excel_file} 未找到")
# except Exception as e:
#     print(f"发生错误: {str(e)}")

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
