# import pandas as pd
# def extract_data_by_id(id_file, data_file, output_file):
#     ids = pd.read_excel(id_file, header=None)[0].tolist()
#     data = pd.read_csv(data_file, header=0)
#     filtered_data = data[data['GeneID'].isin(ids)]
#     filtered_data.to_excel(output_file, index=False)
#     return filtered_data
# id_file = "/Volumes/caca/test_fractionation/02suplement/figure/figure2/common_id.xlsx"
# data_file = "/Volumes/caca/test_fractionation/01figure/sp_expressed.csv"
# output_file = "/Volumes/caca/test_fractionation/02suplement/figure/figure2/common_express.xlsx"
# result = extract_data_by_id(id_file, data_file, output_file)

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
#                         f'{y:,}', ha='center', va='bottom', fontsize=5)
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
#     expression_file = "/Volumes/caca/test_fractionation/02suplement/figure/figure2/common_express.xlsx"
#     output_pdf = "/Volumes/caca/test_fractionation/02suplement/figure/figure2/common_distribution.pdf"
#     results = plot_peptide_distribution(expression_file, output_pdf, threshold=10)
#     if results is not None:
#         print(f"\n{'='*60}")
#         print("详细分布统计结果:")
#         print(f"{'='*60}")
#         print(results.to_string(index=False, formatters={
#             'Number_of_Peptides': '{:,}'.format,
#             'Percentage': '{:.1f}%'.format
#         }))

# import pandas as pd
# def extract_data_by_sample(gene_expression_file):
#     df = pd.read_excel(gene_expression_file, sheet_name="peel")
#     peptide_ids_by_sample = {}
#     for sample in df.columns[1:]:
#         sample_data = df[['GeneID', sample]]
#         valid_gene_ids = sample_data[sample_data[sample] >= 10]['GeneID'].tolist()
#         peptide_ids_by_sample[sample] = valid_gene_ids
#     max_length = max(len(ids) for ids in peptide_ids_by_sample.values())
#     peptide_ids_df = pd.DataFrame({
#         sample: ids + [None] * (max_length - len(ids))
#         for sample, ids in peptide_ids_by_sample.items()
#     })
#     return peptide_ids_df
# gene_expression_file = "/Volumes/caca/test_fractionation/02suplement/figure/figure2/common_express.xlsx"
# output_file = "/Volumes/caca/test_fractionation/02suplement/figure/figure2/peel.xlsx"
# result_df = extract_data_by_sample(gene_expression_file)
# result_df.to_excel(output_file, index=False)
# print(f"✅ 结果已保存至: {output_file}")
