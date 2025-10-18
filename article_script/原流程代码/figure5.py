# # figure1:数量-肽统计
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.backends.backend_pdf import PdfPages
# def plot_peptide_distribution(expression_matrix_path, output_pdf_path, threshold=5):
#     sns.set_theme(style="whitegrid")
#     plt.rcParams.update({
#         'font.size': 12,
#         'axes.titlepad': 20,
#         'axes.labelpad': 10
#     })
#     if expression_matrix_path.endswith('.csv'):
#         df = pd.read_csv(expression_matrix_path, index_col=0)
#     elif expression_matrix_path.endswith('.txt'):
#         df = pd.read_csv(expression_matrix_path, sep='\t', index_col=0)
#     elif expression_matrix_path.endswith(('.xls', '.xlsx')):
#         df = pd.read_excel(expression_matrix_path, index_col=0)
#     else:
#         raise ValueError("不支持的文件格式")
#     binary_df = (df >= threshold).astype(int)
#     peptide_counts = binary_df.sum(axis=1)
#     count_dist = peptide_counts.value_counts().sort_index()
#     plot_data = pd.DataFrame({
#         'Samples': count_dist.index,
#         'Peptides': count_dist.values,
#         'Percentage': (count_dist / len(peptide_counts) * 100).round(1)
#     })
#     colors = ["#1f78b4", "#f0e68c", "#e31a1c"]
#     cmap = LinearSegmentedColormap.from_list("blue_beige_red", colors)
#     norm = plt.Normalize(plot_data['Peptides'].min(), plot_data['Peptides'].max())
#     with PdfPages(output_pdf_path) as pdf:
#         fig1, ax1 = plt.subplots(figsize=(10, 6))
#         bars = ax1.bar(plot_data['Samples'], plot_data['Peptides'], 
#                        width=1.0,
#                        edgecolor='white', linewidth=0.5,
#                        color=[cmap(norm(val)) for val in plot_data['Peptides']])
#         x = plot_data['Samples'].values
#         y = plot_data['Peptides'].values
#         coeffs = np.polyfit(x, y, 10)
#         poly = np.poly1d(coeffs)
#         trend_x = np.linspace(min(x), max(x), 100)
#         ax1.plot(trend_x, poly(trend_x), 
#                  color='black', 
#                  linestyle='--', 
#                  linewidth=2,
#                  label=f'Trend (y={coeffs[0]:.2f}x²{coeffs[1]:+.2f}x{coeffs[2]:+.2f})')
#         ax1.set_title('Peptides Expressed in N Samples', fontweight='bold')
#         ax1.set_xlabel('Number of Samples', fontweight='bold')
#         ax1.set_ylabel('Number of Peptides', fontweight='bold')
#         ax1.set_xticks(plot_data['Samples'])
#         ax1.set_xticklabels(plot_data['Samples'].astype(int))
#         ax1.legend()
#         sns.despine()
#         plt.tight_layout()
#         pdf.savefig(fig1, dpi=300)
#         plt.close()
#         fig2, ax2 = plt.subplots(figsize=(10, 6))
#         bars = ax2.bar(plot_data['Samples'], plot_data['Percentage'],
#                       width=1.0,
#                       edgecolor='white', linewidth=0.5,
#                       color=[cmap(norm(val)) for val in plot_data['Peptides']])
#         ax2.set_title('Percentage of Peptides (%)', fontweight='bold')
#         ax2.set_xlabel('Number of Samples', fontweight='bold')
#         ax2.set_ylabel('Percentage (%)', fontweight='bold')
#         ax2.set_xticks(plot_data['Samples'])
#         ax2.set_xticklabels(plot_data['Samples'].astype(int))
#         ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
#         sns.despine()
#         plt.tight_layout()
#         pdf.savefig(fig2, dpi=300)
#         plt.close()
#     print(f"PDF已保存至: {output_pdf_path}")
#     return plot_data
# if __name__ == "__main__":
#     expression_file = "G:/peptidegenomics/02figure/sp_expressed.csv"
#     output_pdf =  "G:/peptidegenomics/02figure/figure5/gene_expression_distribution_t.pdf"
#     results = plot_peptide_distribution(expression_file, output_pdf)
#     print("\n肽段分布统计结果:")
#     print(results.to_string(index=False))

# # figure2：不同组织中肽表达数量
# import pandas as pd
# gene_expression_df = pd.read_csv("G:/peptidegenomics/02figure/sp_expressed.csv", index_col=0)
# tissue_mapping_df = pd.read_excel('G:/peptidegenomics/01analysis_data/rnaseq/rna-seq.xlsx', sheet_name="Sheet3")
# tissue_to_samples = tissue_mapping_df.groupby('Tissues')['Sample'].apply(list).to_dict()
# peptide_ids_by_tissue = {}
# for tissue, samples in tissue_to_samples.items():
#     tissue_expr = gene_expression_df[samples]
#     expressed_peptides = tissue_expr[(tissue_expr > 5).any(axis=1)].index.tolist()
#     peptide_ids_by_tissue[tissue] = expressed_peptides
# max_length = max(len(ids) for ids in peptide_ids_by_tissue.values())
# peptide_ids_df = pd.DataFrame({
#     tissue: ids + [None]*(max_length - len(ids)) 
#     for tissue, ids in peptide_ids_by_tissue.items()
# })
# output_csv_path = 'G:/peptidegenomics/02figure/figure5a/peptide_expression_ids.csv'
# peptide_ids_df.to_csv(output_csv_path, index=False)
# print("✅ 结果已保存至:", output_csv_path)
# print("\n示例数据：")
# print(peptide_ids_df.head())

# # figure3：不同来源肽chrom-scaffold统计
# import pandas as pd
# input_file = 'G:/peptidegenomics/02figure/sp_matrix_info.xlsx'
# output_file = "G:/peptidegenomics/02figure/figure5a/chromosome_scaffold_stats.csv"
# df = pd.read_excel(input_file)
# def classify_chrom(chrom_value):
#     chrom_str = str(chrom_value)
#     last_three = chrom_str[-3:] if len(chrom_str) >= 3 else chrom_str
#     try:
#         chrom_num = int(''.join(filter(str.isdigit, last_three)))
#         if 485 <= chrom_num <= 501:
#             return "chromosome"
#         else:
#             return "scaffold"
#     except:
#         return "scaffold"
# df['chrom_type'] = df['chrom'].apply(classify_chrom)
# result = df.groupby(['type', 'chrom_type']).size().unstack(fill_value=0)
# result = result[['scaffold', 'chromosome']] if 'chromosome' in result.columns else result
# result = result.reset_index()
# result.columns.name = None
# result.to_csv(output_file, index=False)
# print(f"统计结果已保存到: {output_file}")
# print("\n统计结果预览:")
# print(result)

# figure4: 不同比例染色体肽段分布情况
# import pandas as pd
# import numpy as np
# from Bio import SeqIO
# import os
# def calculate_karyotype(genome_file, output_path):
#     karyotype = []
#     for rec in SeqIO.parse(genome_file, "fasta"):
#         desc_parts = rec.description.split("\t")
#         chr_id = rec.id
#         chr_name = next((x.split("=")[1] for x in desc_parts if x.startswith("OriSeqID=")), chr_id)
#         karyotype.append({"chr": chr_name, "start": 0, "end": len(rec.seq)})
#     df = pd.DataFrame(karyotype)
#     df.to_csv(output_path, sep="\t", index=False)
#     print(f"✅ 染色体结构文件已保存至: {output_path}")
#     return df
# def percentile_bin_density(peptide_df, chr_len):
#     bins = np.linspace(0, chr_len, 11)
#     bins = [int(x) for x in bins]
#     intervals = []
#     for i in range(10):
#         start = bins[i]
#         end = bins[i+1]
#         count = len(peptide_df[
#             (peptide_df['start'] >= start) & 
#             (peptide_df['start'] < end)
#         ])
#         intervals.append([start, end, count])
#     return pd.DataFrame(intervals, columns=["start", "end", "count"])
# def prepare_peptide_density(excel_file, genome_file, density_output, karyo_output):
#     gff_df = pd.read_excel(excel_file)
#     gff_df = gff_df[(gff_df["chrom"] >= "GWHBISF00000485") & (gff_df["chrom"] <= "GWHBISF00000501")]  
#     chrom_name_map = {}
#     for rec in SeqIO.parse(genome_file, "fasta"):
#         gid = rec.id
#         name = next((x.split("=")[1] for x in rec.description.split("\t") if x.startswith("OriSeqID=")), gid)
#         chrom_name_map[gid] = name
#     gff_df["chrom"] = gff_df["chrom"].map(chrom_name_map)
#     karyo_df = calculate_karyotype(genome_file, karyo_output)
#     chr_lengths = dict(zip(karyo_df["chr"], karyo_df["end"]))
#     all_intervals = []
#     for chrom in sorted(gff_df["chrom"].unique()):
#         chr_len = chr_lengths.get(chrom)
#         if not chr_len:
#             continue
#         sub_df = gff_df[gff_df["chrom"] == chrom]
#         density_df = percentile_bin_density(sub_df, chr_len)
#         density_df["chr"] = chrom
#         density_df["percent_range"] = [f"{i*10}%-{(i+1)*10}%" for i in range(10)]
#         all_intervals.append(density_df)
#     result_df = pd.concat(all_intervals, ignore_index=True)
#     result_df = result_df[["chr", "percent_range", "start", "end", "count"]]
#     result_df.to_csv(density_output, sep="\t", index=False)
#     print(f"✅ 肽段密度分布文件已保存至: {density_output}")
#     print("\n示例输出：")
#     print(result_df.head())
# if __name__ == "__main__":
#     excel_file = 'G:/peptidegenomics/02figure/sp_matrix_info.xlsx'
#     genome_file = "D:/Desktop/小肽组/小肽组学/流程/1.构建小肽数据库/杜仲基因组Eu17/GWHBISF00000000.genome.fasta"
#     output_dir = "G:/peptidegenomics/02figure/figure5a"
#     os.makedirs(output_dir, exist_ok=True)
#     density_output = os.path.join(output_dir, "peptide_density_percentile.tsv")
#     karyo_output = os.path.join(output_dir, "karyotype.tsv")
#     prepare_peptide_density(excel_file, genome_file, density_output, karyo_output)

# Figure5_Peptide_TSS_Distance.pdf
import pandas as pd
import numpy as np
from brokenaxes import brokenaxes
import matplotlib.pyplot as plt
excel_file_tss = "G:/peptidegenomics/02figure/figure.xlsx"
excel_file = 'G:/peptidegenomics/02figure/sp_matrix_info.xlsx'
output_pdf = "G:/peptidegenomics/02figure/figure5/Figure_TSS_Gene_Distance_Distribution.pdf"
tss_df = pd.read_excel(excel_file_tss, sheet_name="TSS")
gene_df = pd.read_excel(excel_file)
distances = []
for idx, gene in gene_df.iterrows():
    chrom = gene['chrom']
    gene_start = gene['start']
    strand = gene['strand']
    tss_on_chr = tss_df[tss_df['chrom'] == chrom]
    if strand == "+":
        tss_on_chr = tss_on_chr[tss_on_chr['strand'] == "+"]
    else:
        tss_on_chr = tss_on_chr[tss_on_chr['strand'] == "-"]
    if not tss_on_chr.empty:
        min_distance = np.min(np.abs(tss_on_chr['start'] - gene_start))
        distances.append(min_distance)
    else:
        distances.append(np.nan)
gene_df['min_TSS_distance'] = distances
gene_df = gene_df.dropna(subset=['min_TSS_distance'])
gene_df['min_TSS_distance_kb'] = gene_df['min_TSS_distance'] / 1000

bin_width = 1
max_distance = int(np.ceil(gene_df['min_TSS_distance_kb'].max()))
bins = np.arange(0, max_distance + bin_width, bin_width)
hist_counts, bin_edges = np.histogram(gene_df['min_TSS_distance_kb'], bins=bins)
print("Distance (kb)\tFrequency")
for dist, count in zip(bin_edges[:-1], hist_counts):
    print(f"{dist:.0f}-{dist+bin_width:.0f} kb\t\t{count}")
pd.DataFrame({
    'Distance Range (kb)': [f"{dist:.0f}-{dist+bin_width:.0f}" for dist in bin_edges[:-1]],
    'Frequency': hist_counts
}).to_csv("G:/peptidegenomics/02figure/figure5/TSS_distance_frequency.csv", index=False)

fig = plt.figure(figsize=(8, 5))
bax = brokenaxes(ylims=((0, 60000),(150000, 160000)), hspace=0.05, fig=fig)
bax.hist(gene_df["min_TSS_distance_kb"], bins=50, color="#009E73", edgecolor="black")
bax.set_xlabel("Distance to Nearest TSS (kb)")
bax.set_ylabel("Frequency")
bax.set_title("Peptide to TSS Distance Distribution (Y-axis truncated)")
fig.savefig(output_pdf)
plt.close(fig)
print(f"✅ 截断图已保存为：{output_pdf}")

# # figure6: 统计起始氨基酸数量
# import pandas as pd
# from collections import Counter
# excel_file = 'G:/peptidegenomics/02figure/sp_matrix_info.xlsx'
# sheet_name = "Sheet1"
# sequence_column = "sequence"
# output_csv = "G:/peptidegenomics/02figure/figure5a/amino_acid_start_counts.csv"
# try:
#     df = pd.read_excel(excel_file, sheet_name=sheet_name)
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