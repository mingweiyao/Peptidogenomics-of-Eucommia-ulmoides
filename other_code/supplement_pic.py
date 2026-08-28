# # figureS1: 不同组织染色体分布模式
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
#         karyotype.append({"Chr": chr_name, "Start": 0, "End": len(rec.seq)})
#     df = pd.DataFrame(karyotype)
#     df.to_csv(output_path, sep="\t", index=False)
#     print(f"✅ 染色体结构文件已保存至: {output_path}")
#     return df
# def sliding_window_density(gene_df, chr_len, window_size=2000000, step_size=1000000):
#     windows = []
#     for start in range(0, chr_len, step_size):
#         end = min(start + window_size, chr_len)
#         region = gene_df[(gene_df['start'] <= end) & (gene_df['end'] >= start)]
#         count = len(region)
#         windows.append([start, end, count])
#     return pd.DataFrame(windows, columns=["Start", "End", "Value"])
# def prepare_gene_density(excel_file, genome_file, density_output, karyo_output, window_size=2000000):
#     gff_df = pd.read_excel(excel_file, sheet_name="Leaf_NCPs")
#     chrom_name_map = {}
#     for rec in SeqIO.parse(genome_file, "fasta"):
#         gid = rec.id
#         name = next((x.split("=")[1] for x in rec.description.split("\t") if x.startswith("OriSeqID=")), gid)
#         chrom_name_map[gid] = name
#     gff_df["chrom"] = gff_df["chrom"].map(chrom_name_map)
#     karyo_df = calculate_karyotype(genome_file, karyo_output)
#     chr_lengths = dict(zip(karyo_df["Chr"], karyo_df["End"]))
#     all_windows = []
#     for chrom in sorted(gff_df["chrom"].unique()):
#         chr_len = chr_lengths.get(chrom)
#         if not chr_len:
#             continue
#         sub_df = gff_df[gff_df["chrom"] == chrom]
#         density_df = sliding_window_density(sub_df, chr_len, window_size)
#         density_df["Chr"] = chrom
#         all_windows.append(density_df)
#     result_df = pd.concat(all_windows, ignore_index=True)
#     result_df = result_df[["Chr", "Start", "End", "Value"]]
#     result_df.to_csv(density_output, sep="\t", index=False)
#     print(f"✅ 基因密度文件已保存至: {density_output} (使用原始计数)")
# if __name__ == "__main__":
#     excel_file = "G:/peptidegenomics/02figure/supplenment_pic/Peptides_info.xlsx"
#     genome_file = "D:/Desktop/小肽组/小肽组学/流程/1.构建小肽数据库/杜仲基因组Eu17/GWHBISF00000000.genome.fasta"
#     output_dir = "G:/peptidegenomics/02figure/supplenment_pic"
#     density_output = os.path.join(output_dir, "material_gene_density.tsv")
#     karyo_output = os.path.join(output_dir, "karyotype.tsv")
#     prepare_gene_density(excel_file, genome_file, density_output, karyo_output)

# # figureS2: 物理长度与染色体小肽数量相关性分析
# import pandas as pd
# from Bio import SeqIO
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import pearsonr
# import os
# genome_file = "D:/Desktop/小肽组/小肽组学/流程/1.构建小肽数据库/杜仲基因组Eu17/GWHBISF00000000.genome.fasta"
# excel_file = "G:/peptidegenomics/01analysis_data/sp_loc/Eu_sp_finally.xlsx"
# output_dir = "G:/peptidegenomics/02figure/supplenment_pic"
# os.makedirs(output_dir, exist_ok=True)
# mapping = []
# chr_sizes = {}
# for rec in SeqIO.parse(genome_file, "fasta"):
#     gwh_id = rec.id
#     desc_fields = rec.description.split("\t")
#     chr_name = None
#     for field in desc_fields:
#         if field.startswith("OriSeqID="):
#             chr_name = field.split("=")[1]
#         if field.startswith("Len="):
#             chr_length = int(field.split("=")[1])
#     if chr_name:
#         mapping.append({"chrom": gwh_id, "Chr": chr_name})
#         chr_sizes[chr_name] = chr_length
# mapping_df = pd.DataFrame(mapping)
# mapping_df.to_csv(os.path.join(output_dir, "chr_mapping.tsv"), sep="\t", index=False)
# gff_df = pd.read_excel(excel_file, sheet_name="CPs")
# gff_df = gff_df[(gff_df["chrom"] >= "GWHBISF00000485") & (gff_df["chrom"] <= "GWHBISF00000501")]
# gff_df = gff_df.merge(mapping_df, left_on="chrom", right_on="chrom", how="left")
# materials = gff_df["material"].unique()
# length_df = pd.DataFrame([{"Chr": chr_name, "Physical_Size": size / 1e6} for chr_name, size in chr_sizes.items()])
# for material in materials:
#     material_df = gff_df[gff_df["material"] == material]
#     peptide_count_df = material_df.groupby("Chr").size().reset_index(name="CP")
#     final_df = length_df.merge(peptide_count_df, on="Chr", how="inner").sort_values("Chr")
#     final_df.to_csv(os.path.join(output_dir, f"chrom_length_peptide_count_{material}.tsv"), sep="\t", index=False)
#     plt.figure(figsize=(8, 6))
#     sns.set(style="whitegrid")
#     ax = sns.regplot(
#         data=final_df, 
#         x="Physical_Size", 
#         y="CP", 
#         scatter_kws={'s': 60, 'color': '#117733'}, 
#         line_kws={'color': 'black'}
#     )
#     r, p = pearsonr(final_df["Physical_Size"], final_df["CP"])
#     plt.text(
#         x=max(final_df["Physical_Size"]) * 0.6,
#         y=max(final_df["CP"]) * 0.95,
#         s=f"r = {r:.2f}\np = {p:.2g}",
#         fontsize=12
#     )
#     plt.xlabel("Chromosome Length (Mb)", fontsize=12)
#     plt.ylabel("CP (Peptide Count)", fontsize=12)
#     plt.title(f"Peptide Count vs. Chromosome Length ({material})", fontsize=14)
#     plt.tight_layout()
#     plt.savefig(os.path.join(output_dir, f"Figure2_Length_vs_CP_{material}.pdf"))
#     plt.close()
# print("✅ 分析完成！按material分组的文件已保存。")

# # FigureS3: 相邻肽段距离 CPs & NCPs
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from brokenaxes import brokenaxes
# input_file = "G:/peptidegenomics/01analysis_data/sp_loc/Eu_sp_finally.xlsx"
# output_dir = "G:/peptidegenomics/02figure/supplenment_pic"
# target_material = "Root"
# df = pd.read_excel(input_file, sheet_name="CPs_calc")
# df = df[df["material"] == target_material].copy()
# df_sorted = df.sort_values(by=["chrom", "start"]).reset_index(drop=True)
# merged_peptides = []
# current_chrom, current_start, current_end = None, None, None
# for _, row in df_sorted.iterrows():
#     if current_chrom != row["chrom"]:
#         if current_chrom is not None:
#             merged_peptides.append({"chrom": current_chrom, "start": current_start, "end": current_end})
#         current_chrom, current_start, current_end = row["chrom"], row["start"], row["end"]
#     else:
#         if row["start"] <= current_end:
#             current_end = max(current_end, row["end"])
#         else:
#             merged_peptides.append({"chrom": current_chrom, "start": current_start, "end": current_end})
#             current_start, current_end = row["start"], row["end"]
# if current_chrom is not None:
#     merged_peptides.append({"chrom": current_chrom, "start": current_start, "end": current_end})
# df_merged = pd.DataFrame(merged_peptides)
# df_merged["next_start"] = df_merged.groupby("chrom")["start"].shift(-1)
# df_merged["distance"] = (df_merged["next_start"] - df_merged["end"]).clip(lower=0.1)
# valid_distances = df_merged[df_merged["distance"].notna()]["distance"] / 1000 
# output_csv = f"{output_dir}/CPs_{target_material}.csv"
# df_merged.to_csv(output_csv, index=False)
# print(f"💾 {target_material}数据已保存到: {output_csv}")
# plt.figure(figsize=(10, 6))
# bax = brokenaxes(ylims=((0, 4000), (77000, 78000)), hspace=0.05)
# bax.hist(valid_distances, bins=50, color="#009E73", edgecolor="black")
# bax.set_xlabel("Distance (kb)", fontsize=12)
# bax.set_ylabel("Frequency", fontsize=12)
# bax.set_title(f"Peptide Distance Distribution ({target_material})", pad=20)
# output_plot = f"{output_dir}/CPs_Distance_Histogram_{target_material}.pdf"
# plt.savefig(output_plot, bbox_inches="tight")
# plt.close()
# print(f"✅ {target_material}直方图已保存到: {output_plot}")
# # plt.figure(figsize=(10, 6))
# # plt.hist(valid_distances, bins=50, color="#009E73", edgecolor="black", alpha=0.7)
# # plt.xlabel("Distance (kb)", fontsize=12)
# # plt.ylabel("Frequency", fontsize=12)
# # plt.title(f"Peptide Distance Distribution ({target_material})", fontsize=14)
# # plt.grid(True, linestyle="--", alpha=0.6)
# # ymax = plt.ylim()[1] * 1.1
# # plt.ylim(0, ymax)
# # output_csv = f"{output_dir}/CPs_{target_material}.csv"
# # output_plot = f"{output_dir}/Figure3_Distance_Histogram_{target_material}.pdf"
# # df_merged.to_csv(output_csv, index=False)
# # plt.savefig(output_plot, bbox_inches="tight", dpi=300)
# # plt.close()

# # FigureS3_Peptide_TSS_Distance.pdf
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from brokenaxes import brokenaxes
# excel_file_tss = "G:/peptidegenomics/02figure/figure.xlsx"
# excel_file = "G:/peptidegenomics/01analysis_data/sp_loc/Eu_sp_finally.xlsx"
# output_dir = "G:/peptidegenomics/02figure/supplenment_pic"
# target_material = "Root"
# tss_df = pd.read_excel(excel_file_tss, sheet_name="TSS")
# gene_df = pd.read_excel(excel_file, sheet_name="CPs_calc")
# gene_df = gene_df[gene_df["material"] == target_material].copy()
# distances = []
# for idx, gene in gene_df.iterrows():
#     chrom, gene_start, strand = gene['chrom'], gene['start'], gene['strand']
#     tss_on_chr = tss_df[(tss_df['chrom'] == chrom) & (tss_df['strand'] == strand)]
#     min_distance = np.min(np.abs(tss_on_chr['start'] - gene_start)) if not tss_on_chr.empty else np.nan
#     distances.append(min_distance)
# gene_df['min_TSS_distance'] = distances
# gene_df = gene_df.dropna(subset=['min_TSS_distance'])
# gene_df['min_TSS_distance_kb'] = gene_df['min_TSS_distance'] / 1000
# output_csv = f"{output_dir}/CPs_with_TSS_distance_{target_material}.csv"
# gene_df.to_csv(output_csv, index=False)
# print(f"💾 结果已保存到: {output_csv}")
# bin_width = 1
# max_distance = int(np.ceil(gene_df['min_TSS_distance_kb'].max()))
# bins = np.arange(0, max_distance + bin_width, bin_width)
# hist_counts, bin_edges = np.histogram(gene_df['min_TSS_distance_kb'], bins=bins)
# freq_table = pd.DataFrame({
#     'Distance (kb)': [f"{dist:.0f}-{dist+bin_width:.0f}" for dist in bin_edges[:-1]],
#     'Frequency': hist_counts
# })
# freq_table.to_csv(f"{output_dir}/TSS_distance_frequency_{target_material}.csv", index=False)
# plt.figure(figsize=(10, 6))
# bax = brokenaxes(ylims=((0, 30000), (45000, 50000)), hspace=0.05)
# bax.hist(gene_df["min_TSS_distance_kb"], bins=50, color="#009E73", edgecolor="black")
# bax.set_xlabel("Distance to Nearest TSS (kb)", fontsize=12)
# bax.set_ylabel("Frequency", fontsize=12)
# bax.set_title(f"Peptide-TSS Distance Distribution ({target_material})", fontsize=14)
# output_pdf = f"{output_dir}/CPs_TSS_Distance_{target_material}.pdf"
# plt.savefig(output_pdf, bbox_inches="tight", dpi=300)
# plt.close()
# print(f"✅ 图表已保存到: {output_pdf}")
# # # --- 修改后的绘图部分（移除了Y轴截断）---
# # plt.figure(figsize=(10, 6))
# # plt.hist(
# #     gene_df["min_TSS_distance_kb"], 
# #     bins=50,
# #     color="#009E73",
# #     edgecolor="black",
# #     alpha=0.7
# # )
# # plt.xlabel("Distance to Nearest TSS (kb)", fontsize=12)
# # plt.ylabel("Frequency", fontsize=12)
# # plt.title(f"Peptide-TSS Distance Distribution ({target_material})", fontsize=14)
# # plt.grid(True, linestyle='--', alpha=0.6)
# # ymax = plt.ylim()[1] * 1.05
# # plt.ylim(0, ymax)
# # output_pdf = f"{output_dir}/CPs_TSS_Distance_{target_material}.pdf"
# # plt.savefig(output_pdf, bbox_inches="tight", dpi=300)
# # plt.close()
# # print(f"✅ 结果已保存:\n- 数据文件: {output_csv}\n- 图表文件: {output_pdf}")

import pandas as pd
excel_file = r'G:\peptidegenomics\02figure\figure7\NCP_accessions.xlsx'
csv_file = r'G:\peptidegenomics\02figure\figure7\yield_all\novel_genes_function_predictions.csv'
excel_df = pd.read_excel(excel_file)
csv_df = pd.read_csv(csv_file)
excel_df['NCPs'] = excel_df['NCPs'].str.strip()
csv_df['novel_gene'] = csv_df['novel_gene'].str.strip()
merged_df = pd.merge(csv_df, excel_df, left_on='novel_gene', right_on='NCPs', how='left')
output_file = r'D:\Desktop\yield_merged_output.xlsx'
merged_df.to_excel(output_file, index=False)
print(f"Merge completed and saved to {output_file}")
