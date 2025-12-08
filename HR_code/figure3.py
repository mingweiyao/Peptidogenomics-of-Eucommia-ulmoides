# # figure1: 染色体分布模式
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
#     gff_df = pd.read_excel(excel_file)
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
#     excel_file = "/Volumes/caca/test_fractionation/01figure/sp_express_info.xlsx"
#     genome_file = "/Volumes/caca/Eu_genome/GWHBISF00000000.genome.fasta"
#     output_dir = "/Volumes/caca/test_fractionation/01figure/figure3"
#     density_output = os.path.join(output_dir, "NCP_material_gene_density.tsv")
#     karyo_output = os.path.join(output_dir, "NCP_karyotype.tsv")
#     prepare_gene_density(excel_file, genome_file, density_output, karyo_output)

# # Figure3_Length_vs_NCP.pdf
# import pandas as pd
# from Bio import SeqIO
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import pearsonr
# import os
# genome_file = r"F:\Eu_peptido\00file\00raw\Raw_database\Eu_genome.fasta"
# excel_file = r"F:\Eu_peptido\new_prepare\sp_expression_info.xlsx"
# output_dir = r"F:\Eu_peptido\new_prepare"
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
# gff_df = pd.read_excel(excel_file)
# gff_df = gff_df[gff_df["chrom"] >= "GWHBISF00000485"]
# gff_df = gff_df[gff_df["chrom"] <= "GWHBISF00000501"]
# gff_df = gff_df.merge(mapping_df, left_on="chrom", right_on="chrom", how="left")
# length_df = pd.DataFrame([
#     {"Chr": chr_name, "Physical_Size": size / 1e6}
#     for chr_name, size in chr_sizes.items()
# ])
# peptide_count_df = gff_df.groupby("Chr").size().reset_index(name="NCP")
# final_df = length_df.merge(peptide_count_df, on="Chr", how="inner")
# final_df = final_df.sort_values("Chr")
# final_df.to_csv(os.path.join(output_dir, "chrom_length_peptide_count.tsv"), sep="\t", index=False)
# plt.figure(figsize=(8, 6))
# sns.set(style="whitegrid")
# ax = sns.regplot(data=final_df, x="Physical_Size", y="NCP", scatter_kws={'s': 60, 'color': '#117733'}, line_kws={'color': 'black'})
# r, p = pearsonr(final_df["Physical_Size"], final_df["NCP"])
# plt.text(
#     x=max(final_df["Physical_Size"]) * 0.6,
#     y=max(final_df["NCP"]) * 0.95,
#     s=f"r = {r:.2f}\np = {p:.2g}",
#     fontsize=12
# )
# plt.xlabel("Chromosome Length (Mb)", fontsize=12)
# plt.ylabel("NCP (Peptide Count)", fontsize=12)
# plt.title("Figure 2: Peptide Count vs. Chromosome Length", fontsize=14)
# plt.tight_layout()
# plt.savefig(os.path.join(output_dir, "Figure2_Length_vs_NCP.pdf"))
# print("✅ 输出文件已保存：chrom_length_peptide_count.tsv + Figure2_Length_vs_NCP.pdf")

# # Figure3_Peptide_Distance_Histogram.pdf
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from brokenaxes import brokenaxes
# input_file = r"F:\Eu_peptido\new_prepare\sp_expression_info.xlsx"
# output_csv = r"F:\Eu_peptido\new_prepare\NCPs.csv"
# output_plot = r"F:\Eu_peptido\new_prepare\Figure3_Peptide_Distance_Histogram.pdf"
# df = pd.read_excel(input_file)
# df_sorted = df.sort_values(by=["chrom", "start"]).reset_index(drop=True)
# merged_peptides = []
# current_chrom = None
# current_start = None
# current_end = None
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
# df_merged["distance"] = df_merged["next_start"] - df_merged["end"]
# df_merged["distance"] = np.where(df_merged["distance"] < 0, 0.1, df_merged["distance"])
# valid_distances = df_merged[df_merged["distance"].notna()]["distance"] / 1000
# df_merged.to_csv(output_csv, index=False)
# print(f"💾 合并后的肽段及距离已保存到：{output_csv}")
# fig = plt.figure(figsize=(8, 5))
# bax = brokenaxes(ylims=((0, 4000), (17000, 19000)), hspace=0.05, fig=fig)
# bax.hist(valid_distances, bins=50, color="#009E73", edgecolor="black")
# bax.set_xlabel("Distance (kb)", fontsize=12)
# bax.set_ylabel("Frequency", fontsize=12)
# bax.set_title("Distance Between Adjacent Peptides (Negative Distances Set to 0.1 kb)", pad=20)
# fig.savefig(output_plot, bbox_inches='tight')
# plt.close(fig)
# print(f"✅ 截断图已保存为：{output_plot}")

# # Figure4_Peptide_TSS_Distance.pdf
# import pandas as pd
# import numpy as np
# from brokenaxes import brokenaxes
# import matplotlib.pyplot as plt
# excel_file_tss = r"F:\Eu_peptido\new_prepare\figure3_figure.xlsx"
# excel_file = r"F:\Eu_peptido\new_prepare\sp_expression_info.xlsx"
# output_pdf = r"F:\Eu_peptido\new_prepare\Figure_TSS_Gene_Distance_Distribution.pdf"
# tss_df = pd.read_excel(excel_file_tss, sheet_name="TSS")
# gene_df = pd.read_excel(excel_file)
# distances = []
# for idx, gene in gene_df.iterrows():
#     chrom = gene['chrom']
#     gene_start = gene['start']
#     strand = gene['strand']
#     tss_on_chr = tss_df[tss_df['chrom'] == chrom]
#     if strand == "+":
#         tss_on_chr = tss_on_chr[tss_on_chr['strand'] == "+"]
#     else:
#         tss_on_chr = tss_on_chr[tss_on_chr['strand'] == "-"]
#     if not tss_on_chr.empty:
#         min_distance = np.min(np.abs(tss_on_chr['start'] - gene_start))
#         distances.append(min_distance)
#     else:
#         distances.append(np.nan)
# gene_df['min_TSS_distance'] = distances
# gene_df = gene_df.dropna(subset=['min_TSS_distance'])
# gene_df['min_TSS_distance_kb'] = gene_df['min_TSS_distance'] / 1000
# bin_width = 1
# max_distance = int(np.ceil(gene_df['min_TSS_distance_kb'].max()))
# bins = np.arange(0, max_distance + bin_width, bin_width)
# hist_counts, bin_edges = np.histogram(gene_df['min_TSS_distance_kb'], bins=bins)
# print("Distance (kb)\tFrequency")
# for dist, count in zip(bin_edges[:-1], hist_counts):
#     print(f"{dist:.0f}-{dist+bin_width:.0f} kb\t\t{count}")
# pd.DataFrame({
#     'Distance Range (kb)': [f"{dist:.0f}-{dist+bin_width:.0f}" for dist in bin_edges[:-1]],
#     'Frequency': hist_counts
# }).to_csv(r"F:\Eu_peptido\new_prepare\SS_distance_frequency.csv", index=False)
# fig = plt.figure(figsize=(8, 5))
# bax = brokenaxes(ylims=((0, 6000),(12000, 14000)), hspace=0.05, fig=fig)
# bax.hist(gene_df["min_TSS_distance_kb"], bins=50, color="#009E73", edgecolor="black")
# bax.set_xlabel("Distance to Nearest TSS (kb)")
# bax.set_ylabel("Frequency")
# bax.set_title("Peptide to TSS Distance Distribution (Y-axis truncated)")
# fig.savefig(output_pdf)
# plt.close(fig)
# print(f"✅ 截断图已保存为：{output_pdf}")

# # figure5: chromosome and scaffold
# import pandas as pd
# excel_file = "/Volumes/caca/test_fractionation/01figure/sp_express_info.xlsx"
# df = pd.read_excel(excel_file)
# chrom_range = df[(df['chrom'] >= 'GWHBISF00000485') & (df['chrom'] <= 'GWHBISF00000501')]
# result = df.groupby('type').apply(
#     lambda group: pd.Series({
#         'In_Range': ((group['chrom'] >= 'GWHBISF00000485') & (group['chrom'] <= 'GWHBISF00000501')).sum(),
#         'Out_of_Range': (~((group['chrom'] >= 'GWHBISF00000485') & (group['chrom'] <= 'GWHBISF00000501'))).sum()
#     })
# ).reset_index()
# print(result)

# # figure6: strand
# import pandas as pd
# excel_file = "/Volumes/caca/test_fractionation/01figure/sp_express_info.xlsx"
# df = pd.read_excel(excel_file)
# result = df.groupby('type').apply(
#     lambda group: pd.Series({
#         '+': (group['strand'] == '+').sum(),
#         '-': (~(group['strand'] == '+')).sum()
#     })
# ).reset_index()
# print(result)

# # figure7: strand
# import pandas as pd
# excel_file = "/Volumes/caca/test_fractionation/01figure/sp_express_info.xlsx"
# df = pd.read_excel(excel_file)
# result = df.groupby('type').apply(
#     lambda group: pd.Series({
#         '1': (group['frame'] == 1).sum(),
#         '2': (group['frame'] == 2).sum(),
#         '3': (group['frame'] == 3).sum(),
#         '4': (group['frame'] == 4).sum(),
#         '5': (group['frame'] == 5).sum(),
#         '6': (group['frame'] == 6).sum()
#     })
# ).reset_index()
# print(result)