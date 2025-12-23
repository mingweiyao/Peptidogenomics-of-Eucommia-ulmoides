# # 将Excel中的基因注释信息转换为GFF3格式文件
# from datetime import datetime
# import pandas as pd
# def excel_to_gff3(df, output_gff):
#     required_columns = ['ID', 'start', 'end', 'strand', 'chrom']
#     missing_cols = [col for col in required_columns if col not in df.columns]
#     if missing_cols:
#         print(f"错误: 缺少必要的列: {', '.join(missing_cols)}")
#         return
#     gff_header = f"""##gff-version 3
# ##date {datetime.now().strftime('%Y-%m-%d')}
# ##source {df}
# ##genome-build v1.0
# """
#     gff_lines = []
#     for idx, row in df.iterrows():
#         seqid = row['chrom']
#         source = "EuNCP"
#         start = int(row['start'])
#         end = int(row['end'])
#         strand = row['strand']
#         gene_id = row['ID']
#         gene_line = (
#             f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={gene_id};Name={gene_id}"
#         )
#         gff_lines.append(gene_line)
#         mrna_id = f"{gene_id}.t1"
#         mrna_line = (
#             f"{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={mrna_id};Parent={gene_id};product=predicted protein"
#         )
#         gff_lines.append(mrna_line)
#         exon_id = f"{mrna_id}.exon1"
#         exon_line = (
#             f"{seqid}\t{source}\texon\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={exon_id};Parent={mrna_id}"
#         )
#         gff_lines.append(exon_line)
#         cds_id = f"{mrna_id}.cds"
#         cds_line = (
#             f"{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
#             f"ID={cds_id};Parent={mrna_id}"
#         )
#         gff_lines.append(cds_line)
#     try:
#         with open(output_gff, 'w') as f:
#             f.write(gff_header)
#             f.write("\n".join(gff_lines))
#         print(f"成功生成GFF3文件: {output_gff}")
#         print(f"转换了 {len(df)} 条基因记录，共生成 {len(gff_lines)} 行GFF记录")
#     except Exception as e:
#         print(f"写入GFF文件失败: {e}")
# excel_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_sp_info.xlsx"
# output_gff = r"D:\Desktop\peptidemicro\00file\01figure\figure3\temp_file\unique.gff"
# df = pd.read_excel(excel_file, sheet_name='unique')
# excel_to_gff3(df, output_gff)

# # Length_vs_NCP.pdf
# import pandas as pd
# from Bio import SeqIO
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import pearsonr
# import os
# genome_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\Eu_genome.fasta"
# excel_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_sp_info.xlsx"
# output_dir = r"D:\Desktop\peptidemicro\00file\01figure\figure3"
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
# gff_df = pd.read_excel(excel_file, sheet_name='unique_chrom')
# length_df = pd.DataFrame([
#     {"chrom": chr_name, "Physical_Size": size / 1e6}
#     for chr_name, size in chr_sizes.items()
# ])
# peptide_count_df = gff_df.groupby("chrom").size().reset_index(name="NCP")
# final_df = length_df.merge(peptide_count_df, on="chrom", how="inner")
# final_df = final_df.sort_values("chrom")
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
# plt.savefig(os.path.join(output_dir, "Length_vs_NCP.pdf"))
# print("✅ 输出文件已保存：chrom_length_peptide_count.tsv + Length_vs_NCP.pdf")

# # Peptide_Distance_Histogram.pdf
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from brokenaxes import brokenaxes
# input_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_sp_info.xlsx"
# output_csv = r"D:\Desktop\peptidemicro\00file\01figure\figure3\NCPs.csv"
# output_plot = r"D:\Desktop\peptidemicro\00file\01figure\figure3\Peptide_Distance_Histogram.pdf"
# df = pd.read_excel(input_file, sheet_name="unique")
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
# bax = brokenaxes(ylims=((0, 3000), (3001,3002)), hspace=0.05, fig=fig)
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
# excel_file_tss = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_figure.xlsx"
# excel_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_sp_info.xlsx"
# output_pdf = r"D:\Desktop\peptidemicro\00file\01figure\figure3\TSS_Gene_Distance_Distribution.pdf"
# tss_df = pd.read_excel(excel_file_tss, sheet_name="TSS")
# gene_df = pd.read_excel(excel_file, sheet_name="unique")
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
# }).to_csv(r"D:\Desktop\peptidemicro\00file\01figure\figure3\TSS_distance_frequency.csv", index=False)
# fig = plt.figure(figsize=(8, 5))
# bax = brokenaxes(ylims=((0, 3000),(3001, 3002)), hspace=0.05, fig=fig)
# bax.hist(gene_df["min_TSS_distance_kb"], bins=50, color="#009E73", edgecolor="black")
# bax.set_xlabel("Distance to Nearest TSS (kb)")
# bax.set_ylabel("Frequency")
# bax.set_title("Peptide to TSS Distance Distribution (Y-axis truncated)")
# fig.savefig(output_pdf)
# plt.close(fig)
# print(f"✅ 截断图已保存为：{output_pdf}")