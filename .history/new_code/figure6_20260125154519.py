# # 按照FC和Pvalue筛选表达基因
# import pandas as pd
# def filter_excel_by_log2fc_pvalue(input_excel, output_excel, log2fc_col="log2FC", 
#                                   pvalue_col="PValue", log2fc_cutoff=1, pvalue_cutoff=0.05):
#     xls = pd.ExcelFile(input_excel)
#     result_dfs = []
#     for sheet in xls.sheet_names:
#         filtered = []
#         df = pd.read_excel(xls, sheet_name=sheet)
#         if log2fc_col not in df.columns or pvalue_col not in df.columns:
#             print(f"[跳过] sheet '{sheet}' 缺少 {log2fc_col} 或 {pvalue_col}")
#             continue
#         # 转为数值，避免字符串/空值问题
#         df[log2fc_col] = pd.to_numeric(df[log2fc_col], errors="coerce")
#         df[pvalue_col] = pd.to_numeric(df[pvalue_col], errors="coerce")
#         # 筛选条件
#         df_f = df[(df[log2fc_col].abs() > log2fc_cutoff) & (df[pvalue_col] < pvalue_cutoff)].copy()
#         if df_f.empty: continue
#         filtered = df_f[["name", log2fc_col, pvalue_col]].copy()
#         filtered["source_sheet"] = sheet
#         result_dfs.append(filtered)
#     if not result_dfs:
#         print("没有任何 sheet 满足筛选条件")
#         return
#     final_df = pd.concat(result_dfs, ignore_index=True)
#     # 写入新的 Excel
#     final_df.to_excel(output_excel, index=False)
#     print(f"筛选完成，共输出 {len(final_df)} 行")
#     print(f"结果已写入: {output_excel}")
# input_excel = r"F:\work_mechanism\new_file\02figure\figure6\pep_deseq.xlsx"
# output_excel = r"F:\work_mechanism\new_file\02figure\figure6\pep_deseq_diff.xlsx"
# filter_excel_by_log2fc_pvalue(input_excel, output_excel)

# # FPKM标准化
# import os
# import pandas as pd
# import gffutils
# GFF_FILE = r"F:\work_mechanism\new_file\02figure\figure5\codon\codon_prediction\codon_prediction_v7\sp_codon.gff"
# COUNT_FILE = r"F:\work_mechanism\new_file\02figure\figure6\WGCNA\dt_st_count.xlsx"
# OUT_DIR = r"F:\work_mechanism\new_file\02figure\figure6\WGCNA"
# def prepare_length_data(gff_file):
#     db_path = gff_file + ".db"
#     if not os.path.exists(db_path):
#         print("🔄 正在创建 GFF 数据库...")
#         gffutils.create_db(
#             gff_file,
#             dbfn=db_path,
#             force=True,
#             keep_order=True,
#             merge_strategy="merge",
#             id_spec={"gene": "ID", "mRNA": "ID", "CDS": "Parent"},
#             disable_infer_genes=True,
#             disable_infer_transcripts=True
#         )
#     db = gffutils.FeatureDB(db_path)
#     gene_lengths = {}
#     for gene in db.features_of_type("gene"):
#         total_length = 0
#         for mrna in db.children(gene, featuretype="mRNA", order_by="start"):
#             exons = list(db.children(mrna, featuretype="exon", order_by="start"))
#             if exons:
#                 mrna_len = sum(e.end - e.start + 1 for e in exons)
#                 total_length += mrna_len
#         if total_length > 0:
#             gene_id = gene.id.replace("evm.model.", "evm.TU.")
#             gene_lengths[gene_id] = total_length
#     length_df = pd.DataFrame(gene_lengths.items(), columns=["GeneID", "length"])
#     out_len = os.path.join(OUT_DIR, "gene_lengths.csv")
#     length_df.to_csv(out_len, index=False)
#     print(f"✅ 基因长度表已生成：{out_len}")
#     return length_df
# def read_counts(count_file):
#     if count_file.lower().endswith(".csv"):
#         df = pd.read_csv(count_file)
#     else:
#         df = pd.read_excel(count_file)
#     if df.columns[0] != "GeneID":
#         df = df.rename(columns={df.columns[0]: "GeneID"})
#     # 确保 GeneID 为字符串且去空格
#     df["GeneID"] = df["GeneID"].astype(str).str.strip()
#     return df
# def normalize_fpkm(count_df, length_df, output_file):
#     length_df = length_df.copy()
#     length_df["GeneID"] = length_df["GeneID"].astype(str).str.strip()
#     df = pd.merge(count_df, length_df, on="GeneID", how="inner")
#     df = df[df["length"] > 0].copy()
#     sample_cols = [c for c in df.columns if c not in ["GeneID", "length"]]
#     fpkm_data = {}
#     for sample in sample_cols:
#         counts = pd.to_numeric(df[sample], errors="coerce").fillna(0)
#         N = counts.sum()
#         if N == 0:
#             # 全是0时，直接输出0
#             fpkm = counts * 0.0
#         else:
#             L = df["length"].astype(float)  # bp
#             fpkm = (1e9 * counts) / (N * L)
#         fpkm_data[sample] = fpkm
#     fpkm_df = pd.concat([df[["GeneID"]], pd.DataFrame(fpkm_data)], axis=1)
#     fpkm_df.to_excel(output_file, index=False)
#     print(f"✅ FPKM 输出完成：{output_file}")
# if __name__ == "__main__":
#     os.makedirs(OUT_DIR, exist_ok=True)
#     length_df = prepare_length_data(GFF_FILE)
#     count_df = read_counts(COUNT_FILE)
#     normalize_fpkm(count_df, length_df, os.path.join(OUT_DIR, "dt_st_count_fpkm.xlsx"))

# # id提取
# import pandas as pd
# input_file = r"F:\work_mechanism\new_file\02figure\figure6\veen_analysis.xlsx"
# deseq_file = r"F:\work_mechanism\new_file\02figure\figure6\pep_deseq_diff.xlsx"
# output_file = r"F:\work_mechanism\new_file\02figure\figure6\st_target.xlsx"
# target_id = pd.read_excel(input_file, sheet_name="st")['ID']
# deseq = pd.read_excel(deseq_file, sheet_name="ST")
# sub_deseq = deseq[deseq['name'].isin(target_id)]
# sub_deseq.to_excel(output_file, index=False)

# # 提取相关肽的信息，用于肽分类的构建
# import pandas as pd
# id_file = r"F:\work_mechanism\new_file\02figure\figure6\dt_target.xlsx"
# temp_file = r"F:\work_mechanism\new_file\02figure\figure5\codon\codon_prediction\codon_prediction_v7\candidates_scored.xlsx"
# expression_info = r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_info.xlsx"
# df_id = pd.read_excel(id_file, sheet_name="st")
# temp_file_df = pd.read_excel(temp_file)
# df_express = pd.read_excel(expression_info, sheet_name="unique")
# id_list = df_id['name']
# id_temp = temp_file_df[temp_file_df['ID'].isin(id_list)]
# id_express = df_express[df_express['ID'].isin(id_temp['pep_id'])]
# id_express.to_excel(r"F:\work_mechanism\new_file\02figure\figure6\st_pep_info.xlsx", index=False)

# # 关联基因数提取
# import os
# import numpy as np
# from collections import defaultdict
# from pathlib import Path
# # 1. 提取高TOM值的NCP ID并汇总
# def extract_novel_data(folder_path):
#     novel_data = set()
#     for filename in os.listdir(folder_path):
#         if filename.endswith('.txt'):
#             file_path = os.path.join(folder_path, filename)
#             try:
#                 with open(file_path, 'r', encoding='utf-8') as file:
#                     for line in file:
#                         columns = line.strip().split()
#                         novel_data.add(columns[0])
#             except Exception as e:
#                 print(f"处理文件 {filename} 时出错: {e}")
#     return sorted(novel_data)
# # 2. 统计TOM值分布并输出统计区间
# def count_tom_ranges(folder_path, output_file):
#     all_tom_values = []
#     for filename in os.listdir(folder_path):
#         if filename.endswith('.txt'):
#             file_path = os.path.join(folder_path, filename)
#             with open(file_path, 'r') as file:
#                 next(file)
#                 for line in file:
#                     parts = line.strip().split('\t')
#                     if len(parts) >= 3:
#                         try:
#                             tom = float(parts[2])
#                             all_tom_values.append(tom)
#                         except ValueError:
#                             continue
#     if not all_tom_values:
#         print("警告：未找到有效的TOM值")
#         return
#     min_tom = min(all_tom_values)
#     max_tom = max(all_tom_values)
#     bins = np.linspace(min_tom, max_tom, 11)
#     range_counts = defaultdict(int)
#     # 对每个文件再次统计TOM值所在区间的分布
#     for filename in os.listdir(folder_path):
#         if filename.endswith('.txt'):
#             file_path = os.path.join(folder_path, filename)
#             with open(file_path, 'r') as file:
#                 next(file)
#                 for line in file:
#                     parts = line.strip().split('\t')
#                     if len(parts) >= 3:
#                         try:
#                             tom = float(parts[2])
#                             for i in range(len(bins) - 1):
#                                 if bins[i] <= tom < bins[i + 1]:
#                                     range_counts[(bins[i], bins[i + 1])] += 1
#                                     break
#                             if tom == bins[-1]:
#                                 range_counts[(bins[-2], bins[-1])] += 1
#                         except ValueError:
#                             continue
#     with open(output_file, 'w') as out_file:
#         out_file.write("TOM区间\t数据行数\n")
#         sorted_ranges = sorted(range_counts.items(), key=lambda x: x[0][0])
#         for (start, end), count in sorted_ranges:
#             out_file.write(f"[{start:.4f}, {end:.4f})\t{count}\n")
#     print(f"TOM值区间统计结果已保存到 {output_file}")
#     print(f"TOM值范围: {min_tom:.4f} 到 {max_tom:.4f}")
# # 3. 统计肽连接基因数分布频率
# def process_and_extract_genes(input_folder, output_data_file, output_freq_file, min_count, max_count):
#     gene_counts = defaultdict(int)
#     gene_lines = defaultdict(list)
#     for filename in os.listdir(input_folder):
#         if filename.endswith('.txt'):
#             file_path = os.path.join(input_folder, filename)
#             try:
#                 with open(file_path, 'r', encoding='utf-8') as f:
#                     next(f)  # 跳过第一行
#                     for line in f:
#                         line = line.strip()
#                         if line:
#                             gene = line.split('\t')[0]
#                             gene_counts[gene] += 1
#                             gene_lines[gene].append(line)
#             except Exception as e:
#                 print(f"处理文件 {filename} 时出错: {e}")
#     # 筛选出现次数在指定范围内的基因
#     filtered_genes = {
#         gene: count for gene, count in gene_counts.items()
#         if min_count <= count <= max_count
#     }
#     # 输出筛选后的数据
#     with open(output_data_file, 'w', encoding='utf-8') as f_data:
#         f_data.write("\n".join(
#             line for gene in filtered_genes
#             for line in gene_lines[gene]
#         ))
#     # 频率分布统计
#     freq_dist = defaultdict(int)
#     for count in gene_counts.values():
#         if min_count <= count <= max_count:
#             freq_dist[count] += 1
#     # 输出频率统计结果
#     with open(output_freq_file, 'w', encoding='utf-8') as f_freq:
#         f_freq.write("出现次数\t基因数量\n")
#         for count, num_genes in sorted(freq_dist.items()):
#             f_freq.write(f"{count}\t{num_genes}\n")
# # 4. 提取并合并共有NCP基因对信息
# def merge_files_by_id(folder_path, id_file_path, output_file_path):
#     with open(id_file_path, 'r', encoding='utf-8') as id_file:
#         target_ids = set(line.strip() for line in id_file if line.strip())
#     with open(output_file_path, 'w', encoding='utf-8') as output_file:
#         for filename in os.listdir(folder_path):
#             if filename.endswith('.txt') and filename != os.path.basename(id_file_path):
#                 file_path = os.path.join(folder_path, filename)
#                 try:
#                     with open(file_path, 'r', encoding='utf-8') as input_file:
#                         for line in input_file:
#                             line = line.strip()
#                             if line:
#                                 first_col = line.split()[0]
#                                 if first_col in target_ids:
#                                     output_file.write(line + '\n')
#                 except Exception as e:
#                     print(f"处理文件 {filename} 时出错: {e}")
# # 5. 执行所有功能
# def main():
#     folder_path = "/Volumes/caca/work_mechanism/new_file/02figure/figure6/WGCNA/module_specific_results"
#     output_base = "/Volumes/caca/work_mechanism/new_file/02figure/figure6/WGCNA/connect"
#     # 汇总高TOM值的NCPs
#     novel_data_output = os.path.join(output_base, "wt_tom_all_id.txt")
#     unique_novel_data = extract_novel_data(folder_path)
#     with open(novel_data_output, 'w', encoding='utf-8') as out_file:
#         for item in unique_novel_data:
#             out_file.write(item + '\n')
#     # TOM值分布统计
#     tom_range_output = os.path.join(output_base, "wt_tom_range_statistics.txt")
#     count_tom_ranges(folder_path, tom_range_output)
#     # 统计基因连接频率
#     gene_data_output = os.path.join(output_base, "wt_tom_all_info.txt")
#     gene_freq_output = os.path.join(output_base, "wt_tom_frequency.txt")
#     process_and_extract_genes(folder_path, gene_data_output, gene_freq_output, min_count=1, max_count=30000)
#     # 提取并合并共有的NCP基因对信息
#     id_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure6/WGCNA/common_id.txt"
#     merge_output = os.path.join(output_base, "common_id_tom.txt")
#     merge_files_by_id(folder_path, id_file, merge_output)
#     print("所有任务已完成！")
# if __name__ == "__main__":
#     main()

import pandas as pd
id_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure6/veen_analysis.xlsx"
expression_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure6/WGCNA/dt_st_count_fpkm.xlsx"
df_id = pd.read_excel(id_file, sheet_name="st")
expression_file = pd.read_excel(expression_file)
target_df = expression_file[expression_file["GeneID"].isin(df_id["ID"])]
target_df.to_excel("/Volumes/caca/work_mechanism/new_file/02figure/figure6/st_fpkm.xlsx", index=False)