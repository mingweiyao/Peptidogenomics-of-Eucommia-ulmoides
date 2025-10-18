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
# # 3. 统计肽连接基因数分布频率
# def process_and_extract_genes(input_folder, output_data_file, output_freq_file, min_count, max_count):
#     gene_counts = defaultdict(int)
#     gene_lines = defaultdict(list)
#     for filename in os.listdir(input_folder):
#         if filename.endswith('.txt'):
#             file_path = os.path.join(input_folder, filename)
#             try:
#                 with open(file_path, 'r', encoding='utf-8') as f:
#                     next(f)
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
# # 4. 统计肽连接天然橡胶和绿原酸合成基因数分布频率
# def process_and_extract_genes_rubber_cga(input_folder, output_data_file, min_count, max_count):
#     gene_counts = defaultdict(int)
#     gene_lines = defaultdict(list)
#     for filename in os.listdir(input_folder):
#         if filename.endswith('.txt'):
#             file_path = os.path.join(input_folder, filename)
#             try:
#                 with open(file_path, 'r', encoding='utf-8') as f:
#                     next(f)
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
# # 7. 执行所有功能
# def main():
#     folder_path = "/Volumes/caca/test_fractionation/01figure/WGCNA/WGCNA_Results_yield/module_specific_results"
#     output_base = "/Users/lemon/Desktop"
#     # 汇总高TOM值的NCPs
#     novel_data_output = os.path.join(output_base, "yield_tom_all_id.txt")
#     unique_novel_data = extract_novel_data(folder_path)
#     with open(novel_data_output, 'w', encoding='utf-8') as out_file:
#         for item in unique_novel_data:
#             out_file.write(item + '\n')
#     # TOM值分布统计
#     tom_range_output = os.path.join(output_base, "yield_tom_range_statistics.txt")
#     count_tom_ranges(folder_path, tom_range_output)
#     # 统计基因连接频率
#     gene_data_output = os.path.join(output_base, "yield_tom_all_info.txt")
#     gene_freq_output = os.path.join(output_base, "yield_tom_frequency.txt")
#     process_and_extract_genes(folder_path, gene_data_output, gene_freq_output, min_count=1, max_count=30000)
#     # 统计天然橡胶合成和绿原酸调控基因
#     rubber_cga_path = "/Volumes/caca/test_fractionation/01figure/WGCNA/WGCNA_Results_yield/module_rubber_gca"
#     gene_data_output_rubber_cga = os.path.join(output_base, "rubber_cga_tom_all_info.txt")
#     process_and_extract_genes_rubber_cga(rubber_cga_path, gene_data_output_rubber_cga, min_count=1, max_count=30000)
#     print("所有任务已完成！")
# if __name__ == "__main__":
#     main()

# tom_file = "/Volumes/caca/test_fractionation/01figure/figure6/rubber_cga_tom_all_info.txt"
# output_file = "/Volumes/caca/test_fractionation/01figure/figure6/rubber_cga_ids.txt"
# sp_id = set()
# with open(tom_file, "r", encoding='utf-8') as tom:
#     for line in tom:
#         sp_id.add(line.split()[0])
# with open(output_file, "w", encoding='utf-8') as f:
#     for item in sp_id:
#         f.write(f"{item}\n")
# print(f"唯一ID数量: {len(sp_id)}")
# print(f"已保存到: {output_file}")

# import pandas as pd
# def rubber_cga_info_extract(input_file, output_file):
#     rubber_cga_info = pd.read_excel(input_file, sheet_name="Sheet2")
#     id_list = pd.read_excel(input_file, sheet_name="Sheet1")["ID"]
#     filtered = rubber_cga_info[rubber_cga_info.iloc[:, 0].isin(id_list)]
#     filtered.to_excel(output_file, index=False)
# input_file = "/Volumes/caca/test_fractionation/01figure/figure6/yield_CGA.xlsx"
# output_file = "/Volumes/caca/test_fractionation/01figure/figure6/rubber_cga_info.xlsx"
# rubber_cga_info_extract(input_file, output_file)

# import pandas as pd
# id_file = "/Volumes/caca/test_fractionation/01figure/figure6/yield_CGA.xlsx"
# GO_file = "/Volumes/caca/test_fractionation/01figure/GO/yield/novel_genes_function_predictions.csv"
# output_file = "/Volumes/caca/test_fractionation/01figure/figure6/yield_CGA_GO.xlsx"
# id = pd.read_excel(id_file, sheet_name="Sheet1")["ID"]
# GO = pd.read_csv(GO_file, dtype=str)
# filter_data = GO[GO.iloc[:, 0].isin(id)]
# filter_data.to_excel(output_file, index=False)

# import pandas as pd
# input_file = "/Volumes/caca/test_fractionation/01figure/figure6/yield_CGA_GO.xlsx"
# df = pd.read_excel(input_file, sheet_name="Sheet1")
# counts = df.iloc[:, 8].value_counts()
# print(counts)

import pandas as pd
input_file = '/Volumes/caca/test_fractionation/01figure/figure6/rubber_cga_info.xlsx'
df = pd.read_excel(input_file, sheet_name='Sheet1')
gene_ncp_map = df.groupby('Gene')['NCP'].apply(list).reset_index()
max_ncp_count = gene_ncp_map['NCP'].apply(len).max()
for i in range(max_ncp_count):
    gene_ncp_map[f'NCP_{i+1}'] = gene_ncp_map['NCP'].apply(lambda x: x[i] if i < len(x) else None)
gene_ncp_map_final = gene_ncp_map.drop('NCP', axis=1)
output_file = '/Volumes/caca/test_fractionation/01figure/figure6/rubber_cga_info_veen.xlsx'
gene_ncp_map_final.to_excel(output_file, index=False)


# import os
# import numpy as np
# from collections import defaultdict
# from pathlib import Path
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
# def main():
#     folder_path = r"G:\test\03WGCNA\WGCNA_Results_yield\module_specific_results"
#     id_file = r"G:\test\02figure\figure7\yield_common_id.txt"
#     merge_output = r"G:\test\02figure\figure7\yield_common_id_tom.txt"
#     merge_files_by_id(folder_path, id_file, merge_output)
# if __name__ == "__main__":
#     main()

# tom_file = r"G:\test\02figure\figure7\yield_common_id_tom.txt"
# sp_id = set()
# with open(tom_file, "r", encoding='utf-8') as tom:
#     for line in tom:
#         sp_id.add(line.split()[0])
# print(len(sp_id))

# # 合并SRR信息
# import pandas as pd
# def merge_csv_info_into_excel(excel_path, csv_path, output_path, sheet_name="Sheet2", csv_sep=",", csv_encoding="utf-8"):
#     df_x = pd.read_excel(excel_path, sheet_name=sheet_name, dtype=str)
#     df_x = df_x.copy()
#     df_x.rename(columns={df_x.columns[0]: "SRR"}, inplace=True)
#     df_x["SRR"] = df_x["SRR"].astype(str).str.strip()
#     df_c = pd.read_csv(csv_path, sep=csv_sep, encoding=csv_encoding, dtype=str)
#     df_c = df_c.copy()
#     df_c.rename(columns={df_c.columns[0]: "SRR"}, inplace=True)
#     df_c["SRR"] = df_c["SRR"].astype(str).str.strip()
#     dup_count = df_c["SRR"].duplicated(keep="first").sum()
#     if dup_count > 0:
#         print(f"CSV 中检测到 {dup_count} 个重复 SRR，已保留首次出现的记录。")
#         df_c = df_c.drop_duplicates(subset=["SRR"], keep="first")
#     merged = df_x.merge(df_c, on="SRR", how="left", suffixes=("", "_csv"))
#     merged.to_excel(output_path, index=False)
#     print(f"合并完成，已保存到：{output_path}")
# if __name__ == "__main__":
#     excel_path = r"D:\Desktop\test.xlsx"
#     csv_path   = r"D:\Desktop\SraRunInfo.csv"
#     output_path = r"D:\Desktop\merged_output.xlsx"
#     merge_csv_info_into_excel(excel_path, csv_path, output_path)

# # Table S7
# import pandas as pd
# def generate_id_presence_matrix(excel_ids_path, excel_accessions_path, excel_mapping_path, output_path):
#     # 读 Excel A
#     df_a = pd.read_excel(excel_ids_path, dtype=str)
#     df_a = df_a.apply(lambda col: col.str.strip())
#     df_a = df_a.fillna("")
#     # 为每列ID创建集合
#     id_sets = {col: set(df_a[col].dropna().astype(str).str.strip()) - {""}
#                for col in df_a.columns}
#     # 读 Excel B
#     df_b = pd.read_excel(excel_accessions_path, dtype=str, sheet_name="Sheet2")
#     accession_col = df_b.columns[0]
#     df_b[accession_col] = df_b[accession_col].str.strip()
#     # 读映射表
#     df_m = pd.read_excel(excel_mapping_path, dtype=str)
#     map_acc_col = df_m.columns[0]
#     map_id_col = df_m.columns[1]
#     df_m[map_acc_col] = df_m[map_acc_col].str.strip()
#     df_m[map_id_col] = df_m[map_id_col].str.strip()
#     # 合并 mapping 到 B
#     df_merged = df_b.merge(df_m, how="left",
#                            left_on=accession_col, right_on=map_id_col)
#     # 对 A 的每一列生成 0/1 列
#     for col_name, id_set in id_sets.items():
#         df_merged[col_name] = df_merged[map_acc_col].apply(
#             lambda x: 1 if x in id_set else 0
#         )
#     # 写出
#     df_merged.to_excel(output_path, index=False)
#     print(f"完成！结果已保存到 {output_path}")
# if __name__ == "__main__":
#     excel_ids_path = r"D:\Desktop\yield_veen_pre.xlsx"
#     excel_accessions_path = r"D:\Desktop\test.xlsx"
#     excel_mapping_path = r"D:\Desktop\NCP_accessions.xlsx"
#     output_path = r"D:\Desktop\B_with_matrix.xlsx"
#     generate_id_presence_matrix(
#         excel_ids_path,
#         excel_accessions_path,
#         excel_mapping_path,
#         output_path
#     )
