# import os
# import numpy as np
# from collections import defaultdict
# from pathlib import Path
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
#     filtered_genes = {
#         gene: count for gene, count in gene_counts.items()
#         if min_count <= count <= max_count
#     }
#     with open(output_data_file, 'w', encoding='utf-8') as f_data:
#         f_data.write("\n".join(
#             line for gene in filtered_genes
#             for line in gene_lines[gene]
#         ))
#     freq_dist = defaultdict(int)
#     for count in gene_counts.values():
#         if min_count <= count <= max_count:
#             freq_dist[count] += 1
#     with open(output_freq_file, 'w', encoding='utf-8') as f_freq:
#         f_freq.write("出现次数\t基因数量\n")
#         for count, num_genes in sorted(freq_dist.items()):
#             f_freq.write(f"{count}\t{num_genes}\n")
# def main():
#     folder_path = "/Volumes/caca/test_fractionation/01figure/WGCNA/WGCNA_Results_DT_ST/module_specific_results"
#     output_base = "/Users/lemon/Desktop"
#     novel_data_output = os.path.join(output_base, "dt_st_tom_all_id.txt")
#     unique_novel_data = extract_novel_data(folder_path)
#     with open(novel_data_output, 'w', encoding='utf-8') as out_file:
#         for item in unique_novel_data:
#             out_file.write(item + '\n')
#     tom_range_output = os.path.join(output_base, "dt_st_tom_range_statistics.txt")
#     count_tom_ranges(folder_path, tom_range_output)
#     gene_data_output = os.path.join(output_base, "dt_st_tom_all_info.txt")
#     gene_freq_output = os.path.join(output_base, "dt_st_tom_frequency.txt")
#     process_and_extract_genes(folder_path, gene_data_output, gene_freq_output, min_count=1, max_count=30000)
#     print("所有任务已完成！")
# if __name__ == "__main__":
#     main()

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
#     folder_path = r"G:\test\03WGCNA\WGCNA_Results_DT_ST\module_specific_results"
#     id_file = r"G:\test\02figure\figure6\ST_DT\st_common_id.txt"
#     merge_output = r"G:\test\02figure\figure6\ST_DT\st_common_id_tom.txt"
#     merge_files_by_id(folder_path, id_file, merge_output)
# if __name__ == "__main__":
#     main()

# tom_file = r"G:\test\02figure\figure6\ST_DT\dt_common_id_tom.txt"
# sp_id = set()
# with open(tom_file, "r", encoding='utf-8') as tom:
#     for line in tom:
#         sp_id.add(line.split()[0])
# print(len(sp_id))

import pandas as pd
express_info = r"G:\test_fractionation\express_info_115.xlsx"
output_file = r"G:\test_fractionation\down.xlsx"
df_info = pd.read_excel(express_info, sheet_name="Sheet1")
df_id = pd.read_excel(express_info, sheet_name="down")
filter_df = df_info[df_info['ID'].isin(df_id["GeneID"])]
filter_df.to_excel(output_file, index=False)