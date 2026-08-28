import os
import numpy as np
from collections import defaultdict
from pathlib import Path
# 1. 提取高TOM值的NCP ID并汇总
def extract_novel_data(folder_path):
    novel_data = set()
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(folder_path, filename)
            try:
                with open(file_path, 'r', encoding='utf-8') as file:
                    for line in file:
                        columns = line.strip().split()
                        novel_data.add(columns[0])
            except Exception as e:
                print(f"处理文件 {filename} 时出错: {e}")
    return sorted(novel_data)
# 2. 统计TOM值分布并输出统计区间
def count_tom_ranges(folder_path, output_file):
    all_tom_values = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r') as file:
                next(file)
                for line in file:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        try:
                            tom = float(parts[2])
                            all_tom_values.append(tom)
                        except ValueError:
                            continue
    if not all_tom_values:
        print("警告：未找到有效的TOM值")
        return
    min_tom = min(all_tom_values)
    max_tom = max(all_tom_values)
    bins = np.linspace(min_tom, max_tom, 11)
    range_counts = defaultdict(int)
    # 对每个文件再次统计TOM值所在区间的分布
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r') as file:
                next(file)
                for line in file:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        try:
                            tom = float(parts[2])
                            for i in range(len(bins) - 1):
                                if bins[i] <= tom < bins[i + 1]:
                                    range_counts[(bins[i], bins[i + 1])] += 1
                                    break
                            if tom == bins[-1]:
                                range_counts[(bins[-2], bins[-1])] += 1
                        except ValueError:
                            continue
    with open(output_file, 'w') as out_file:
        out_file.write("TOM区间\t数据行数\n")
        sorted_ranges = sorted(range_counts.items(), key=lambda x: x[0][0])
        for (start, end), count in sorted_ranges:
            out_file.write(f"[{start:.4f}, {end:.4f})\t{count}\n")
    print(f"TOM值区间统计结果已保存到 {output_file}")
    print(f"TOM值范围: {min_tom:.4f} 到 {max_tom:.4f}")
# 3. 统计肽连接基因数分布频率
def process_and_extract_genes(input_folder, output_data_file, output_freq_file, min_count, max_count):
    gene_counts = defaultdict(int)
    gene_lines = defaultdict(list)
    for filename in os.listdir(input_folder):
        if filename.endswith('.txt'):
            file_path = os.path.join(input_folder, filename)
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    next(f)  # 跳过第一行
                    for line in f:
                        line = line.strip()
                        if line:
                            gene = line.split('\t')[0]
                            gene_counts[gene] += 1
                            gene_lines[gene].append(line)
            except Exception as e:
                print(f"处理文件 {filename} 时出错: {e}")
    # 筛选出现次数在指定范围内的基因
    filtered_genes = {
        gene: count for gene, count in gene_counts.items()
        if min_count <= count <= max_count
    }
    # 输出筛选后的数据
    with open(output_data_file, 'w', encoding='utf-8') as f_data:
        f_data.write("\n".join(
            line for gene in filtered_genes
            for line in gene_lines[gene]
        ))
    # 频率分布统计
    freq_dist = defaultdict(int)
    for count in gene_counts.values():
        if min_count <= count <= max_count:
            freq_dist[count] += 1
    # 输出频率统计结果
    with open(output_freq_file, 'w', encoding='utf-8') as f_freq:
        f_freq.write("出现次数\t基因数量\n")
        for count, num_genes in sorted(freq_dist.items()):
            f_freq.write(f"{count}\t{num_genes}\n")
# 4. 提取并合并共有NCP基因对信息
def merge_files_by_id(folder_path, id_file_path, output_file_path):
    with open(id_file_path, 'r', encoding='utf-8') as id_file:
        target_ids = set(line.strip() for line in id_file if line.strip())
    with open(output_file_path, 'w', encoding='utf-8') as output_file:
        for filename in os.listdir(folder_path):
            if filename.endswith('.txt') and filename != os.path.basename(id_file_path):
                file_path = os.path.join(folder_path, filename)
                try:
                    with open(file_path, 'r', encoding='utf-8') as input_file:
                        for line in input_file:
                            line = line.strip()
                            if line:
                                first_col = line.split()[0]
                                if first_col in target_ids:
                                    output_file.write(line + '\n')
                except Exception as e:
                    print(f"处理文件 {filename} 时出错: {e}")
# 5. 执行所有功能
def main():
    folder_path = r"G:\peptidegenomics\test\WGCNA_Results_wt\module_specific_results"
    output_base = r"D:\Desktop"
    # 汇总高TOM值的NCPs
    novel_data_output = os.path.join(output_base, "wt_tom_all_id.txt")
    unique_novel_data = extract_novel_data(folder_path)
    with open(novel_data_output, 'w', encoding='utf-8') as out_file:
        for item in unique_novel_data:
            out_file.write(item + '\n')
    # TOM值分布统计
    tom_range_output = os.path.join(output_base, "wt_tom_range_statistics.txt")
    count_tom_ranges(folder_path, tom_range_output)
    # 统计基因连接频率
    gene_data_output = os.path.join(output_base, "wt_tom_all_info.txt")
    gene_freq_output = os.path.join(output_base, "wt_tom_frequency.txt")
    process_and_extract_genes(folder_path, gene_data_output, gene_freq_output, min_count=1, max_count=30000)
    # 提取并合并共有的NCP基因对信息
    id_file = r"G:\peptidegenomics\test\wt\common_id.txt"
    merge_output = os.path.join(output_base, "common_id_wt_tom.txt")
    merge_files_by_id(folder_path, id_file, merge_output)
    print("所有任务已完成！")
if __name__ == "__main__":
    main()
