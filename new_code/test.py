import pandas as pd
import numpy as np
input_file = r"F:\work_mechanism\new_file\02figure\figure4\SMD_analysis.xlsx"
# Step 1: 读取文件
id_df = pd.read_excel(input_file, sheet_name="ID")
group_df = pd.read_excel(input_file, sheet_name="Group")
quant_matrix_df = pd.read_excel(input_file, sheet_name="CPM")
# Step 2: 提取表达量数据
filter_df = quant_matrix_df[quant_matrix_df['ID'].isin(id_df['ID'])]
# Step 3: 提取分组数据
group = {}
for _, row in group_df.iterrows():
    if row['Group'] not in group:
        group[row['Group']] = []
    group[row['Group']].append(row['Sample'])
# Step 4: 计算每组的均值和标准差
group_stats = {}
for group_name, samples in group.items():
    group_expr = filter_df[samples]
    group_mean = group_expr.mean(axis=1)
    group_std = group_expr.std(axis=1)
    group_stats[group_name] = (group_mean, group_std)
# Step 5: 计算标准化差异（SMD）并追加到原始数据中
temp_group1 = []
for group1 in group_stats:
    temp_group1.append(group1)
    for group2 in group_stats:
        if group2 not in temp_group1:
            mean_diff = abs(group_stats[group1][0] - group_stats[group2][0])
            pooled_std = np.sqrt((group_stats[group1][1]**2 + group_stats[group2][1]**2) / 2)
            standardized_diff = mean_diff / pooled_std
            column_name = f'SMD_{group1}_vs_{group2}'
            filter_df[column_name] = standardized_diff
output_file = r"F:\work_mechanism\new_file\02figure\figure4\SMD_output.xlsx"
filter_df.to_excel(output_file, index=False)