# import pandas as pd
# from Bio import SeqIO
# df = pd.read_excel('/Volumes/caca/test_fractionation/NCP_with_M.xlsx')
# genome = SeqIO.to_dict(SeqIO.parse('/Volumes/caca/test_fractionation/00raw/Eu_genome.fasta', 'fasta'))
# sequences = []
# for index, row in df.iterrows():
#     chrom = str(row['chrom'])
#     start = int(row['start'])
#     end = int(row['end'])
#     if chrom in genome:
#         sequence = genome[chrom].seq[start-1:end]
#         sequences.append(str(sequence))
#     else:
#         sequences.append(None)
# df['sequence'] = sequences
# df.to_excel('/Volumes/caca/test_fractionation/output.xlsx', index=False)
# print("序列提取完成！")

# import os
# import re
# import pandas as pd
# # ========= 配置区（按需修改） =========
# corr_file = "/Volumes/caca/test_fractionation/01figure/figure6/rubber_CGA_NCP.xlsx"
# corr_sheet = "pearson_filter"
# out_file = "/Volumes/caca/test_fractionation/01figure/figure6/heatmaps.xlsx"
# abs_cutoff = 0.8      
# # ====================================
# corr_df = pd.read_excel(corr_file, sheet_name=corr_sheet)
# num_cols = corr_df.columns[1:]
# expr_df = pd.read_excel(corr_file, sheet_name="Sheet1")
# expr_df["ID"] = expr_df["ID"].astype(str)
# with pd.ExcelWriter(out_file, engine="xlsxwriter") as writer:
#     for col in num_cols:
#         mask = corr_df[col].abs() >= abs_cutoff
#         ids = corr_df.loc[mask, "ID"].dropna().astype(str).drop_duplicates()
#         filter_expr = expr_df[expr_df['ID'].isin(ids)]
#         filter_expr.to_excel(writer, sheet_name=f"{col}", index=False)
# print(f"完成：已按列将表达量子矩阵写入多sheet → {out_file}")

# # CPM标准化
# import pandas as pd
# import numpy as np
# import os
# in_excel = "/Volumes/caca/test_fractionation/01figure/figure6/heatmap_counts/heatmap_count.xlsx"
# id_col_name = "GeneID"
# out_excel = "/Volumes/caca/test_fractionation/01figure/figure6/heatmap_counts/heatmap_count_CPM.xlsx"
# write_log2 = True
# df = pd.read_excel(in_excel, sheet_name="Sheet1")
# id_col = df[[id_col_name]].copy()
# expr = df.drop(columns=[id_col_name]).copy()
# expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)
# libsize = expr.sum(axis=0)
# denom = libsize.replace(0, np.nan)
# cpm = expr.div(denom, axis=1) * 1e6
# cpm = cpm.fillna(0)
# cpm_out = pd.concat([id_col, cpm], axis=1)
# if write_log2:
#     log2cpm = np.log2(cpm + 1.0)
#     log2cpm_out = pd.concat([id_col, log2cpm], axis=1)
# with pd.ExcelWriter(out_excel, engine="xlsxwriter") as writer:
#     cpm_out.to_excel(writer, sheet_name="CPM", index=False)
#     if write_log2:
#         log2cpm_out.to_excel(writer, sheet_name="log2CPM", index=False)
# print(f"[完成] 已写出：{out_excel}")
# print("样本总reads（library size）概览：")
# print(libsize.describe())

# # 提取数据
# import pandas as pd
# input_file = "/Volumes/caca/test_fractionation/01figure/figure6/heatmap_counts/Eu_genome_quan.xlsx"
# output_file = "/Volumes/caca/test_fractionation/01figure/figure6/heatmap_counts/rubber_CGA_CPM_log.xlsx"
# df_id = pd.read_excel(input_file, sheet_name="Sheet2",usecols=["ID"])
# df_express = pd.read_excel(input_file, sheet_name="Sheet1")
# filtered_df = df_express[df_express['GeneID'].isin(df_id['ID'])]
# filtered_df.to_excel(output_file, index=False)


# import pandas as pd
# # ===== 配置区 =====
# # r/p 矩阵
# input_file    = "/Volumes/caca/test_fractionation/01figure/figure6/heatmap_counts/pearson.xlsx"
# r_sheet       = "R_Value"   # r 矩阵
# p_sheet       = "P_Value"   # p 矩阵
# # 表达量矩阵（第一列为 ID）
# expr_file     = "/Volumes/caca/test_fractionation/01figure/figure6/heatmap_counts/Eu_genome_quan.xlsx"
# expr_sheet    = "Sheet1"    # 表达量所在表
# expr_id_col   = "GeneID"        # 表达量文件中 ID 列名
# # ——新增：列名→第二名称 的映射（放在表达量 Excel 的另一个 sheet）——
# map_sheet     = "Sheet3"          # 存放映射关系的 sheet 名
# map_key_col   = "primary_name"      # 映射表中“原列名”的列名（= r/p 的列名）
# map_value_col = "secondary_name"    # 映射表中“第二名称（别名）”的列名
# # 输出
# output_file   = "/Volumes/caca/test_fractionation/01figure/figure6/heatmap_counts/rubber_CGA_NCP_filter_0.5.xlsx"
# # 阈值
# r_threshold   = 0.8
# p_threshold   = 0.05
# # 若首行是额外标题/说明，需要丢弃
# drop_first_row = True
# # ==================
# # 读入 r/p
# r_df_raw = pd.read_excel(input_file, sheet_name=r_sheet)
# p_df_raw = pd.read_excel(input_file, sheet_name=p_sheet)
# if drop_first_row:
#     r_df_raw = r_df_raw.iloc[1:].reset_index(drop=True)
#     p_df_raw = p_df_raw.iloc[1:].reset_index(drop=True)
# # 默认第一列为 ID
# id_col = r_df_raw.columns[0]
# # 对齐 ID 与列
# id_intersect = pd.Index(r_df_raw[id_col].astype(str)).intersection(
#     pd.Index(p_df_raw[id_col].astype(str))
# )
# r_sub = r_df_raw[r_df_raw[id_col].astype(str).isin(id_intersect)].copy()
# p_sub = p_df_raw[p_df_raw[id_col].astype(str).isin(id_intersect)].copy()
# r_sub = r_sub.set_index(id_col).loc[id_intersect]
# p_sub = p_sub.set_index(id_col).loc[id_intersect]
# col_intersect = r_sub.columns.intersection(p_sub.columns)
# r_mat = r_sub[col_intersect].apply(pd.to_numeric, errors="coerce")
# p_mat = p_sub[col_intersect].apply(pd.to_numeric, errors="coerce")
# # 读入表达量矩阵与映射
# expr_df = pd.read_excel(expr_file, sheet_name=expr_sheet)
# expr_df[expr_id_col] = expr_df[expr_id_col].astype(str)
# alias_df = pd.read_excel(expr_file, sheet_name=map_sheet)
# alias_df[map_key_col] = alias_df[map_key_col].astype(str)
# alias_df[map_value_col] = alias_df[map_value_col].astype(str)
# name_map = dict(zip(alias_df[map_key_col], alias_df[map_value_col]))  # {原列名: 第二名称}
# summary = []
# with pd.ExcelWriter(output_file, engine="xlsxwriter") as writer:
#     for col in col_intersect:
#         # 该列的筛选（|r|>=阈值 且 p<阈值）
#         mask = r_mat[col].abs().ge(r_threshold) & p_mat[col].lt(p_threshold)
#         ids = r_mat.index[mask].astype(str)
#         # 目标行（sheet 对应的列名本身）
#         target_row = expr_df[expr_df[expr_id_col] == str(col)].copy()
#         if not target_row.empty:
#             alias_name = name_map.get(str(col), str(col))
#             target_row.loc[:, expr_id_col] = alias_name  # 仅替换目标行 ID
#         # 筛选到的 ID 的表达量
#         expr_sel = expr_df[expr_df[expr_id_col].isin(ids)].copy()
#         # 合并：目标行置顶，再接筛选到的 ID（去重）
#         combined = pd.concat(
#             [target_row, expr_sel[~expr_sel[expr_id_col].isin(target_row[expr_id_col])]],
#             axis=0
#         )
#         # ===== 行过滤：所有数值列都 < 0.5 的行，在写出前删除 =====
#         if not combined.empty:
#             num_cols = [c for c in combined.columns if c != expr_id_col]
#             vals = combined[num_cols].apply(pd.to_numeric, errors="coerce")
#             # NaN 视为“未知”，不计入“全都<0.5”
#             drop_mask = vals.fillna(float('inf')).lt(0.5).all(axis=1)
#             combined = combined.loc[~drop_mask].copy()
#         # 写出（sheet 名保持原列名）
#         sheet_name = col if len(col) <= 31 else col[:28] + "..."
#         if combined.empty:
#             # 仍写出空表头，便于检查
#             pd.DataFrame(columns=expr_df.columns).to_excel(writer, sheet_name=sheet_name, index=False)
#         else:
#             ordered_cols = [expr_id_col] + [c for c in combined.columns if c != expr_id_col]
#             combined[ordered_cols].to_excel(writer, sheet_name=sheet_name, index=False)
#         summary.append({
#             "Column": col,
#             "Alias_used": name_map.get(str(col), str(col)),
#             "Target_row_present": bool(len(target_row) > 0),
#             "Passed_ID_count_before_filter": int(len(ids) + (1 if not target_row.empty else 0)),
#             "Dropped_all<0.5_rows": int(drop_mask.sum()) if 'drop_mask' in locals() and not combined.empty else 0,
#             "Rows_after_filter": int(len(combined)) if not combined.empty else 0
#         })
#     pd.DataFrame(summary).sort_values("Rows_after_filter", ascending=False)\
#         .to_excel(writer, sheet_name="summary", index=False)
# print(f"完成：每列已输出（目标行改用别名）+（筛选ID表达量），并在写出前移除了“全列<0.5”的行 → {output_file}")


import pandas as pd
from intervaltree import Interval, IntervalTree
input_file = "/Volumes/caca/test_fractionation/01figure/figure3/portion_chromo.xlsx"
NCP_info = pd.read_excel(input_file, sheet_name="NCP_validated")
genome_info = pd.read_excel(input_file, sheet_name="genome")
def calculate_total_length(df):
    chrom_lengths = {}
    for chrom, group in df.groupby('chrom'):
        sorted_intervals = sorted(group[['start', 'end']].values, key=lambda x: x[0])
        merged_intervals = []
        current_start, current_end = sorted_intervals[0]
        for start, end in sorted_intervals[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged_intervals.append((current_start, current_end))
                current_start, current_end = start, end
        merged_intervals.append((current_start, current_end))
        total_length = sum(end - start for start, end in merged_intervals)
        chrom_lengths[chrom] = total_length
    return pd.DataFrame(list(chrom_lengths.items()), columns=['chrom', 'total_length'])
NCP_length = calculate_total_length(NCP_info)
genome_length = calculate_total_length(genome_info)
final_length = pd.merge(NCP_length, genome_length, on="chrom", how="outer", suffixes=('_NCP', '_genome'))
final_length.to_csv("/Volumes/caca/test_fractionation/01figure/figure3/chrom_total_lengths.csv", index=False)
print(final_length)
