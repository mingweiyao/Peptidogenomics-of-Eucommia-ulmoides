# # 根据group提取并筛选表达的sp
# import pandas as pd
# gene_expression_df = pd.read_csv(r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_cpm.csv", index_col=0)
# tissue_mapping_df = pd.read_excel(r"D:\Desktop\Total_rna_seq.xlsx", sheet_name="extract")
# tissue_group_to_samples = (tissue_mapping_df.groupby(['Tissues', 'Group'])['Sample'].apply(list).to_dict())
# peptide_ids_by_tissue = {}
# for (tissue, group), samples in tissue_group_to_samples.items():
#     samples = [s for s in samples if s in gene_expression_df.columns]
#     group_expr = gene_expression_df.loc[:, samples]
#     group_mask = (group_expr >= 1).all(axis=1)
#     if tissue not in peptide_ids_by_tissue:
#         peptide_ids_by_tissue[tissue] = set()
#     peptide_ids_by_tissue[tissue].update(group_expr.index[group_mask].tolist())
# peptide_ids_by_tissue = {tissue: list(ids) for tissue, ids in peptide_ids_by_tissue.items()}
# max_length = max((len(ids) for ids in peptide_ids_by_tissue.values()), default=0)
# peptide_ids_df = pd.DataFrame({
#     tissue: ids + [None] * (max_length - len(ids))
#     for tissue, ids in peptide_ids_by_tissue.items()
# })
# output_excel_path = r"D:\Desktop\rnaseq_veen_leaf_location.xlsx"
# peptide_ids_df.to_excel(output_excel_path, index=False)

# # 转录组数据基础的NCP分布veen图
# import pandas as pd
# gene_expression_df = pd.read_csv(r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_cpm.csv", index_col=0)
# tissue_mapping_df = pd.read_excel(r"D:\Desktop\Total_rna_seq.xlsx", sheet_name="extract")
# tissue_group_to_samples = (tissue_mapping_df.groupby(['Tissues', 'Group'])['Sample'].apply(list).to_dict())
# peptide_ids_by_tissue = {}
# for (tissue, group), samples in tissue_group_to_samples.items():
#     samples = [s for s in samples if s in gene_expression_df.columns]
#     group_expr = gene_expression_df.loc[:, samples]
#     group_mask = (group_expr >= 1).all(axis=1)
#     if tissue not in peptide_ids_by_tissue:
#         peptide_ids_by_tissue[tissue] = set()
#     peptide_ids_by_tissue[tissue].update(group_expr.index[group_mask].tolist())
# peptide_ids_by_tissue = {tissue: list(ids) for tissue, ids in peptide_ids_by_tissue.items()}
# max_length = max((len(ids) for ids in peptide_ids_by_tissue.values()), default=0)
# peptide_ids_df = pd.DataFrame({
#     tissue: ids + [None] * (max_length - len(ids))
#     for tissue, ids in peptide_ids_by_tissue.items()
# })
# output_excel_path = r"D:\Desktop\rnaseq_veen_leaf_location.xlsx"
# peptide_ids_df.to_excel(output_excel_path, index=False)

# # 数量-肽统计
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# from brokenaxes import brokenaxes
# import warnings
# warnings.filterwarnings('ignore')
# def plot_peptide_distribution_bar_broken_y(expression_matrix_path, output_pdf_path, output_csv_path, threshold=1, ylims=((0, 400), (600, 800), (2800, 3000)), bar_color="#4C72B0", edge_color="white"):
#     if expression_matrix_path.endswith('.csv'):
#         df = pd.read_csv(expression_matrix_path, index_col=0)
#     elif expression_matrix_path.endswith(('.xls', '.xlsx')):
#         df = pd.read_excel(expression_matrix_path, index_col=0)
#     else:
#         raise ValueError("仅支持 csv / xlsx")
#     print(f"读取矩阵：{df.shape[0]} 个ID × {df.shape[1]} 个列")
#     # ---------- 2. 二值化 & 统计 ----------
#     binary_df = (df >= threshold).astype(int)
#     counts_per_id = binary_df.sum(axis=1)
#     count_dist = counts_per_id.value_counts().sort_index()
#     total_ids = len(counts_per_id)
#     plot_data = pd.DataFrame({
#         "Number_of_Columns": count_dist.index.astype(int),
#         "Number_of_IDs": count_dist.values.astype(int),
#         "Percentage": (count_dist.values / total_ids * 100).round(2)
#     })
#     # ---------- 3) 输出 CSV ----------
#     plot_data.to_csv(output_csv_path, index=False)
#     print(f"✅ 统计结果已保存为：{output_csv_path}")
#     # ---------- 4) 画柱状图（broken y-axis） ----------
#     fig = plt.figure(figsize=(10, 8))
#     bax = brokenaxes(ylims=ylims, hspace=0.04, fig=fig)
#     x = plot_data["Number_of_Columns"].values
#     y = plot_data["Number_of_IDs"].values
#     bax.bar(x, y, width=0.8, color=bar_color, edgecolor=edge_color, linewidth=1.0)
#     bax.set_xlabel(f"Number of columns with CPM ≥ {threshold}")
#     bax.set_ylabel("Count of IDs")
#     bax.set_title("Peptide/ID Expression Distribution (Y-axis truncated)")
#     fig.savefig(output_pdf_path, bbox_inches="tight")
#     plt.close(fig)
#     print(f"✅ 截断柱状图已保存为：{output_pdf_path}")
#     return plot_data
# if __name__ == "__main__":
#     expression_file = r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_cpm.csv"
#     output_pdf = r"F:\work_mechanism\new_file\02figure\figure4\test_peptide_distribution_bar_truncated.pdf"
#     output_csv = r"F:\work_mechanism\new_file\02figure\figure4\test_peptide_distribution_stat.csv"
#     plot_peptide_distribution_bar_broken_y(
#         expression_file,
#         output_pdf,
#         output_csv,
#         threshold=1,
#         ylims=((0, 400), (600, 800), (2800, 3000))
#     )
