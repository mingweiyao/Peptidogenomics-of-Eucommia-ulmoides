# 转录组数据基础的NCP分布veen图
# figure4b：不同组织中肽表达数量
import pandas as pd
gene_expression_df = pd.read_csv(r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\03ouput\finally_ncps_filtered_cpm.xlsx", index_col=0)
tissue_mapping_df = pd.read_excel(r"F:\Eu_peptido\00file\00raw\rnaseq\Total_rna_seq.xlsx", sheet_name="Sheet4")
tissue_to_samples = tissue_mapping_df.groupby('Tissues')['Sample'].apply(list).to_dict()
peptide_ids_by_tissue = {}
for tissue, samples in tissue_to_samples.items():
    tissue_expr = gene_expression_df[samples]
    expressed_peptides = tissue_expr[(tissue_expr >= 1).any(axis=1)].index.tolist()
    peptide_ids_by_tissue[tissue] = expressed_peptides
max_length = max(len(ids) for ids in peptide_ids_by_tissue.values())
peptide_ids_df = pd.DataFrame({
    tissue: ids + [None]*(max_length - len(ids)) 
    for tissue, ids in peptide_ids_by_tissue.items()
})
output_csv_path = r"F:\Eu_peptido\new_prepare\figure4_peptide_expression_ids.csv"
peptide_ids_df.to_csv(output_csv_path, index=False)
print("✅ 结果已保存至:", output_csv_path)
print("\n示例数据：")
print(peptide_ids_df.head())