# 根据group提取并筛选表达的sp
import pandas as pd
gene_expression_df = pd.read_csv(r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_cpm.csv", index_col=0)
tissue_mapping_df = pd.read_excel(r"D:\Desktop\Total_rna_seq.xlsx", sheet_name="extract")
tissue_group_to_samples = (tissue_mapping_df.groupby(['Tissues', 'Group'])['Sample'].apply(list).to_dict())
peptide_ids_by_tissue = {}
for (tissue, group), samples in tissue_group_to_samples.items():
    samples = [s for s in samples if s in gene_expression_df.columns]
    group_expr = gene_expression_df.loc[:, samples]
    group_mask = (group_expr >= 1).all(axis=1)
    if tissue not in peptide_ids_by_tissue:
        peptide_ids_by_tissue[tissue] = set()
    peptide_ids_by_tissue[tissue].update(group_expr.index[group_mask].tolist())
peptide_ids_by_tissue = {tissue: list(ids) for tissue, ids in peptide_ids_by_tissue.items()}
max_length = max((len(ids) for ids in peptide_ids_by_tissue.values()), default=0)
peptide_ids_df = pd.DataFrame({
    tissue: ids + [None] * (max_length - len(ids))
    for tissue, ids in peptide_ids_by_tissue.items()
})
output_excel_path = r"D:\Desktop\rnaseq_veen_leaf_location.xlsx"
peptide_ids_df.to_excel(output_excel_path, index=False)