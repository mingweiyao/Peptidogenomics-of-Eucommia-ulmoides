# # 如果在某个样本中表达，那么筛选在该样本所在的所有重复中都表达的NCP
# import pandas as pd
# expression_csv = r"G:\Eu_peptido\20251113 horticulture research\file\sp_expressed.csv"
# group_excel = r"G:\Eu_peptido\20251113 horticulture research\file\Total_rna_seq.xlsx"
# group_sheet_name = "group"
# output_excel = r"G:\Eu_peptido\20251113 horticulture research\file\expression_filtered.xlsx"
# expr = pd.read_csv(expression_csv, index_col=0)
# group_df = pd.read_excel(group_excel, sheet_name=group_sheet_name)
# srr_col = group_df.columns[0]
# group_col = group_df.columns[2]
# group_df = group_df[group_df[srr_col].isin(expr.columns)]
# group_to_srrs = (
#     group_df
#     .groupby(group_col)[srr_col]
#     .apply(list)
#     .to_dict()
# )
# row_keep_mask = pd.Series(False, index=expr.index)
# for grp, srr_list in group_to_srrs.items():
#     srr_in_expr = [s for s in srr_list if s in expr.columns]
#     if not srr_in_expr:
#         continue
#     grp_mask = (expr[srr_in_expr] >= 10).all(axis=1)
#     row_keep_mask = row_keep_mask | grp_mask
# expr_filtered = expr[row_keep_mask]
# expr_filtered.to_excel(output_excel, index_label="Gene")
# print(f"筛选完成！共保留 {expr_filtered.shape[0]} 行，已保存到: {output_excel}")

from datetime import datetime
import pandas as pd

def excel_to_gff3(df, output_gff):
    required_columns = ['ID', 'start', 'end', 'strand', 'chrom']
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        print(f"错误: 缺少必要的列: {', '.join(missing_cols)}")
        return
    gff_header = f"""##gff-version 3
##date {datetime.now().strftime('%Y-%m-%d')}
##source {df}
##genome-build v1.0
"""
    gff_lines = []
    for idx, row in df.iterrows():
        seqid = row['chrom']
        source = "EuNCP"
        start = int(row['start'])
        end = int(row['end'])
        strand = row['strand']
        gene_id = row['ID']
        gene_line = (
            f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={gene_id};Name={gene_id}"
        )
        gff_lines.append(gene_line)
        mrna_id = f"{gene_id}.t1"
        mrna_line = (
            f"{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={mrna_id};Parent={gene_id};product=predicted protein"
        )
        gff_lines.append(mrna_line)
        exon_id = f"{mrna_id}.exon1"
        exon_line = (
            f"{seqid}\t{source}\texon\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={exon_id};Parent={mrna_id}"
        )
        gff_lines.append(exon_line)
        cds_id = f"{mrna_id}.cds"
        cds_line = (
            f"{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
            f"ID={cds_id};Parent={mrna_id}"
        )
        gff_lines.append(cds_line)
    try:
        with open(output_gff, 'w') as f:
            f.write(gff_header)
            f.write("\n".join(gff_lines))
        print(f"成功生成GFF3文件: {output_gff}")
        print(f"转换了 {len(df)} 条基因记录，共生成 {len(gff_lines)} 行GFF记录")
    except Exception as e:
        print(f"写入GFF文件失败: {e}")

excel_file = "D:/Desktop/sp_express_info.xlsx"
output_gff = "D:/Desktop/junction.gff"
df = pd.read_excel(excel_file, sheet_name="junction")
excel_to_gff3(df, output_gff)