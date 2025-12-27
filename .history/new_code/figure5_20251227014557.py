# from datetime import datetime
# import pandas as pd
# def excel_to_gff3(df, output_gff):
#     required_columns = ['ID', 'start', 'end', 'strand', 'chrom']
#     missing_cols = [col for col in required_columns if col not in df.columns]
#     if missing_cols:
#         print(f"错误: 缺少必要的列: {', '.join(missing_cols)}")
#         return
#     gff_header = f"""##gff-version 3
# ##date {datetime.now().strftime('%Y-%m-%d')}
# ##source {df}
# ##genome-build v1.0
# """
#     gff_lines = []
#     for idx, row in df.iterrows():
#         seqid = row['chrom']
#         source = "EuNCP"
#         start = int(row['start'])
#         end = int(row['end'])
#         strand = row['strand']
#         gene_id = row['ID']
#         gene_line = (
#             f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={gene_id};Name={gene_id}"
#         )
#         gff_lines.append(gene_line)
#         mrna_id = f"{gene_id}.t1"
#         mrna_line = (
#             f"{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={mrna_id};Parent={gene_id};product=predicted protein"
#         )
#         gff_lines.append(mrna_line)
#         exon_id = f"{mrna_id}.exon1"
#         exon_line = (
#             f"{seqid}\t{source}\texon\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={exon_id};Parent={mrna_id}"
#         )
#         gff_lines.append(exon_line)
#         cds_id = f"{mrna_id}.cds"
#         cds_line = (
#             f"{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
#             f"ID={cds_id};Parent={mrna_id}"
#         )
#         gff_lines.append(cds_line)
#     try:
#         with open(output_gff, 'w') as f:
#             f.write(gff_header)
#             f.write("\n".join(gff_lines))
#         print(f"成功生成GFF3文件: {output_gff}")
#         print(f"转换了 {len(df)} 条基因记录，共生成 {len(gff_lines)} 行GFF记录")
#     except Exception as e:
#         print(f"写入GFF文件失败: {e}")
# excel_file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\output_best.xlsx"
# output_gff = r"D:\Desktop\peptidemicro\00file\01figure\figure5\NCP_codon.gff"
# df = pd.read_excel(excel_file, sheet_name="codon_1.0")
# excel_to_gff3(df, output_gff)

# # 合并文件
# import pandas as pd
# import os
# from tqdm import tqdm
# def merge_count_files(input_dir, RNA_info_file, output_file, gene_id_col="Geneid"):
#     count_files = pd.read_excel(RNA_info_file, sheet_name="group")
#     merged_df = None
#     for _, row in tqdm(count_files.iterrows(), desc="合并进度"):
#         file = row['Sample']
#         sample_name = f"{file}_counts.txt"
#         file_path = os.path.join(input_dir, sample_name)
#         try:
#             df = pd.read_csv(file_path, sep='\t', comment='#')
#             counts = df[[gene_id_col, df.columns[-1]]]
#             counts.columns = ['GeneID', file]
#             if merged_df is None:
#                 merged_df = counts
#             else:
#                 merged_df = pd.merge(merged_df, counts, on='GeneID', how='outer')
#         except Exception as e:
#             print(f"\n处理失败 {file}: {str(e)}")
#             continue
#     if merged_df is not None:
#         print(f"\n合并后数据维度: {merged_df.shape}")
#         if merged_df.duplicated('GeneID').any():
#             print(f"警告：存在重复基因ID，将取第一个出现的值")
#             merged_df = merged_df.drop_duplicates('GeneID')
#         merged_df.to_csv(output_file, index=False)
#         return merged_df
#     else:
#         raise ValueError(f"错误：未成功合并任何数据")
# def main():
#     count_dir = r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\02count_gene"
#     RNA_info_file = r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\Total_rna_seq.xlsx"
#     output_file = r"D:\Desktop\peptidemicro\00file\01figure\figure5\total_gene_count_matrix.csv"
#     merge_count_files(count_dir, RNA_info_file, output_file)
# if __name__ == "__main__":
#     main()

# 提取特定基因表达量
import pandas as pd
id_mapping_df = pd.read_excel("/Users/lemon/Desktop/rubber.xlsx", sheet_name="Sheet2")
data_df = pd.read_excel("/Users/lemon/Desktop/rubber_deseq2.xlsx", sheet_name="all_deseq")
non_evm_df = data_df[~data_df['ID'].str.startswith('evm')]
mapped_ids = id_mapping_df['ID']
mapped_df = data_df[data_df['name'].isin(mapped_ids)]
final_df = pd.concat([non_evm_df, mapped_df])
final_df.to_excel("filtered_data.xlsx", index=False)
print(final_df)
