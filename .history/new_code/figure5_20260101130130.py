# # 筛选肽基因与基因完全不重叠的候选肽基因
# import pandas as pd
# from collections import defaultdict
# def intervals_overlap_closed(a_start, a_end, b_start, b_end):
#     return max(a_start, b_start) <= min(a_end, b_end)
# def build_gene_buckets(gene_df):
#     buckets = defaultdict(list)
#     for _, g in gene_df.iterrows():
#         chrom = str(g["chrom"])
#         s = int(g["start"])
#         e = int(g["end"])
#         if s > e:
#             s, e = e, s
#         buckets[chrom].append((s, e))
#     return buckets
# def row_has_overlap_with_any_gene(row, gene_buckets):
#     chrom = str(row["chrom"])
#     s = int(row["phy_start"])
#     e = int(row["phy_end"])
#     if s > e:
#         s, e = e, s
#     intervals = gene_buckets.get((chrom), [])
#     for gs, ge in intervals:
#         if intervals_overlap_closed(s, e, gs, ge):
#             return True
#     return False
# def main(excel_path, sheet1_name, sheet2_name, out_path):
#     gene_df = pd.read_excel(excel_path, sheet_name=sheet1_name, engine="openpyxl")
#     s2 = pd.read_excel(excel_path, sheet_name=sheet2_name, engine="openpyxl")
#     gene_buckets = build_gene_buckets(gene_df)
#     s2 = s2.copy()
#     s2["_overlap_with_gene"] = s2.apply(lambda r: row_has_overlap_with_any_gene(r, gene_buckets), axis=1)
#     s2["_no_overlap_with_gene"] = ~s2["_overlap_with_gene"]
#     g = s2.groupby("accession")["_no_overlap_with_gene"] 
#     acc_all_no = set(g.all()[lambda x: x].index)
#     acc_any_no = set(g.any()[lambda x: x].index)
#     out_task1 = s2[s2["accession"].isin(acc_all_no)].drop(columns=["_overlap_with_gene", "_no_overlap_with_gene"])
#     out_task2 = s2[s2["accession"].isin(acc_any_no)].drop(columns=["_overlap_with_gene", "_no_overlap_with_gene"])
#     tmp = s2[s2["accession"].isin(acc_all_no)].copy()
#     tmp["total_score"] = pd.to_numeric(tmp["total_score"], errors="coerce")
#     best_rows = tmp.sort_values(["accession", "total_score"], ascending=[True, False]) \
#                    .drop_duplicates(subset=["accession"], keep="first") \
#                    .drop(columns=["_overlap_with_gene", "_no_overlap_with_gene"], errors="ignore")
#     with pd.ExcelWriter(out_path, engine="openpyxl") as w:
#         out_task1.to_excel(w, sheet_name="task1_all_no_overlap", index=False)
#         out_task2.to_excel(w, sheet_name="task2_any_no_overlap", index=False)
#         best_rows.to_excel(w, sheet_name="task1_best_total_score", index=False)
# if __name__ == "__main__":
#     main(
#         excel_path=r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\output_candidates.xlsx",
#         sheet1_name="gene_coor",
#         sheet2_name="output_candidates",
#         out_path=r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\output_candidates_filtered.xlsx",
#     )

# # 将Excel中的基因坐标转换为GFF3格式
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
# excel_file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\output_candidates_filtered.xlsx"
# output_gff = r"D:\Desktop\peptidemicro\00file\01figure\figure5\NCP_codon.gff"
# df = pd.read_excel(excel_file, sheet_name="task1_best_total_score")
# excel_to_gff3(df, output_gff)

# # 合并文件
# import pandas as pd
# import os
# from tqdm import tqdm
# def merge_count_files(input_dir, RNA_info_file, output_file, gene_id_col="Geneid"):
#     count_files = pd.read_excel(RNA_info_file, sheet_name="Sheet5")
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
#     count_dir = r"D:\Desktop\peptidemicro\00file\01figure\figure5\rnaseq"
#     RNA_info_file = r"D:\Desktop\peptidemicro\00file\01figure\Total_rna_seq.xlsx"
#     output_file = r"D:\Desktop\peptidemicro\00file\01figure\figure5\total_sp_count_matrix.csv"
#     merge_count_files(count_dir, RNA_info_file, output_file)
# if __name__ == "__main__":
#     main()

# # TPM标准化
# import os
# import pandas as pd
# import gffutils
# GFF_FILE = r"D:\Desktop\peptidemicro\00file\01figure\figure5\NCP_codon.gff"
# COUNT_FILE = r"D:\Desktop\peptidemicro\00file\01figure\total_all_matrix.xlsx"
# OUT_DIR = r"D:\Desktop\peptidemicro\00file\01figure"
# def prepare_length_data(gff_file):
#     db_path = gff_file + ".db"
#     if not os.path.exists(db_path):
#         print("🔄 正在创建 GFF 数据库...")
#         gffutils.create_db(
#             gff_file,
#             dbfn=db_path,
#             force=True,
#             keep_order=True,
#             merge_strategy="merge",
#             id_spec={"gene": "ID", "mRNA": "ID", "CDS": "Parent"},
#             disable_infer_genes=True,
#             disable_infer_transcripts=True
#         )
#     db = gffutils.FeatureDB(db_path)
#     gene_lengths = {}
#     for gene in db.features_of_type("gene"):
#         total_length = 0
#         for mrna in db.children(gene, featuretype="mRNA", order_by="start"):
#             exons = list(db.children(mrna, featuretype="exon", order_by="start"))
#             if exons:
#                 mrna_len = sum(e.end - e.start + 1 for e in exons)
#                 total_length += mrna_len
#         if total_length > 0:
#             gene_id = gene.id.replace("evm.model.", "evm.TU.")
#             gene_lengths[gene_id] = total_length
#     length_df = pd.DataFrame(
#         gene_lengths.items(),
#         columns=["GeneID", "length"]
#     )
#     out_len = os.path.join(OUT_DIR, "gene_lengths.csv")
#     length_df.to_csv(out_len, index=False)
#     print(f"✅ 基因长度表已生成：{out_len}")
#     return length_df
# def read_counts(count_file):
#     if count_file.lower().endswith(".csv"):
#         df = pd.read_csv(count_file)
#     else:
#         df = pd.read_excel(count_file)
#     if df.columns[0] != "GeneID":
#         df = df.rename(columns={df.columns[0]: "GeneID"})
#     return df
# def normalize_tpm(count_df, length_df, output_file):
#     df = pd.merge(count_df, length_df, on='GeneID', how='inner')
#     df = df[df['length'] > 0]
#     sample_cols = [col for col in df.columns if col not in ['GeneID', 'length']]
#     tpm_data = {}
#     for sample in sample_cols:
#         rpk = (df[sample] * 10**3) / df['length']
#         per_million_scaling_factor = rpk.sum() / 10**6
#         tpm = rpk / per_million_scaling_factor
#         tpm_data[sample] = tpm
#     tpm_df = pd.concat([df[['GeneID']], pd.DataFrame(tpm_data)], axis=1)
#     tpm_df.to_excel(output_file, index=False)
# if __name__ == "__main__":
#     length_df = prepare_length_data(GFF_FILE)
#     count_df = read_counts(COUNT_FILE)
#     normalize_tpm(count_df, length_df, os.path.join(OUT_DIR, "total_all_matrix_tpm.xlsx"))

# # 提取特定基因表达量
# import pandas as pd
# id_mapping_df = pd.read_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure5\rubber\rubber.xlsx", sheet_name="Sheet2")
# data_df = pd.read_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure5\rubber\Eu_tissue.xlsx")
# mapped_ids = id_mapping_df['ID']
# mapped_df = data_df[data_df['GeneID'].isin(mapped_ids)]
# mapped_df.to_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure5\rubber\Eu_tissue_mapped_gene.xlsx", index=False)

# # 替换ID
# import pandas as pd
# map_file = r"D:\Desktop\peptidemicro\00file\01figure\figure5\rubber\rubber.xlsx"
# map_df = pd.read_excel(map_file, sheet_name="Sheet2")
# map_df = map_df[["ID", "name"]]
# id_to_name = dict(zip(map_df["ID"], map_df["name"]))
# data_file = r"D:\Desktop\peptidemicro\00file\01figure\figure5\rubber\Eu_tissue_mapped_gene_pearson.xlsx"
# data_df = pd.read_excel(data_file, sheet_name="Sheet1")
# data_df["Var2"] = data_df["Var2"].map(id_to_name).fillna(data_df["Var2"])
# out_file = r"D:\Desktop\peptidemicro\00file\01figure\figure5\rubber\Var2_replaced_with_name.xlsx"
# data_df.to_excel(out_file, index=False)
# print("替换完成，结果已保存为：", out_file)

# # 提取|r|>0.8, p<0.05的基因的tpm表达量，然后计算相关性
# import pandas as pd
# id_mapping_df = pd.read_excel("/Volumes/caca/work_mechanism/figure5/rubber/Eu_tissue_mapped_gene_pearson.xlsx", sheet_name="Sheet2")
# data_df = pd.read_excel("/Volumes/caca/work_mechanism/figure5/rubber/Eu_tissue.xlsx")
# mapped_ids = id_mapping_df['ID']
# mapped_df = data_df[data_df['GeneID'].isin(mapped_ids)]
# mapped_df.to_excel("/Volumes/caca/work_mechanism/figure5/rubber/Eu_tissue_mapped_gene_pearson_tpm.xlsx", index=False)

# 提取不同Group的量
import pandas as pd
IN_FILE = "/Volumes/caca/work_mechanism/figure5/rubber/Eu_tissue_mapped_gene_pearson.xlsx"
OUT_FILE = "/Volumes/caca/work_mechanism/figure5/rubber/Eu_tissue_mapped_gene_pearson_group.xlsx"
ID_SHEET = "id_group"
DATA_SHEET = "Sheet1"
def unique_keep_order(seq):
    seen = set()
    out = []
    for x in seq:
        if pd.isna(x):
            continue
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out
def main():
    id_df = pd.read_excel(IN_FILE, sheet_name=ID_SHEET)
    id_df = id_df[["Var2", "Group"]].dropna(subset=["Var2", "Group"])
    var2_to_group = dict(zip(id_df["Var2"], id_df["Group"]))
    df = pd.read_excel(IN_FILE, sheet_name=DATA_SHEET)
    df = df[["Var1", "Var2"]].dropna(subset=["Var1", "Var2"])
    df["Group"] = df["Var2"].map(var2_to_group)
    unmapped = df[df["Group"].isna()].copy()
    df_mapped = df.dropna(subset=["Group"]).copy()
    grouped = (
        df_mapped.groupby("Group")["Var1"]
        .apply(lambda s: unique_keep_order(s.tolist()))
        .to_dict()
    )
    out_df = pd.DataFrame({g: pd.Series(v) for g, v in grouped.items()})
    with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
        out_df.to_excel(writer, sheet_name="Group_to_Var1", index=False)
        if not unmapped.empty:
            unmapped.to_excel(writer, sheet_name="Unmapped", index=False)
if __name__ == "__main__":
    main()

# # 提取候选肽基因的表达量
# import pandas as pd
# accession_file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\output_candidates_filtered.xlsx"
# count_file = r"D:\Desktop\peptidemicro\00file\01figure\figure5\total_all_matrix.xlsx"
# out_file = r"D:\Desktop\peptidemicro\00file\01figure\figure5\total_all_matrix_filtered.xlsx"
# accession_df = pd.read_excel(accession_file, sheet_name="task1_all_no_overlap", engine="openpyxl")
# count_df = pd.read_excel(count_file, engine="openpyxl")
# accessions = set(accession_df["accession"].dropna().astype(str))
# filtered_count_df = count_df[count_df["GeneID"].astype(str).isin(accessions)]
# filtered_count_df.to_excel(out_file, index=False)
