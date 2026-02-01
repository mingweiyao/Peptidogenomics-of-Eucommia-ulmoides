"""
这个脚本是用来读取GFF文件并进行TPM标准化的
"""
# TPM标准化
import os
import pandas as pd
import gffutils
GFF_FILE = r"F:\CF\Tks_2021.gff"
COUNT_FILE = r"F:\CF\merged_counts.xlsx"
OUT_DIR = r"F:\CF\normalized_output"
def prepare_length_data(gff_file):
    db_path = gff_file + ".db"
    if not os.path.exists(db_path):
        print("🔄 正在创建 GFF 数据库...")
        gffutils.create_db(
            gff_file,
            dbfn=db_path,
            force=True,
            keep_order=True,
            merge_strategy="merge",
            id_spec={"gene": "ID", "mRNA": "ID", "CDS": "Parent"},
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
    db = gffutils.FeatureDB(db_path)
    gene_lengths = {}
    for gene in db.features_of_type("gene"):
        total_length = 0
        for mrna in db.children(gene, featuretype="mRNA", order_by="start"):
            exons = list(db.children(mrna, featuretype="exon", order_by="start"))
            if exons:
                mrna_len = sum(e.end - e.start + 1 for e in exons)
                total_length += mrna_len
        if total_length > 0:
            gene_id = gene.id.replace("evm.model.", "evm.TU.")
            gene_lengths[gene_id] = total_length
    length_df = pd.DataFrame(
        gene_lengths.items(),
        columns=["GeneID", "length"]
    )
    out_len = os.path.join(OUT_DIR, "gene_lengths.csv")
    length_df.to_csv(out_len, index=False)
    print(f"✅ 基因长度表已生成：{out_len}")
    return length_df
def read_counts(count_file):
    if count_file.lower().endswith(".csv"):
        df = pd.read_csv(count_file)
    else:
        df = pd.read_excel(count_file)
    if df.columns[0] != "GeneID":
        df = df.rename(columns={df.columns[0]: "GeneID"})
    return df
def normalize_tpm(count_df, length_df, output_file):
    df = pd.merge(count_df, length_df, on='GeneID', how='inner')
    df = df[df['length'] > 0]
    sample_cols = [col for col in df.columns if col not in ['GeneID', 'length']]
    tpm_data = {}
    for sample in sample_cols:
        rpk = (df[sample] * 10**3) / df['length']
        per_million_scaling_factor = rpk.sum() / 10**6
        tpm = rpk / per_million_scaling_factor
        tpm_data[sample] = tpm
    tpm_df = pd.concat([df[['GeneID']], pd.DataFrame(tpm_data)], axis=1)
    tpm_df.to_excel(output_file, index=False)
if __name__ == "__main__":
    length_df = prepare_length_data(GFF_FILE)
    count_df = read_counts(COUNT_FILE)
    normalize_tpm(count_df, length_df, os.path.join(OUT_DIR, "tpm_normalized.xlsx"))

# # FPKM标准化
# import os
# import pandas as pd
# import gffutils
# GFF_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/codon/codon_prediction/codon_prediction_v7/sp_codon.gff"
# COUNT_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/rubber_count_5.xlsx"
# OUT_DIR = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/"
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
#     length_df = pd.DataFrame(gene_lengths.items(), columns=["GeneID", "length"])
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
#     # 确保 GeneID 为字符串且去空格
#     df["GeneID"] = df["GeneID"].astype(str).str.strip()
#     return df
# def normalize_fpkm(count_df, length_df, output_file):
#     length_df = length_df.copy()
#     length_df["GeneID"] = length_df["GeneID"].astype(str).str.strip()
#     df = pd.merge(count_df, length_df, on="GeneID", how="inner")
#     df = df[df["length"] > 0].copy()
#     sample_cols = [c for c in df.columns if c not in ["GeneID", "length"]]
#     fpkm_data = {}
#     for sample in sample_cols:
#         counts = pd.to_numeric(df[sample], errors="coerce").fillna(0)
#         N = counts.sum()
#         if N == 0:
#             # 全是0时，直接输出0
#             fpkm = counts * 0.0
#         else:
#             L = df["length"].astype(float)  # bp
#             fpkm = (1e9 * counts) / (N * L)
#         fpkm_data[sample] = fpkm
#     fpkm_df = pd.concat([df[["GeneID"]], pd.DataFrame(fpkm_data)], axis=1)
#     fpkm_df.to_excel(output_file, index=False)
#     print(f"✅ FPKM 输出完成：{output_file}")
# if __name__ == "__main__":
#     os.makedirs(OUT_DIR, exist_ok=True)
#     length_df = prepare_length_data(GFF_FILE)
#     count_df = read_counts(COUNT_FILE)
#     normalize_fpkm(count_df, length_df, os.path.join(OUT_DIR, "rubber_count_5_fpkm.xlsx"))


# import pandas as pd
# import os
# def main():
#     base_dir = r"F:\CF"
#     all_files = os.listdir(base_dir)
#     gene_id_col = 'Geneid'
#     merge_df = None
#     for file in all_files:
#         if file.endswith(".txt"):
#             df = pd.read_csv(os.path.join(base_dir, file), sep="\t", comment='#')
#             counts = df[[gene_id_col, df.columns[-1]]]
#             counts.columns = ['Geneid', file]
#             if merge_df is None:
#                 merge_df = counts
#             else:
#                 merge_df = pd.merge(merge_df, counts, on='Geneid', how='outer')
#     merge_df.to_excel(os.path.join(base_dir, "merged_counts.xlsx"), index=False)
#     print("合并完成，输出文件: merged_counts.xlsx")
# if __name__ == "__main__":
#     main()