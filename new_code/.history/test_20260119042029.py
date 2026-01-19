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
