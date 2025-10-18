# import pandas as pd
# from Bio.SeqUtils.ProtParam import ProteinAnalysis
# from tqdm import tqdm
# import logging
# logging.basicConfig(level=logging.WARNING)
# logger = logging.getLogger()
# def analyze_peptide_properties(sequence):
#     if not sequence or not isinstance(sequence, str):
#         return {
#             "length": 0,
#             "molecular_weight": None,
#             "isoelectric_point": None,
#             "gravy": None,
#             "aromaticity": None,
#             "instability_index": None
#         }
#     try:
#         analysis = ProteinAnalysis(sequence)
#         properties = {
#             "length": len(sequence),
#             "molecular_weight": analysis.molecular_weight(),
#             "isoelectric_point": analysis.isoelectric_point(),
#             "gravy": analysis.gravy(),
#         }
#         try:
#             properties.update({
#                 "instability_index": analysis.instability_index()
#             })
#         except Exception as e:
#             logger.warning(f"部分性质计算失败: {str(e)}")
#         return properties  
#     except Exception as e:
#         logger.error(f"肽段分析失败: {str(e)}")
#         return {
#             "length": len(sequence),
#             "molecular_weight": None,
#             "isoelectric_point": None,
#             "gravy": None,
#             "aromaticity": None,
#             "instability_index": None
#         }
# def batch_analyze_peptides_from_excel(file_path):
#     df = pd.read_excel(file_path, sheet_name="NCP")
#     if 'sequence' not in df.columns:
#         logger.error("Excel 文件中未找到 'sequence' 列")
#         return None
#     results = []
#     for _, row in tqdm(df.iterrows(), total=len(df), desc="分析肽段性质"):
#         seq = row.get('sequence', '')
#         results.append(analyze_peptide_properties(seq))
#     result_df = pd.DataFrame(results)
#     return result_df
# file_path = r'G:\test_fractionation\00raw\output\Eu_sp_finally.xlsx'
# output_df = batch_analyze_peptides_from_excel(file_path)
# if output_df is not None:
#     output_df.to_excel(r'G:\test_fractionation\00raw\output\output_with_peptide_properties.xlsx', index=False)
#     print("分析结果已保存到 'output_with_peptide_properties.xlsx'")

# 去除重复的数据行
import pandas as pd
def remove_duplicates(input_file, output_file, columns):
    df = pd.read_excel(input_file, sheet_name="NCP")
    df_unique = df.drop_duplicates(subset=columns, keep='first')
    df_unique.to_excel(output_file, index=False)
    print(f"去重完成! 原始行数: {len(df)}, 去重后行数: {len(df_unique)}")
if __name__ == "__main__":
    input_file = r"G:\test_fractionation\00raw\output\Eu_sp_finally.xlsx"
    output_file = r"G:\test_fractionation\00raw\output\output.xlsx"
    columns = ["sequence", "start", "end", "strand", "chrom"]
    remove_duplicates(input_file, output_file, columns)