"""
这个脚本包括了去除质谱数据中不同组织检测到的重复的肽段的代码
"""
# # 提取表达的NCP的信息
# import pandas as pd
# csv_file = r"F:\Eu_peptido\new_prepare\sp_expression.csv"
# excel_file = r"F:\Eu_peptido\00file\00raw\sp_loc\Eu_sp_finally.xlsx"
# out_file = r"F:\Eu_peptido\new_prepare\sp_expression_info.xlsx"
# csv_df = pd.read_csv(csv_file)
# id_list = csv_df.iloc[:, 0].dropna().astype(str).unique()
# excel_df = pd.read_excel(excel_file, sheet_name="NCP")
# excel_id_col = "ID"
# filtered_df = excel_df[excel_df[excel_id_col].astype(str).isin(id_list)]
# filtered_df.to_excel(out_file, index=False)
# print(f"完成！共筛到 {len(filtered_df)} 行，已保存到 {out_file}")

# import pandas as pd
# def extract_data_by_id(id_file, data_file, output_file):
#     ids = pd.read_csv(id_file, header=None)[0].tolist()
#     data = pd.read_excel(data_file, sheet_name="NCP")
#     filtered_data = data[data['ID'].isin(ids)]
#     filtered_data.to_excel(output_file, index=False)
#     return filtered_data
# id_file = "F:\Eu_peptido\new_preparetotal_expressed.csv"
# data_file = "/Volumes/caca/test_fractionation/00raw/sp_loc/Eu_sp_finally.xlsx"
# output_file = "/Volumes/caca/test_fractionation/00raw/rnaseq/03output/sp_express_info_notdup.xlsx"
# result = extract_data_by_id(id_file, data_file, output_file)

