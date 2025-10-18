# # figureA: 上调下调基因
# import pandas as pd
# def extract_data_by_id(id_file, data_file, output_file):
#     ids = pd.read_excel(id_file, sheet_name="drought")
#     data = pd.read_excel(data_file, sheet_name="DT_down")
#     filtered_data = data[data['name'].isin(ids['ID'])]
#     filtered_data.to_excel(output_file, index=False)
# id_file = "/Volumes/caca/test_fractionation/01figure/figure5/response_veen_analysis.xlsx"
# data_file = "/Volumes/caca/test_fractionation/01figure/figure5/DEGs_ST_DT/st_dt_DEGs.xlsx"
# output_file = "/Volumes/caca/test_fractionation/01figure/figure5/dt_down.xlsx"
# extract_data_by_id(id_file, data_file, output_file)

# import pandas as pd
# def extract_data_by_id(id_file, data_file, output_file):
#     ids = pd.read_excel(id_file, sheet_name="Sheet1")
#     data = pd.read_excel(data_file, sheet_name="DT")
#     filtered_data = data[data['name'].isin(ids['drought'])]
#     filtered_data.to_excel(output_file, index=False)
# id_file = r"G:\test_fractionation\01figure\figure5\response_veen.xlsx"
# data_file = r"G:\test_fractionation\01figure\figure5\DEGs_ST_DT\st_dt_DEGs.xlsx"
# output_file = r"G:\test_fractionation\01figure\figure5\drought_DEGs.xlsx"
# extract_data_by_id(id_file, data_file, output_file)

# import pandas as pd
# def count_values_in_columns(input_csv, output_excel, columns_to_analyze):
#     df = pd.read_csv(input_csv)
#     with pd.ExcelWriter(output_excel) as writer:
#         for column in columns_to_analyze:
#             if column in df.columns:
#                 value_counts = df[column].value_counts().reset_index()
#                 value_counts.columns = [column, 'Count']
#                 value_counts.to_excel(writer, sheet_name=column, index=False)
#             else:
#                 print(f"警告: 列 {column} 不存在于数据中")
# input_csv = "/Volumes/caca/test_fractionation/01figure/GO/dt_st/novel_genes_function_predictions.csv"
# output_excel = "/Volumes/caca/test_fractionation/01figure/figure5/GO.xlsx" 
# columns_to_analyze = ["BP_Description", "MF_Description", "CC_Description"]
# count_values_in_columns(input_csv, output_excel, columns_to_analyze)

import pandas as pd
import os

def extract_and_summarize(excel_file, target_path, output_file):
    df = pd.read_excel(excel_file, sheet_name="Sheet1")
    folder_names = df["ID"].tolist()
    all_data = pd.DataFrame()
    for folder in folder_names:
        folder_path = os.path.join(target_path, folder)
        if os.path.isdir(folder_path):
            csv_file = os.path.join(folder_path, 'full_go_results.csv')
            if os.path.exists(csv_file):
                df_go = pd.read_csv(csv_file)
                df_go['NCP_ID'] = folder
                all_data = pd.concat([df_go, all_data], ignore_index=True)
            else:
                print(f"文件 '{csv_file}' 不存在.")
        else:
            print(f"文件夹 '{folder_path}' 不存在.")
    all_data.to_csv(output_file, index=False)
    print(f"数据已汇总到 {output_file}")

excel_file = r'G:\test_fractionation\01figure\figure6\yield_CGA.xlsx'
target_path = r'G:\test_fractionation\01figure\GO\yield\gene_results'
output_file = r'G:\test_fractionation\01figure\figure6\yield_CGA_GO.csv'

extract_and_summarize(excel_file, target_path, output_file)
