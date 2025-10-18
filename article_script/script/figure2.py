import pandas as pd

def extract_data_by_id(id_file, data_file, output_file):
    ids = pd.read_csv(id_file, header=None)[0].tolist()
    data = pd.read_excel(data_file, sheet_name="NCP")
    filtered_data = data[data['ID'].isin(ids)]
    filtered_data.to_excel(output_file, index=False)
    return filtered_data

id_file = "/Volumes/caca/test_fractionation/00raw/rnaseq/03output/total_expressed.csv"
data_file = "/Volumes/caca/test_fractionation/00raw/sp_loc/Eu_sp_finally.xlsx"
output_file = "/Volumes/caca/test_fractionation/00raw/rnaseq/03output/sp_express_info_notdup.xlsx"

result = extract_data_by_id(id_file, data_file, output_file)

