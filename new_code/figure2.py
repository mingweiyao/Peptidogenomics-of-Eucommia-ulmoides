import pandas as pd
import os

def data_analysis(info_file, output_dir):
    info_df = pd.read_excel(info_file, sheet_name="unique")
    # figure1: strand analysis
    group_counts = info_df.groupby(['type', 'strand']).size().unstack(fill_value=0)
    group_counts.to_excel(os.path.join(output_dir, "figure_a_strand.xlsx"), index=True)
    # figure2: length distribution
    type_length_df = info_df[['type', 'length']]
    type_length_pivot = type_length_df.pivot(columns='type', values='length')   
    type_length_pivot.to_excel(os.path.join(output_dir, "figure_b_length.xlsx"), index=False)
    # figure3: mw distribution
    type_mw_df = info_df[['type', 'molecular_weight']]
    type_mw_pivot = type_mw_df.pivot(columns='type', values='molecular_weight')
    type_mw_pivot.to_excel(os.path.join(output_dir, "figure_c_mw.xlsx"), index=False)
    # figure4: pI distribution
    type_pi_df = info_df[['type', 'isoelectric_point']]
    type_pi_pivot = type_pi_df.pivot(columns='type', values='isoelectric_point')
    type_pi_pivot.to_excel(os.path.join(output_dir, "figure_d_pi.xlsx"), index=False)
    # figure5: GRAVY distribution
    type_gravy_df = info_df[['type', 'gravy']]
    type_gravy_pivot = type_gravy_df.pivot(columns='type', values='gravy')
    type_gravy_pivot['NCP'] = 'NCP'
    cols = ['NCP'] + [col for col in type_gravy_pivot.columns if col != 'NCP']
    type_gravy_pivot = type_gravy_pivot[cols]
    type_gravy_pivot.columns = ['Type'] + list(type_gravy_pivot.columns[1:])
    type_gravy_pivot.to_excel(os.path.join(output_dir, "figure_e_gravy.xlsx"), index=False)
    # figure6: instability_index distribution
    type_instability_df = info_df[['type', 'instability_index']]
    type_instability_pivot = type_instability_df.pivot(columns='type', values='instability_index')
    type_instability_pivot['NCP'] = 'NCP'
    cols = ['NCP'] + [col for col in type_instability_pivot.columns if col != 'NCP']
    type_instability_pivot = type_instability_pivot[cols]
    type_instability_pivot.columns = ['Type'] + list(type_instability_pivot.columns[1:])
    type_instability_pivot.to_excel(os.path.join(output_dir, "figure_f_instability.xlsx"), index=False)

def main():
    info_file = r"D:\Desktop\peptidemicro\00file\01figure\finally_expressed_sp_info.xlsx"
    output_dir = r"D:\Desktop\peptidemicro\00file\01figure\figure2"
    data_analysis(info_file, output_dir)

if __name__ == "__main__":
    main()
