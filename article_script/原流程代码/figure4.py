import pandas as pd
import os
df = pd.read_excel("G:/peptidegenomics/01analysis_data/sp_loc/Eu_sp_finally.xlsx", sheet_name="NCPs")
result1 = df.pivot_table(
    index="strand", 
    columns="material", 
    values="accessions", 
    aggfunc="count", 
    fill_value=0
)
result1.to_csv("G:/peptidegenomics/02figure/figure4a/data1_strand_material_counts.txt", sep="\t")
print("数据一已保存 → G:/peptidegenomics/02figure/figure4a/data1_strand_material_counts.txt")

result2 = df.pivot_table(
    index="type", 
    columns="material", 
    values="accessions", 
    aggfunc="count", 
    fill_value=0
)
result2.to_csv("G:/peptidegenomics/02figure/figure4a/data2_type_material_counts.txt", sep="\t")
print("数据二已保存 → G:/peptidegenomics/02figure/figure4a/data2_type_material_counts.txt")

for material in df["material"].unique():
    material_data = df[df["material"] == material]
    length_data = material_data[["type", "length"]].sort_values(["type", "length"], ascending=[True, False])
    length_data.to_csv(f"G:/peptidegenomics/02figure/figure4a/data3_{material}_lengths.txt", sep="\t", index=False, header=False)
    print(f"数据三已保存 → G:/peptidegenomics/02figure/figure4a/data3_{material}_lengths.txt")
    mw_pivot = material_data.pivot_table(
        index=material_data.groupby("type").cumcount(),
        columns="type", 
        values="molecular_weight", 
        aggfunc="first"
    )
    mw_pivot.to_csv(f"G:/peptidegenomics/02figure/figure4a/MW_{material}.txt", sep="\t", index=False)
    print(f"数据四已保存 → G:/peptidegenomics/02figure/figure4a/MW_{material}.txt")
    pi_pivot = material_data.pivot_table(
        index=material_data.groupby("type").cumcount(), 
        columns="type", 
        values="isoelectric_point", 
        aggfunc="first"
    )
    pi_pivot.to_csv(f"G:/peptidegenomics/02figure/figure4a/PI_{material}.txt", sep="\t", index=False)
    print(f"数据五已保存 → G:/peptidegenomics/02figure/figure4a/PI_{material}.txt")

print("\n所有分析完成！请检查 G:/peptidegenomics/02figure/figure4a/ 目录下的文件。")