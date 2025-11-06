import pandas as pd
import os
from tqdm import tqdm
import numpy as np
from Bio import SeqIO
import re

def merge_count_files(input_dir, rna_info_file, gene_id_col="Geneid"):
    count_files = pd.read_excel(rna_info_file, sheet_name="Sheet2")
    merged_df = None
    for _, row in tqdm(count_files.iterrows(), desc="合并进度"):
        file = row['Sample']
        sample_name = f"{file}_counts.txt"
        file_path = os.path.join(input_dir, sample_name)
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#')
            counts = df[[gene_id_col, df.columns[-1]]]
            counts.columns = ['GeneID', file]
            if merged_df is None:
                merged_df = counts
            else:
                merged_df = pd.merge(merged_df, counts, on='GeneID', how='outer')
        except Exception as e:
            print(f"\n处理失败 {file}: {str(e)}")
            continue
    if merged_df is not None:
        print(f"\n合并后数据维度: {merged_df.shape}")
        return merged_df
    else:
        raise ValueError(f"错误：未成功合并任何数据")    

def filter_codon_gene(expression_matrix, cds_file):
    kept_id = []
    for rec in SeqIO.parse(cds_file, "fasta"):
        kept_id.append(re.sub(r'\.[^.]*$', '', rec.id))
    return expression_matrix[expression_matrix['GeneID'].isin(kept_id)]

def filter_expressed_genes(count_df, output_prefix, mean_threshold=5):
    expr_matrix = count_df.drop(columns=['GeneID'])
    nonzero_mean = expr_matrix.replace(0, np.nan).mean(axis=1, skipna=True)
    has_expression = (expr_matrix > 0).any(axis=1)
    condition = (nonzero_mean >= mean_threshold) & (has_expression)
    expressed_genes = count_df[condition]
    expressed_file = f"{output_prefix}_expressed.csv"
    expressed_genes.to_csv(expressed_file, index=False)
    print(f"\n基因筛选结果：")
    print(f"总基因数: {len(count_df)}")
    print(f"表达基因数（均值≥{mean_threshold}且至少1样本表达）: {len(expressed_genes)}")
    print(f"结果保存至: {expressed_file}")
    return expressed_genes

def prepare_length_data(cds_file):
    gene_lengths = {}
    for rec in SeqIO.parse(cds_file, "fasta"):
        _, _, start, end = rec.description.split('\t')
        rec_id = re.sub(r'\.[^.]*$', '', rec.id)
        gene_lengths[rec_id] = int(end) - int(start) + 1
    gene_length_df = pd.DataFrame(list(gene_lengths.items()), columns=['GeneID', 'length'])
    return gene_length_df

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
    tpm_df.to_csv(output_file, index=False)
    print(f"TPM标准化完成: {output_file} (总条目数: {len(tpm_df)})")

def extract_gene_locations(expression_df, cds_file):
    gene_info = {}
    for rec in SeqIO.parse(cds_file, "fasta"):
        chrom_temp, strand, start, end = rec.description.split('\t')
        chrom = chrom_temp.split(" ")[-1]
        rec_id = re.sub(r'\.[^.]*$', '', rec.id)
        if rec_id in expression_df['GeneID'].values:
            gene_info[rec_id] = (chrom, strand, start, end)
    return gene_info

def gene_location(expression_df, element_file, cds_file, output_file):
    gene_locations = extract_gene_locations(expression_df, cds_file)
    element_df = pd.read_excel(element_file, sheet_name="Genomic_Features")
    results = []
    for gene_id, gene_info in gene_locations.items():
        chrom, strand, start, end = gene_info
        chrom_elements = element_df[(element_df['chrom'] == chrom) & (element_df['strand'] == strand)].sort_values('start').reset_index(drop=True)
        overlapping_elements = []
        for _, elem in chrom_elements.iterrows():
            if not (elem['end'] < int(start) or elem['start'] > int(end)):
                overlapping_elements.append(elem['type'])
        unique_elements = []
        for elem_type in overlapping_elements:
            if elem_type not in unique_elements:
                unique_elements.append(elem_type)
        elements_string = "-".join(unique_elements) if unique_elements else None
        results.append({
            'gene_id': gene_id,
            'chrom': chrom,
            'strand': strand,
            'start': start,
            'end': end,
            'spanning_elements': elements_string,
            'element_count': len(unique_elements)
        })
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_file, index=False)

def main():
    input_dir = r"G:\Eu_peptido\20251018imeta\new_file_no_predict\01new_gene\quantitative"
    rna_info_file = r"G:\Eu_peptido\20251018imeta\new_file_no_predict\00raw\Total_rna_seq.xlsx"
    cds_file = r"G:\Eu_peptido\20251018imeta\new_file_no_predict\01new_gene\analysis_results\codon_filter_by_ms_gtf_cds.fasta"
    element_file = r"G:\Eu_peptido\20251018imeta\new_file_no_predict\00raw\peptide_analysis_results.xlsx"
    output_dir = r"G:\Eu_peptido\20251018imeta\new_file_no_predict\01new_gene"
    os.makedirs(output_dir, exist_ok=True)
    print("=== 步骤1/4: 合并计数文件 ===")
    expression_matrix = merge_count_files(input_dir, rna_info_file)
    # 筛选序列中没有终止密码子的gene的表达量
    expression_matrix_filter = filter_codon_gene(expression_matrix, cds_file)
    print("\n=== 步骤2/4: 筛选表达基因 ===")
    # 筛选序列中没有终止密码子的真实表达的gene
    expression_genes = filter_expressed_genes(expression_matrix_filter, os.path.join(output_dir, "total"))
    print("\n=== 步骤3/4: TPM标准化 ===")
    # 对序列中没有终止密码子的真实表达的gene的表达量TPM标准化
    length_df = prepare_length_data(cds_file)
    normalize_tpm(
        expression_genes,
        length_df,
        os.path.join(output_dir, "total_expressed_tpm.csv")
    )
    print("\n=== 步骤4/4: 基因位置分类 ===")
    # 对序列中没有终止密码子的真实表达的gene的位置进行分类
    gene_location(expression_genes, element_file, cds_file, os.path.join(output_dir, "total_expression_matrix.csv"))


if __name__ == "__main__":
    main()