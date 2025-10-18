import pandas as pd
import os
import gffutils
from tqdm import tqdm
import numpy as np

def merge_count_files(input_dir, RNA_info_file, output_file, gene_id_col="Geneid"):
    count_files = pd.read_excel(RNA_info_file, sheet_name="Sheet4")
    merged_df = None
    for _, row in tqdm(count_files.iterrows(), desc="合并进度"):
        file = row['ID']
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
        if merged_df.duplicated('GeneID').any():
            print(f"警告：存在重复基因ID，将取第一个出现的值")
            merged_df = merged_df.drop_duplicates('GeneID')
        merged_df.to_csv(output_file, index=False)
        return merged_df
    else:
        raise ValueError(f"错误：未成功合并任何数据")

def extract_and_filter_data(csv_df, excel_file_path, output_file_path):
    ids_to_match = csv_df["GeneID"].tolist()
    gff_make_df = pd.read_excel(excel_file_path, sheet_name='gff_make')
    print(f"从Excel中读取到 {len(gff_make_df)} 行gff_make数据")
    matched_gff = gff_make_df[gff_make_df['GeneID'].isin(ids_to_match)]
    print(f"在gff_make表中匹配到 {len(matched_gff)} 行数据")
    final_matches = matched_gff[matched_gff['type'] != "CDS(out of frame)"]
    ncps_df = pd.read_excel(excel_file_path, sheet_name='NCPs')
    match_columns = ['accessions', 'sequence', 'chrom', 'strand', 'type']
    match_info = pd.merge(final_matches[match_columns + ['GeneID']], ncps_df, on=match_columns, how='inner')
    print(f"在NCPs表中匹配到 {len(match_info)} 行数据")  
    ids_none_cds = match_info["GeneID"].tolist()
    csv_none_cds_df = csv_df[csv_df["GeneID"].isin(ids_none_cds)]
    csv_none_cds_df.to_csv(output_file_path, index=False) 
    return csv_none_cds_df

def combine_gene_sp_data(gene_df, sp_df, output_file):
    combined_df = pd.concat([gene_df, sp_df], axis=0)
    combined_df.to_csv(output_file, index=False)
    print(f"\n合并后总数据维度: {combined_df.shape}")
    return combined_df

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

def prepare_length_data(gff_file, length_excel):
    if not os.path.exists(gff_file + '.db'):
        print("🔄 正在创建GFF数据库...")
        gffutils.create_db(
            gff_file,
            dbfn=gff_file + '.db',
            force=True,
            keep_order=True,
            merge_strategy='merge',
            id_spec={'gene': 'ID', 'mRNA': 'ID', 'CDS':'Parent'},
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
    db = gffutils.FeatureDB(gff_file + '.db')
    gene_lengths = {}
    for gene in db.features_of_type('gene'):
        total_length = 0
        for mRNA in db.children(gene, featuretype='mRNA'):
            exons = list(db.children(mRNA, featuretype='exon'))
            if exons:
                mRNA_length = sum(e.end - e.start + 1 for e in exons)
                total_length += mRNA_length
        if total_length > 0:
            gene_id = gene.id.replace('evm.model.', 'evm.TU.')
            gene_lengths[gene_id] = total_length
    gene_length_df = pd.DataFrame(list(gene_lengths.items()), columns=['GeneID', 'length'])   
    sp_length_df = pd.read_excel(length_excel, sheet_name='gff_make', usecols=['GeneID', 'length'])
    sp_length_df['length'] = sp_length_df['length'] * 3
    sp_length_df.rename(columns={'GeneID': 'GeneID'}, inplace=True)
    length_df = pd.concat([gene_length_df, sp_length_df], axis=0)
    return length_df

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
    return tpm_df

def extract_condition_samples(count_df, tpm_df, RNA_info_file, output_prefix, sheet_name="Sheet5"):
    condition_data = pd.read_excel(RNA_info_file, sheet_name=sheet_name)
    for condition in condition_data.columns:
        sample_ids = condition_data[condition].dropna().tolist()
        available_samples = [col for col in count_df.columns if col in sample_ids]
        condition_count_df = count_df[['GeneID'] + available_samples]
        count_output = f"{output_prefix}_{condition}_counts.csv"
        condition_count_df.to_csv(count_output, index=False)
        condition_tpm_df = tpm_df[['GeneID'] + available_samples]
        tpm_output = f"{output_prefix}_{condition}_tpm.csv"
        condition_tpm_df.to_csv(tpm_output, index=False)
        log2_tpm_df = condition_tpm_df.copy()
        log2_tpm_df[available_samples] = np.log2(log2_tpm_df[available_samples] + 0.1)
        log2_output = f"{output_prefix}_{condition}_log2tpm.csv"
        log2_tpm_df.to_csv(log2_output, index=False)
        print(f"\n条件 [{condition}] 提取结果:")
        print(f"  样本数: {len(available_samples)}/{len(sample_ids)}")
        print(f"  Count文件: {count_output}")
        print(f"  TPM文件: {tpm_output}")
        if len(available_samples) < len(sample_ids):
            missing = set(sample_ids) - set(available_samples)
            print(f"  缺失样本: {missing}")

def extract_data(csv_df, excel_file_path, match_ncps_file):
    ids_to_match = csv_df["GeneID"].tolist()
    gff_make_df = pd.read_excel(excel_file_path, sheet_name='gff_make')
    print(f"从Excel中读取到 {len(gff_make_df)} 行gff_make数据")
    matched_gff = gff_make_df[gff_make_df['GeneID'].isin(ids_to_match)]
    print(f"在gff_make表中匹配到 {len(matched_gff)} 行数据")
    ncps_df = pd.read_excel(excel_file_path, sheet_name='NCPs')
    match_columns = ['accessions', 'sequence', 'chrom', 'strand', 'type']
    match_info = pd.merge(
        matched_gff[match_columns + ['GeneID']],
        ncps_df,
        on=match_columns,
        how='inner'
    )
    print(f"在NCPs表中匹配到 {len(matched_gff)} 行数据")  
    match_info.to_excel(match_ncps_file, index=False)

if __name__ == "__main__":
    base_dir = "/media/wanglab/caca/peptidegenomics"
    gene_input_dir = os.path.join(base_dir, "00raw_data/STAR_bam")
    sp_input_dir = os.path.join(base_dir, "00raw_data/STAR_bam_sp")
    RNA_info_file = os.path.join(base_dir, "01analysis_data/rnaseq/rna-seq.xlsx")
    gff_file = os.path.join(base_dir, "00raw_data/GWHBISF00000000.gff")
    length_excel = os.path.join(base_dir, "01analysis_data/sp_loc/Eu_sp_finally.xlsx")
    output_dir = os.path.join(base_dir, "01analysis_data/rnaseq")
    os.makedirs(output_dir, exist_ok=True)

    print("=== 步骤1/6: 合并计数文件 ===")
    gene_matrix = merge_count_files(
        gene_input_dir, RNA_info_file,
        output_file=os.path.join(output_dir, "gene_matrix.csv")
    )
    sp_matrix = merge_count_files(
        sp_input_dir, RNA_info_file,
        output_file=os.path.join(output_dir, "sp_matrix_all.csv")
    )

    print("=== 步骤2/6: 剔除CDS(out of frame)肽 ===")
    sp_matrix_use = extract_and_filter_data(
        sp_matrix, length_excel, 
        output_file_path=os.path.join(output_dir, "sp_matrix.csv")
    )

    print("\n=== 步骤3/6: 合并并筛选表达基因/小肽 ===")
    combined_matrix = combine_gene_sp_data(
        gene_matrix, sp_matrix_use,
        os.path.join(output_dir, "combined_matrix.csv")
    )
    expressed_genes_peptides = filter_expressed_genes(
        combined_matrix, 
        os.path.join(output_dir, "combined")
    )

    print("\n=== 步骤4/6: TPM标准化 ===")
    length_df = prepare_length_data(gff_file, length_excel)
    tpm_matrix = normalize_tpm(
        expressed_genes_peptides,
        length_df,
        os.path.join(output_dir, "combined_expressed_tpm.csv")
    )

    print("\n=== 步骤5/6: 提取条件样本数据 ===")
    extract_condition_samples(
        expressed_genes_peptides,
        tpm_matrix,
        RNA_info_file,
        os.path.join(output_dir, "condition")
    )
    print("=== 步骤6/6: 提取匹配肽信息===")
    sp_matrix_use = extract_data(
        expressed_genes_peptides, length_excel, 
        match_ncps_file=os.path.join(output_dir, "sp_matrix_info.xlsx")
    )

    print("\n所有处理完成！最终结果保存在:", output_dir)
