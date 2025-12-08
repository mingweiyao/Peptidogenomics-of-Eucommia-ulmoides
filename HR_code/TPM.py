"""
这个脚本是用来读取GFF文件并进行TPM标准化的
"""
# import gffutils
# def prepare_length_data(gff_file):
#     if not os.path.exists(gff_file + '.db'):
#         print("🔄 正在创建GFF数据库...")
#         gffutils.create_db(
#             gff_file,
#             dbfn=gff_file + '.db',
#             force=True,
#             keep_order=True,
#             merge_strategy='merge',
#             id_spec={'gene': 'ID', 'mRNA': 'ID', 'CDS':'Parent'},
#             disable_infer_genes=True,
#             disable_infer_transcripts=True
#         )
#     db = gffutils.FeatureDB(gff_file + '.db')
#     gene_lengths = {}
#     for gene in db.features_of_type('gene'):
#         total_length = 0
#         for mRNA in db.children(gene, featuretype='mRNA'):
#             exons = list(db.children(mRNA, featuretype='exon'))
#             if exons:
#                 mRNA_length = sum(e.end - e.start + 1 for e in exons)
#                 total_length += mRNA_length
#         if total_length > 0:
#             gene_id = gene.id.replace('evm.model.', 'evm.TU.')
#             gene_lengths[gene_id] = total_length
#     gene_length_df = pd.DataFrame(list(gene_lengths.items()), columns=['GeneID', 'length'])   
#     return gene_length_df
# def normalize_tpm(count_df, length_df, output_file):
#     df = pd.merge(count_df, length_df, on='GeneID', how='inner')
#     df = df[df['length'] > 0]
#     sample_cols = [col for col in df.columns if col not in ['GeneID', 'length']]
#     tpm_data = {}
#     for sample in sample_cols:
#         rpk = (df[sample] * 10**3) / df['length']
#         per_million_scaling_factor = rpk.sum() / 10**6
#         tpm = rpk / per_million_scaling_factor
#         tpm_data[sample] = tpm
#     tpm_df = pd.concat([df[['GeneID']], pd.DataFrame(tpm_data)], axis=1)
#     tpm_df.to_csv(output_file, index=False)
#     print(f"TPM标准化完成: {output_file} (总条目数: {len(tpm_df)})")
#     return tpm_df