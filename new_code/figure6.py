# # 提取特定基因表达量
# import pandas as pd
# id_mapping_df = pd.read_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure6\deseq.xlsx", sheet_name="dt")
# data_df = pd.read_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure6\total_matrix_tpm.xlsx")
# mapped_ids = id_mapping_df['ID']
# mapped_df = data_df[data_df['GeneID'].isin(mapped_ids)]
# mapped_df.to_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure6\dt_deseq_tpm.xlsx", index=False)
# print(mapped_df)

# # 提取序列对应的基因组序列
# from Bio import SeqIO
# import pandas as pd
# sp_info_file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\output_best.xlsx"
# sp_info_df = pd.read_excel(sp_info_file, sheet_name="codon_1.0")
# id_file = r"D:\Desktop\peptidemicro\00file\01figure\figure6\引物设计.xlsx"
# id_df = pd.read_excel(id_file)
# id_list = id_df['ID_dup'].tolist()
# genome_fasta_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\Eu_genome.fasta"
# genome_dict = {}
# for record in SeqIO.parse(genome_fasta_file, "fasta"):
#     genome_dict[record.id] = record.seq
# mapped_id_df = sp_info_df[sp_info_df['ID'].isin(id_list)]
# for index, row in mapped_id_df.iterrows():
#     chrom = row['chrom']
#     start = row['start']
#     end = row['end']0
#     strand = row['strand']
#     total_score = row['total_score']
#     prior = row['prior']
#     gene_seq = genome_dict[chrom][start-1:end]
#     if strand == '-':
#         gene_seq = gene_seq.reverse_complement()
#         extern_seq = genome_dict[chrom][end:end+50].reverse_complement()
#     else:
#         extern_seq = genome_dict[chrom][start-51:start-1]
#     mapped_id_df.at[index, 'gene_seq'] = str(gene_seq)
#     mapped_id_df.at[index, 'extern_seq'] = str(extern_seq)
# mapped_id_df.to_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure6\引物设计序列.xlsx", index=False)
