# # 转录组数据基础的NCP分布veen图
# import pandas as pd
# gene_expression_df = pd.read_csv(r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\03ouput\finally_filtered_genes.csv", index_col=0)
# tissue_mapping_df = pd.read_excel(r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\Total_rna_seq.xlsx", sheet_name="group")
# tissue_group_to_samples = (tissue_mapping_df.groupby(['Tissues', 'Group'])['Sample'].apply(list).to_dict())
# peptide_ids_by_tissue = {}
# for (tissue, group), samples in tissue_group_to_samples.items():
#     samples = [s for s in samples if s in gene_expression_df.columns]
#     group_expr = gene_expression_df.loc[:, samples]
#     group_mask = (group_expr >= 1).all(axis=1)
#     if tissue not in peptide_ids_by_tissue:
#         peptide_ids_by_tissue[tissue] = set()
#     peptide_ids_by_tissue[tissue].update(group_expr.index[group_mask].tolist())
# peptide_ids_by_tissue = {tissue: list(ids) for tissue, ids in peptide_ids_by_tissue.items()}
# max_length = max((len(ids) for ids in peptide_ids_by_tissue.values()), default=0)
# peptide_ids_df = pd.DataFrame({
#     tissue: ids + [None] * (max_length - len(ids))
#     for tissue, ids in peptide_ids_by_tissue.items()
# })
# output_excel_path = r"D:\Desktop\peptidemicro\00file\01figure\figure4\rnaseq_veen.xlsx"
# peptide_ids_df.to_excel(output_excel_path, index=False)

# # 数量-肽统计
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# from brokenaxes import brokenaxes
# import warnings
# warnings.filterwarnings('ignore')
# def plot_peptide_distribution_bar_broken_y(expression_matrix_path, output_pdf_path, output_csv_path, threshold=1, ylims=((0, 400), (600, 800), (2800, 3000)), bar_color="#4C72B0", edge_color="white"):
#     if expression_matrix_path.endswith('.csv'):
#         df = pd.read_csv(expression_matrix_path, index_col=0)
#     elif expression_matrix_path.endswith(('.xls', '.xlsx')):
#         df = pd.read_excel(expression_matrix_path, index_col=0)
#     else:
#         raise ValueError("仅支持 csv / xlsx")
#     print(f"读取矩阵：{df.shape[0]} 个ID × {df.shape[1]} 个列")
#     # ---------- 2. 二值化 & 统计 ----------
#     binary_df = (df >= threshold).astype(int)
#     counts_per_id = binary_df.sum(axis=1)
#     count_dist = counts_per_id.value_counts().sort_index()
#     total_ids = len(counts_per_id)
#     plot_data = pd.DataFrame({
#         "Number_of_Columns": count_dist.index.astype(int),
#         "Number_of_IDs": count_dist.values.astype(int),
#         "Percentage": (count_dist.values / total_ids * 100).round(2)
#     })
#     # ---------- 3) 输出 CSV ----------
#     plot_data.to_csv(output_csv_path, index=False)
#     print(f"✅ 统计结果已保存为：{output_csv_path}")
#     # ---------- 4) 画柱状图（broken y-axis） ----------
#     fig = plt.figure(figsize=(10, 8))
#     bax = brokenaxes(ylims=ylims, hspace=0.04, fig=fig)
#     x = plot_data["Number_of_Columns"].values
#     y = plot_data["Number_of_IDs"].values
#     bax.bar(x, y, width=0.8, color=bar_color, edgecolor=edge_color, linewidth=1.0)
#     bax.set_xlabel(f"Number of columns with CPM ≥ {threshold}")
#     bax.set_ylabel("Count of IDs")
#     bax.set_title("Peptide/ID Expression Distribution (Y-axis truncated)")
#     fig.savefig(output_pdf_path, bbox_inches="tight")
#     plt.close(fig)
#     print(f"✅ 截断柱状图已保存为：{output_pdf_path}")
#     return plot_data
# if __name__ == "__main__":
#     expression_file = r"D:\Desktop\peptidemicro\00file\01figure\finally_filtered_genes.csv"
#     output_pdf = r"D:\Desktop\peptidemicro\00file\01figure\figure4\peptide_distribution_bar_truncated.pdf"
#     output_csv = r"D:\Desktop\peptidemicro\00file\01figure\figure4\peptide_distribution_stat.csv"
#     plot_peptide_distribution_bar_broken_y(
#         expression_file,
#         output_pdf,
#         output_csv,
#         threshold=1,
#         ylims=((0, 400), (600, 800), (2800, 3000))
#     )

# 起始密码子预测
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
from multiprocessing import Pool
START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
STOP_CODONS  = {"TAA", "TAG", "TGA"}
MINUS_START_CODONS = {"CAT","CAG","CAC","CAA","CGT"}
MINUS_STOP_CODONS = {"TTA","CTA","TCA"}
CODON_PRIOR = {
    "ATG": 0.0,

    "CTG": -1.0,
    "GTG": -1.0,
    "TTG": -1.0,

    "ACG": -1.5,
    "AUA": -1.5,
    "AUU": -1.5,
    "AUC": -1.5,

    "AAG": -2.0,
    "AGG": -2.0,
    "CGU": -2.0,
    "CGC": -2.0,
    "CGG": -2.0,
    "CAG": -2.0
}
MINUS_CODON_PRIOR = {
    "CAT": 0.0,

    "CAG": -1.0,
    "CAC": -1.0,
    "CAA": -1.0,

    "CGT": -1.5,
    "TAT": -1.5,
    "AAT": -1.5,
    "GAT": -1.5,

    "CTT": -2.0,
    "CCT": -2.0,
    "ACG": -2.0,
    "GCG": -2.0,
    "CCG": -2.0,
    "CTG": -2.0
}
_genome_dict = None
_max_scan_nt = None
def init_worker(genome_file, max_scan_nt):
    global _genome_dict, _max_scan_nt
    _max_scan_nt = max_scan_nt
    _genome_dict = {}
    for rec in SeqIO.parse(genome_file, "fasta"):
        _genome_dict[rec.id] = rec.seq
def merge_segments(segments):
    min_start = min(segments, key=lambda x: x[0])[0]
    max_end = max(segments, key=lambda x: x[1])[1]  
    return min_start, max_end
def filter_peptide_seq_cal_cov(peptide_file, database_file):
    database_dict = {}
    for rec in SeqIO.parse(database_file, "fasta"):
        database_dict[rec.id] = rec.seq
    df_NCP = pd.read_excel(peptide_file, sheet_name="unique")
    df_NCP['proteins'] = df_NCP['proteins'].astype(int)
    df_NCP_filter = df_NCP[(df_NCP['proteins'] == 1)]
    accession_segments = {}
    for _, row in df_NCP_filter.iterrows():
        accession = row['accessions']
        start = row['start']
        end = row['end']
        chrom = row['chrom']
        strand = row['strand']
        if accession not in accession_segments:
            accession_segments[accession] = []
        accession_segments[accession].append((start, end, chrom, strand))
    stats = []
    for accession, segments in accession_segments.items():
        if accession in database_dict:
            seq = database_dict[accession]
            length = len(seq)
            min_start, max_end = merge_segments(segments)
            cov_length = (max_end-min_start+1) / 3
            coverage = cov_length / length
            chrom = segments[0][2]
            strand = segments[0][3]
            peptide_count = len(segments)
            stats.append({
                'accession': accession,
                'chrom': chrom,
                'strand': strand,
                'min_start': min_start,
                'max_end': max_end,
                'coverage': coverage,
                'peptide_count': peptide_count,
                'sequence_length_aa': length
            })
    return stats
def codon_test(seq, min_start, max_end, max_scan_nt):
    seq = str(Seq(str(seq)).upper())
    max_length = int(max_scan_nt / 3)
    for i in range(max_length):
        if seq[min_start-1-3*i:min_start+2-3*i] in START_CODONS:
            phy_start = min_start-3*i
            phy_value = CODON_PRIOR[seq[min_start-1-3*i:min_start+2-3*i]]
            break
    else:
        phy_start = None
        phy_value = None
    for i in range(max_length):
        if seq[max_end+i*3 : max_end+3*(i+1)] in STOP_CODONS:
            phy_end = max_end + i*3
            break
    else:
        phy_end = None
    return phy_start, phy_value, phy_end
def codon_test_minus(seq, min_start, max_end, max_scan_nt):
    seq = str(Seq(str(seq)).upper())
    max_length = int(max_scan_nt / 3)
    for i in range(max_length):
        if seq[min_start-1-3*(i+1):min_start-1-3*i] in MINUS_STOP_CODONS:
            phy_start = min_start-3*i
            break
    else:
        phy_start = None
    for i in range(max_length):
        if seq[max_end+3*(i-1):max_end+3*i] in MINUS_START_CODONS:
            phy_end = max_end+3*i
            phy_value = MINUS_CODON_PRIOR[seq[max_end+3*(i-1):max_end+3*i]]
            break
    else:
        phy_end = None
        phy_value = None
    return phy_start, phy_value, phy_end
def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank):
    genome_seq = Seq(str(genome_seq).upper())
    if strand == '+':
        ctx = genome_seq[phy_start-1-flank:phy_start+flank+2]
    else:
        ctx = genome_seq[phy_end-3-flank:phy_end+flank].reverse_complement()
    codon_pos = flank
    core = 0.0
    if ctx[codon_pos-3] in ('A', 'G'): core += 2.0
    if ctx[codon_pos+3] == 'G':       core += 2.0
    extended = 0.0
    for rel, ref in { -6:'G', -5:'C', -4:'C', -2:'C', +4:'G' }.items():
        if ctx[codon_pos+rel] == ref: extended += 0.5
    return {"core": core, "extended": extended, "total": core+extended, "context": ctx}
def run_scan_and_output_for_item(item):
    global _genome_dict, _max_scan_nt
    chrom = item['chrom']
    strand = item['strand']
    min_start = int(item['min_start'])
    max_end = int(item['max_end'])
    gseq = _genome_dict[chrom]

    if strand == '+':
        phy_start, phy_value, phy_end = codon_test(gseq, min_start, max_end, _max_scan_nt)
        prior_triplet = str(gseq[phy_start-1:phy_start+2]) if phy_start else None
    else:
        phy_start, phy_value, phy_end = codon_test_minus(gseq, min_start, max_end, _max_scan_nt)
        prior_triplet = str(gseq[phy_end-3:phy_end].reverse_complement()) if phy_end else None
    item['phy_start'] = phy_start
    item['phy_end'] = phy_end
    item['prior'] = prior_triplet
    item['total_score'] = None
    if phy_start and phy_end:
        kozak_results = kozak_score_cal(gseq, phy_start, phy_end, strand, flank=6)
        item['total_score'] = phy_value + kozak_results['total']
        item['kozak_seq'] = kozak_results['context']
    return item
def run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=100):
    with Pool(processes=nproc, initializer=init_worker, initargs=(genome_file, max_scan_nt)) as pool:
        stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
    return stats_update
def main():
    peptide_file = r"D:\Desktop\peptidemicro\00file\01figure\finally_expressed_sp_info.xlsx"
    database_file = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_peptide_database_customized_5.fa"
    genome_file = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_genome.fasta"
    output_file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\output.csv"    
    stats = filter_peptide_seq_cal_cov(peptide_file, database_file)
    stats_update = run_scan_and_output(stats, genome_file, max_scan_nt=300)    
    df = pd.DataFrame(stats_update)
    df.to_csv(output_file, index=False)
if __name__ == "__main__":
    main()