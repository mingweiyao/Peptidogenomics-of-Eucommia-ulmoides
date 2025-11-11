import pandas as pd
from Bio import SeqIO

def merge_segments(segments):
    min_start = min(segments, key=lambda x: x[0])[0]
    max_end = max(segments, key=lambda x: x[1])[1]  
    return (max_end - min_start + 1) / 3

def analysis_peptide(peptide_file, database_file, output_file):
    database_dict = {}
    for rec in SeqIO.parse(database_file, "fasta"):
        database_dict[rec.id] = rec.seq
    df_NCP = pd.read_excel(peptide_file, sheet_name="NCP")
    df_NCP_filter = df_NCP[df_NCP['type'] != 'exon_diff']

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
            len_merge_seq = merge_segments(segments)
            coverage = len_merge_seq / length
            chrom = segments[0][2]
            strand = segments[0][3]
            peptide_count = len(segments)
            stats.append({
                'accession': accession,
                'chrom': chrom,
                'strand': strand,
                'coverage': coverage,
                'peptide_count': peptide_count,
                'sequence_length': length
            })
    stats_df = pd.DataFrame(stats)
    stats_df.to_excel(output_file, index=False)
    print(f"分析结果已保存至: {output_file}")

def main():
    peptide_file = r"D:\Desktop\peptidemicro\00raw\sp_loc\Eu_sp_finally.xlsx"
    database_file = r"D:\Desktop\peptidemicro\00raw\Eu_peptide_database_customized_5.fa"
    output_file = r"D:\Desktop\output_file.xlsx"
    analysis_peptide(peptide_file, database_file, output_file)

if __name__ == "__main__":
    main()
