import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio.Seq import Seq
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
    accession_segments = {}
    for _, row in df_NCP.iterrows():
        accession = row['accessions']
        start = int(row['start'])
        end = int(row['end'])
        chrom = row['chrom']
        strand = row['strand']
        accession_segments.setdefault(accession, []).append((start, end, chrom, strand))
    stats = []
    for accession, segments in accession_segments.items():
        if accession not in database_dict:
            continue
        seq = database_dict[accession]
        length_aa = len(seq)
        min_start, max_end = merge_segments(segments)
        cov_length_aa = (max_end - min_start + 1) / 3.0
        coverage = cov_length_aa / length_aa
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
                'sequence_length_aa': length_aa
            })
    return stats

def find_best_orf_plus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
    seq_str = str(Seq(str(genome_seq)).upper())
    best = {
        "phy_start": None,
        "phy_end": None,
        "prior_triplet": None,
        "prior_value": None,
        "kozak_seq": None,
        "kozak_total": None,
        "total_score": None
    }
    for phy_start, triplet, prior_value in _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
        pass
def run_scan_and_output_for_item(item):
    global _genome_dict, _max_scan_nt
    chrom = item['chrom']
    strand = item['strand']
    min_start = int(item['min_start'])
    max_end = int(item['max_end'])
    if chrom not in _genome_dict:
        item['phy_start'] = None
        item['phy_end'] = None
        item['prior'] = None
        item['prior_value'] = None
        item['kozak_seq'] = None
        item['kozak_total'] = None
        item['total_score'] = None
        item['note'] = 'chrom_not_found'
        return item
    gseq = _genome_dict[chrom]
    if strand == '+':
        best = find_best_orf_plus(gseq, min_start, max_end, _max_scan_nt, flank=6)
    else:
        best = find_best_orf_minus(gseq, min_start, max_end, _max_scan_nt, flank=6)

def run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=20):
    with Pool(processes=nproc, initializer=init_worker, initargs=(genome_file, max_scan_nt)) as pool:
        stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
    return stats_update
def main():
    peptide_file  = r"D:\Desktop\peptidemicro\00file\01figure\finally_expressed_sp_info.xlsx"
    database_file = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_peptide_database_customized_5.fa"
    genome_file   = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_genome.fasta" # 是否更换
    output_file   = r"D:\Desktop\peptidemicro\00file\01figure\figure4\output.csv"
    stats = filter_peptide_seq_cal_cov(peptide_file, database_file)
    stats_update = run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=20)
if __name__ == "__main__":
    main()