"""Scan for non-canonical start codons around MS-supported peptide spans (per accession).

For every six-frame database accession supported by unique peptides (type != exon_diff,
proteins == 1), compute the supported genomic span and coverage, then scan in-frame
(strand-aware) up to --max-scan-nt for candidate start codons (ATG/CTG/GTG/TTG/ACG).
Each candidate is scored with a codon prior plus a Kozak (core + extended) score;
candidate positions, triplets, priors and scores are written as semicolon-joined columns.

Inputs:  --peptide-file (aa-pipeline results xlsx, sheet "NCP"), --database-file
         (six-frame peptide database FASTA providing accession sequences), --genome-file
Output:  --output-file CSV, one row per accession
Notes:   runs a multiprocessing pool (--processes, default 100); prior table
         ATG 0.0 / CTG -1.0 / GTG -1.4 / TTG -1.5 / ACG -1.1.
Usage:
  py 04_scan_start_codons.py --peptide-file ./data/raw/peptide_analysis_results.xlsx \
      --database-file ./data/reference/Eu_peptide_database_customized_5.fa \
      --genome-file ./data/reference/Eu_genome.fasta \
      --output-file ./output/11_transcript_model_construction/start_codon_scan.csv
"""
import argparse
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
    "GTG": -1.4,
    "TTG": -1.5,
    "ACG": -1.1
}
MINUS_CODON_PRIOR = {
    "CAT": 0.0,
    "CAG": -1.0,
    "CAC": -1.4,
    "CAA": -1.5,
    "CGT": -1.1
}

def merge_segments(segments):
    min_start = min(segments, key=lambda x: x[0])[0]
    max_end = max(segments, key=lambda x: x[1])[1]  
    return min_start, max_end

def filter_peptide_seq_cal_cov(peptide_file, database_file):
    database_dict = {}
    for rec in SeqIO.parse(database_file, "fasta"):
        database_dict[rec.id] = rec.seq
    df_NCP = pd.read_excel(peptide_file, sheet_name="NCP")
    df_NCP['proteins'] = df_NCP['proteins'].astype(int)
    df_NCP_filter = df_NCP[(df_NCP['type'] != 'exon_diff') & (df_NCP['proteins'] == 1)]
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
    phy_end = None
    for i in range(max_length):
        start = max_end + i * 3
        end = max_end + (i + 1) * 3
        codon = seq[start:end]
        if codon in STOP_CODONS:
            phy_end = start
            break
    candidates = []
    for i in range(max_length):
        pos = min_start - 3 * i
        codon = seq[pos - 1:pos + 2]
        if codon in STOP_CODONS:
            break
        if codon in START_CODONS:
            candidates.append({
                "start":pos, "triplet":codon, "prior": CODON_PRIOR[codon]
            })
    if candidates:
        main_start = candidates[0]["start"]
        main_prior = candidates[0]["prior"]
    else:
        main_start = None
        main_prior = None
    return main_start, main_prior, phy_end, candidates

def codon_test_minus(seq, min_start, max_end, max_scan_nt):
    seq = str(Seq(str(seq)).upper())
    max_length = int(max_scan_nt / 3)
    phy_start = None
    for i in range(max_length):
        left = min_start-1-3*(i+1)
        right = min_start-1-3*i
        codon = seq[left:right]
        if codon in MINUS_STOP_CODONS:
            phy_start = min_start - 3 * i
            break
    candidates = []
    phy_end_main = None
    main_prior = None
    for i in range(max_length):
        start = max_end + 3 * (i - 1)
        end = max_end + 3 * i
        codon = seq[start:end]
        if codon in MINUS_STOP_CODONS and i != 0:
            break
        if codon in MINUS_START_CODONS:
            coord = max_end + 3 * i
            prior = MINUS_CODON_PRIOR[codon]
            candidates.append({
                "end": coord,
                "triplet": codon,
                "prior": prior,
            })
            if phy_end_main is None:
                phy_end_main = coord
                main_prior = prior
    return phy_start, main_prior, phy_end_main, candidates        

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

def run_scan_and_output_for_item(item, genome_dict, max_scan_nt):
    chrom = item['chrom']
    strand = item['strand']
    min_start = int(item['min_start'])
    max_end = int(item['max_end'])
    item['phy_start'] = None
    item['phy_end'] = None
    item['prior'] = None
    item['total_score'] = None
    item['all_codon_pos'] = None
    item['all_codon_triplet'] = None
    item['all_codon_prior'] = None
    item['all_codon_kozak_total'] = None
    item['all_codon_total_score'] = None

    gseq = genome_dict[chrom]
    if strand == '+':
        phy_start, phy_value, phy_end, candidates = codon_test(gseq, min_start, max_end, max_scan_nt)
        prior_triplet = str(gseq[phy_start-1:phy_start+2]) if phy_start else None
    else:
        phy_start, phy_value, phy_end, candidates = codon_test_minus(gseq, min_start, max_end, max_scan_nt)
        prior_triplet = str(gseq[phy_end-3:phy_end].reverse_complement()) if phy_end else None
    item['phy_start'] = phy_start
    item['phy_end'] = phy_end
    item['prior'] = prior_triplet
    if candidates:
        all_pos = []
        all_triplet = []
        all_prior = []
        all_kozak_total = []
        all_overall = []
        for cand in candidates:
            if strand == '+':
                pos = cand['start'] 
                trip = cand['triplet']
                prior = cand['prior']
                k_res = kozak_score_cal(gseq, pos, pos, strand, flank=6)
            else:
                pos = cand['end']
                trip = cand['triplet']
                prior = cand['prior']
                k_res = kozak_score_cal(gseq, 0, pos, strand, flank=6)
            k_total = k_res['total']
            overall = prior + k_total
            all_pos.append(pos)
            all_triplet.append(trip)
            all_prior.append(prior)
            all_kozak_total.append(k_total)
            all_overall.append(overall)
        item['total_score'] = all_overall[0]
        item['all_codon_pos'] = ';'.join(str(x) for x in all_pos)
        item['all_codon_triplet'] = ';'.join(all_triplet)
        item['all_codon_prior'] = ';'.join(str(x) for x in all_prior)
        item['all_codon_kozak_total'] = ';'.join(str(x) for x in all_kozak_total)
        item['all_codon_total_score'] = ';'.join(str(x) for x in all_overall)
    return item

def run_scan_and_output(stats, genome_file, max_scan_nt, processes=100):
    genome_dict = {}
    for rec in SeqIO.parse(genome_file, "fasta"):
        genome_dict[rec.id] = rec.seq
    with Pool(processes=processes) as pool:
        stats_update = pool.starmap(run_scan_and_output_for_item, [(item, genome_dict, max_scan_nt) for item in stats])
    return stats_update

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--peptide-file", required=True, help="aa-pipeline results xlsx (sheet 'NCP')")
    ap.add_argument("--database-file", required=True, help="six-frame peptide database FASTA")
    ap.add_argument("--genome-file", required=True, help="genome FASTA")
    ap.add_argument("--output-file", required=True, help="output CSV path")
    ap.add_argument("--processes", type=int, default=100, help="multiprocessing pool size")
    ap.add_argument("--max-scan-nt", type=int, default=300, help="max nt to scan for start codons")
    args = ap.parse_args()

    stats = filter_peptide_seq_cal_cov(args.peptide_file, args.database_file)
    stats_update = run_scan_and_output(stats, args.genome_file,
                                       max_scan_nt=args.max_scan_nt, processes=args.processes)

    df = pd.DataFrame(stats_update)
    df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
