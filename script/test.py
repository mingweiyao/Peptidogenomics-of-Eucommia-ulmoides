# import os
# from Bio import SeqIO

# def test_transcript(directory):
#     for f in os.listdir(directory):
#         if f.endswith('.gff3'):
#             with open(os.path.join(directory, f), 'r') as file:
#                 for line in file:
#                     if 'transcript' in line:
#                         print(f"{f}")
#                         break
#         elif f.endswith('.fasta'):
#             for record in SeqIO.parse(os.path.join(directory, f), "fasta"):
#                 if 'transcript' in record.description:
#                     print(f"{f}")
#                     break
# def main():
#     gff3_dictorty = r"G:\Eu_peptido\20251018 imeta\file\00raw\GFF3"
#     fasta_directory = r"G:\Eu_peptido\20251018 imeta\file\00raw\fasta"
#     test_transcript(gff3_dictorty)
#     test_transcript(fasta_directory)

# if __name__ == "__main__":
#     main()

# from Bio import SeqIO
# count = 0
# for rec in SeqIO.parse(r"G:\Eu_peptido\20251018 imeta\file\00raw\output\file2_nonredundant_coding_transcript_pep.fasta", "fasta"):
#     count += 1
# print(count)

# count = 0
# with open(r"G:\Eu_peptido\20251018 imeta\file\00raw\output\file1_nonredundant_coding_transcripts.gff3") as f:
#     for line in f:
#         if "gene" in line:
#             count+=1
# print(count)

import pandas as pd
import re
from datetime import datetime
from tqdm import tqdm
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor

RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
RE_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_CANDIDATE_INDEX = None

def _init_worker(candidate_index):
    global _CANDIDATE_INDEX
    _CANDIDATE_INDEX = candidate_index

def _worker_match_one_peptide(peptide_tuple):
    chrom, pstart, pend, strand = peptide_tuple
    hit_ids = []
    for tid, exon_regions in _CANDIDATE_INDEX.get((chrom, strand), []):
        tr_start = exon_regions[0][0]
        tr_end = exon_regions[-1][1]
        if not (tr_start <= pstart <= pend <= tr_end):
            continue
        matched = False
        for (exon_start, exon_end) in exon_regions:
            if pstart >= exon_start and pend <= exon_end:
                if (pstart - exon_start) % 3 == 0:
                    matched = True
                break
        if not matched and len(exon_regions) >= 2:
            for idx in range(len(exon_regions) - 1):
                e1_start, e1_end = exon_regions[idx]
                e2_start, e2_end = exon_regions[idx + 1]
                if (e1_start <= pstart <= e1_end) and (e2_start <= pend <= e2_end):
                    if (pstart - e1_start) % 3 == 0:
                        matched = True
                    break
        if matched:
            hit_ids.append(tid)             
    return hit_ids

def filter_by_ms(ms_file, gtf_file, workers=None):
    # 读取 GTF 文件
    transcript_lines = defaultdict(list)
    transcript_info = {}    
    current_gene_lines = []
    current_transcript_id = None

    with open(gtf_file, 'r') as f:
        for line in tqdm(f, desc="Loading GTF"):
            if not line or line.startswith('#'):
                continue
            line = line.rstrip('\n')
            chrom, _, ftype, start, end, _, strand, _, attrs = line.split('\t')
            if ftype == 'transcript':
                if current_transcript_id and current_gene_lines:
                    transcript_lines[current_transcript_id] = current_gene_lines
                tid_match = RE_TRANSCRIPT_ID.search(attrs)
                gid_match = RE_GENE_ID.search(attrs)
                if not (tid_match and gid_match):
                    continue
                tid = tid_match.group(1)
                gid = gid_match.group(1)
                current_transcript_id = tid
                current_gene_lines = [line]
                transcript_info[tid] = {
                    'gene_id': gid, 'chrom': chrom, 'strand': strand,
                    'start': int(start), 'end': int(end),
                }
            elif ftype == 'exon':
                if current_transcript_id:
                    current_gene_lines.append(line)
        if current_transcript_id and current_gene_lines:
            transcript_lines[current_transcript_id] = current_gene_lines

    transcript_exon_region = {}
    for tid, lines in tqdm(transcript_lines.items(), desc="Parsing exon regions"):
        exon = []
        for line in lines:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if f[2] == 'exon':
                exon.append((int(f[3]), int(f[4])))
        transcript_exon_region[tid] = sorted(exon)

    candidate_index = defaultdict(list)
    for tid, info in transcript_info.items():
        exons_region = transcript_exon_region[tid]
        if not exons_region:
            continue
        key = (info['chrom'], info['strand'])
        candidate_index[key].append((tid, exons_region))

    df = pd.read_excel(ms_file, sheet_name="NCP", engine="openpyxl")
    peptides = list(df[['chrom', 'start', 'end', 'strand']].itertuples(index=False, name=None))
    counter = Counter()
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(candidate_index,)) as ex:
        for hit_list in tqdm(ex.map(_worker_match_one_peptide, peptides), total=len(peptides), desc="Matching peptides (parallel)"):
            if hit_list:
                counter.update(hit_list)
    transcript_peptide_count = {tid: counter.get(tid, 0) for tid in transcript_lines.keys()}
    return transcript_lines, transcript_peptide_count

def generate_outputs(transcript_lines, transcript_peptide_count, output_prefix):
    gff_output_file = f"{output_prefix}_transcripts_with_peptides.gtf"
    with open(gff_output_file, 'w') as f:
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("##source PeptideSupportedTranscripts\n")
        for tid, lines in transcript_lines.items():
            pc = transcript_peptide_count.get(tid, 0)
            if pc > 0:
                for line in lines:
                    if (not line.startswith('#')) and ('\ttranscript\t' in line):
                        if line.endswith(';') or ';' in line:
                            line = line.rstrip() + f' peptide_count "{pc}";'
                        else:
                            line = line.rstrip() + f' peptide_count "{pc}";'
                    f.write(line + '\n')
                f.write("###\n")

    # 2) 统计
    stats_output_file = f"{output_prefix}_peptide_count_distribution.csv"
    cnt = Counter(transcript_peptide_count.values())
    pd.DataFrame(
        [{'peptide_count': k, 'transcript_count': v} for k, v in sorted(cnt.items())]
    ).to_csv(stats_output_file, index=False)

def main():
    ms_file = "/media/wanglab/caca/Eu_peptido/20251018imeta/file/00raw/Eu_sp_finally.xlsx"
    gtf_file = "/media/wanglab/caca/Eu_peptido/20251018imeta/file/00raw/output/filter_file1_nonredundant_coding_transcripts.gtf"
    output_prefix = "/media/wanglab/caca/Eu_peptido/20251018imeta/file/01new_gene/analysis_results/filter_by_ms"
    transcript_lines, transcript_peptide_count = filter_by_ms(
        ms_file, gtf_file, workers=100
    )
    generate_outputs(transcript_lines, transcript_peptide_count, output_prefix)

if __name__ == "__main__":
    main()