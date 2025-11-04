# # 转录本长度筛选
# import pandas as pd
# input_file = "/Users/lemon/Desktop/output_with_no_predict/filter_by_ms_transcripts_with_peptides.gtf"
# df = pd.read_csv(
#         input_file, sep='\t', comment='#', header=None,
#         names=['chrom','source','feature','start','end','score','strand','frame','attrs'],
#         dtype={'chrom':'string','source':'string','feature':'string','start':'int64','end':'int64','score':'string','strand':'string','frame':'string','attrs':'string'},
#         engine='c'
#     )
# feature_type = df[df['feature'] == 'transcript']
# print(len(feature_type))

# # 外显子长度统计
# import pandas as pd
# from collections import defaultdict
# import re
# RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
# RE_GENE_ID       = re.compile(r'gene_id "([^"]+)"')
# def parse_gtf_exons(gtf_file):
#     exon_lengths = []
#     with open(gtf_file, 'r') as f:
#         for line in f:
#             if line.startswith('#'):
#                 continue
#             fields = line.strip().split('\t')
#             if len(fields) < 9:
#                 continue
#             chrom, source, feature, start, end, score, strand, frame, attrs = fields
#             if feature == 'exon':
#                 start = int(start)
#                 end = int(end)
#                 exon_length = end - start + 1
#                 exon_lengths.append(exon_length)
#                 if exon_length == 1:
#                     print(f"{attrs}")
#     return exon_lengths
# def plot_exon_length_distribution(exon_lengths, output_file):
#     length_distribution = defaultdict(int)
#     for length in exon_lengths:
#         length_distribution[length] += 1
#     df = pd.DataFrame(length_distribution.items(), columns=['Exon Length', 'Frequency'])
#     df = df.sort_values(by='Exon Length')  # 按长度排序
#     df.to_csv(output_file, index=False)
#     print(f"Exon length distribution saved to {output_file}")
# def main():
#     gtf_file = '/Volumes/caca/Eu_peptido/20251018imeta/new_file_no_predict/01new_gene/analysis_results/filter_by_ms_transcripts_with_peptides.gtf'  # 替换为你的 GTF 文件路径
#     output_file = '/Volumes/caca/Eu_peptido/20251018imeta/new_file_no_predict/01new_gene/analysis_results/mapped_transcript_exon_length_distribution.csv'  # 输出的 CSV 文件路径
#     exon_lengths = parse_gtf_exons(gtf_file)
#     plot_exon_length_distribution(exon_lengths, output_file)
# if __name__ == "__main__":
#     main()

import os
import re
import pandas as pd
from datetime import datetime
from tqdm import tqdm
from bisect import bisect_left
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
RE_GENE_ID       = re.compile(r'gene_id "([^"]+)"')
_PEP_BUCKET = None
_START_BUCKET = None
_TRANSCRIPT_EXONS = None
_TID2KEY = None

def _init_worker(pep_bucket, start_bucket, transcript_exons, tid2key):
    global _PEP_BUCKET, _START_BUCKET, _TRANSCRIPT_EXONS, _TID2KEY
    _PEP_BUCKET = pep_bucket
    _START_BUCKET = start_bucket
    _TRANSCRIPT_EXONS = transcript_exons
    _TID2KEY = tid2key

def _worker_count_transcript_chunk(tid_chunk):
    counts = Counter()
    mapped_ids = set()
    map_pairs = []
    for tid in tid_chunk:
        key = _TID2KEY.get(tid)
        if key is None:
            continue
        peps = _PEP_BUCKET.get(key, [])
        if not peps:
            continue
        starts = _START_BUCKET[key]
        exons  = _TRANSCRIPT_EXONS.get(tid, [])
        if not exons:
            continue
        seed_ids = set()
        for (exon_start, exon_end) in exons:
            pos = bisect_left(starts, exon_start)
            j = pos
            while j < len(peps) and peps[j]['start'] <= exon_end:
                p = peps[j]
                if p['end'] <= exon_end:
                    pid = p['peptide_id']
                    if pid not in seed_ids:
                        counts[tid] += 1
                        seed_ids.add(pid)
                        mapped_ids.add(p['peptide_id'])
                        map_pairs.append((pid, tid))
                j += 1
    return counts, mapped_ids, map_pairs

def parse_gtf_transcript(gtf_file):
    transcript_lines = defaultdict(list)
    transcript_info = {}
    transcript_exon_region = {}
    current_transcript_id = None
    current_lines = []
    with open(gtf_file, 'r') as f:
        for line in tqdm(f, desc="Loading GTF"):
            if not line or line.startswith('#'):
                continue
            line = line.rstrip('\n')
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs = fields
            if ftype == 'transcript':
                if current_transcript_id and current_lines:
                    if len(current_lines) == 2:
                        transcript_lines[current_transcript_id] = current_lines
                tid_match = RE_TRANSCRIPT_ID.search(attrs)
                gid_match = RE_GENE_ID.search(attrs)
                tid = tid_match.group(1)
                gid = gid_match.group(1)
                current_transcript_id = tid
                current_lines = [line]
                transcript_info[tid] = {
                    'gene_id': gid,
                    'chrom': str(chrom),
                    'strand': str(strand),
                    'start': int(start),
                    'end': int(end)
                }
                transcript_exon_region.setdefault(tid, [])
            elif ftype == 'exon':
                if current_transcript_id:
                    current_lines.append(line)
                    parts = line.split('\t')
                    exon_start = int(parts[3])
                    exon_end   = int(parts[4])
                    transcript_exon_region[current_transcript_id].append((exon_start, exon_end))
    if current_transcript_id and current_lines:
        transcript_lines[current_transcript_id] = current_lines
    for tid, exons in transcript_exon_region.items():
        transcript_exon_region[tid] = sorted(exons)
    single_exon_tids = [tid for tid, exons in transcript_exon_region.items() if len(exons) == 1]
    transcript_lines_single = {tid: transcript_lines[tid] for tid in single_exon_tids if tid in transcript_lines}
    transcript_exon_region_single = {tid: transcript_exon_region[tid] for tid in single_exon_tids}
    tid2key = {}
    for tid in single_exon_tids:
        info = transcript_info.get(tid)
        if not info:
            continue
        tid2key[tid] = (info['chrom'], info['strand'])
    return transcript_lines_single, transcript_exon_region_single, tid2key

def build_peptide_buckets(ms_file):
    df = pd.read_excel(ms_file, sheet_name="NCP")
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    pep_bucket = defaultdict(list)
    for rec in df[['peptide_id', 'chrom', 'start', 'end', 'strand']].to_dict(orient='records'):
        key = (str(rec['chrom']), str(rec['strand']))
        pep_bucket[key].append({
            'start': int(rec['start']),
            'end':   int(rec['end']),
            'peptide_id': rec['peptide_id']
        })
    start_bucket = {}
    for key, lst in pep_bucket.items():
        lst.sort(key=lambda x: x['start'])
        start_bucket[key] = [p['start'] for p in lst]
    return df, pep_bucket, start_bucket

def filter_by_ms_and_extract_cds(ms_file, gtf_file, genome_file, workers=40, chunk_size=5000):
    transcript_lines, transcript_exon_regions, tid2key = parse_gtf_transcript(gtf_file)
    peptide_df, pep_bucket, start_bucket = build_peptide_buckets(ms_file)
    tids = list(transcript_lines.keys())
    chunks = [tids[i:i+chunk_size] for i in range(0, len(tids), chunk_size)]
    total_count = Counter()
    mapped_peptide_ids = set()
    mapped_pairs_all = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(pep_bucket, start_bucket, transcript_exon_regions, tid2key)) as pool:
        futs = [pool.submit(_worker_count_transcript_chunk, ch) for ch in chunks]



def main():
    ms_file = "/data/Eu/Eu_rnaseq/output_test/Eu_sp_finally.xlsx"
    gtf_file = "/data/Eu/Eu_rnaseq/output_test/file1_nonredundant_transcripts.gtf"
    genome_file = ""
    output_prefix = "/data/Eu/Eu_rnaseq/output_test/filter_by_ms"
    workers = 40
    chunk_size = 600
    filter_by_ms_and_extract_cds(ms_file, gtf_file, genome_file, output_prefix, workers, chunk_size)    

if __name__ == "__main__":
    main()