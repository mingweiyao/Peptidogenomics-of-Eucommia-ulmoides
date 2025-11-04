import os
import re
import pandas as pd
from datetime import datetime
from tqdm import tqdm
from bisect import bisect_left
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
RE_GENE_ID       = re.compile(r'gene_id "([^"]+)"')

_PEP_BUCKET = None
_START_BUCKET = None
_TID2KEY = None
_TRANSCRIPT_EXON = None

COVERAGE_THRESHOLD = 0.10

def _init_worker(pep_bucket, start_bucket, tid2key, transcript_exon):
    global _PEP_BUCKET, _START_BUCKET, _TID2KEY, _TRANSCRIPT_EXON
    _PEP_BUCKET  = pep_bucket
    _START_BUCKET = start_bucket
    _TID2KEY = tid2key
    _TRANSCRIPT_EXON = transcript_exon

def _worker_count_transcript_chunk(tid_chunk):
    counts = Counter()
    mapped_ids = set()
    map_pairs = []
    map_intervals = []
    for tid in tid_chunk:
        key = _TID2KEY.get(tid)
        peps = _PEP_BUCKET.get(key, [])
        starts = _START_BUCKET.get(key, [])
        exon_start, exon_end = _TRANSCRIPT_EXON.get(tid, (None, None))
        pos = bisect_left(starts, exon_start)
        seen_pid = set()
        j = pos
        while j < len(peps) and peps[j]['start'] <= exon_end:
            p = peps[j]
            if p['end'] <= exon_end:
                pid = p['peptide_id']
                if pid not in seen_pid:
                    counts[tid] += 1
                    seen_pid.add(pid)
                    mapped_ids.add(pid)
                    map_pairs.append((pid, tid))
                    map_intervals.append((tid, pid, p['start'], p['end']))
            j += 1
    return counts, mapped_ids, map_pairs, map_intervals

def parse_gtf_transcript(gtf_file):
    transcript_lines = defaultdict(list)
    transcript_info = {}
    transcript_exons = defaultdict(list)
    with open(gtf_file, 'r') as f:
        for line in tqdm(f, desc="Loading GTF"):
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs = fields
            if ftype == 'transcript':
                m_tid = RE_TRANSCRIPT_ID.search(attrs)
                m_gid = RE_GENE_ID.search(attrs)
                tid = m_tid.group(1)
                gid = m_gid.group(1)
                transcript_lines[tid].append(line.rstrip('\n'))
                transcript_info[tid] = {
                    'gene_id': gid,
                    'chrom': str(chrom),
                    'strand': str(strand),
                    'start': int(start),
                    'end': int(end)
                }                
            elif ftype == 'exon':
                m_tid = RE_TRANSCRIPT_ID.search(attrs)
                tid = m_tid.group(1)
                transcript_lines[tid].append(line.rstrip('\n'))
                exon_start = int(fields[3])
                exon_end   = int(fields[4])
                transcript_exons[tid].append((exon_start, exon_end))
    single_exon_tids = [tid for tid, exons in transcript_exons.items() if len(exons)==1]
    transcript_lines_single = {tid: transcript_lines[tid] for tid in single_exon_tids if tid in transcript_lines}
    transcript_exon = {tid: tuple(sorted(exs)[0]) for tid, exs in transcript_exons.items() if len(exs)==1}
    tid2key = {}
    for tid in single_exon_tids:
        info = transcript_info.get(tid)
        tid2key[tid] = (info['chrom'], info['strand'])
    return transcript_lines_single, transcript_exon, tid2key

def build_peptide_buckets(ms_file):
    df = pd.read_excel(ms_file, sheet_name="NCP")
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    pep_bucket = defaultdict(list)
    for _, rec in df.iterrows():
        key = (str(rec['chrom']), str[rec['strand']])
        record = {
            'start': int(rec['start']),
            'end': int(rec['end']),
            'peptide_id': str(rec['peptide_id']),
            'sequence': str(rec['sequence'])
        }
        pep_bucket[key].append(record)
    start_bucket = {}
    for key, lst in pep_bucket.items():
        lst.sort(key=lambda x: x['start'])
        start_bucket[key] = [p['start'] for p in lst]
    return df, pep_bucket, start_bucket

def parse_genome_file(genome_file):
    genome_sequence = {}
    for rec in SeqIO.parse(genome_file, "fasta"):
        genome_sequence[rec.id] = str(rec.seq).upper()
    return genome_sequence

def filter_and_extract_cds(ms_file, gtf_file, genome_file, workers, chunk_size):
    transcript_lines, transcript_exon, tid2key = parse_gtf_transcript(gtf_file)
    peptide_df, pep_bucket, start_bucket = build_peptide_buckets(ms_file)
    tids = list(transcript_lines.keys())
    chunks = [tids[i:i+chunk_size] for i in range(0, len(tids), chunk_size)]
    total_counts = Counter()
    mapped_peptide_ids = set()
    mapped_pairs_all = []
    per_tid_intervals = defaultdict(list)
    
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(pep_bucket, start_bucket, tid2key, transcript_exon)) as pool:
        futs = [pool.submit(_worker_count_transcript_chunk, ch) for ch in chunks]
        for fu in tqdm(as_completed(futs), total=len(futs), desc="Mapping: "):
            c_part, ids_part, pairs_part, intervals_part = fu.result()
            total_counts.update(c_part)
            mapped_peptide_ids |= ids_part
            mapped_pairs_all.extend(pairs_part)
            for tid, pid, s, e in intervals_part:
                per_tid_intervals[tid].append((s, e, pid))

    for tid, intervals in per_tid_intervals.items():
        frame_filter = {}
        for s,e,pid in intervals:
            test_len = (s - transcript_exon[tid][0]) % 3
            frame_filter[test_len] += 1
        

    def union_len(intervals):
        if not intervals:
            return 0
        ivs = sorted([(s,e) for s,e,_ in intervals])
        total = 0
        cur_s, cur_e = ivs[0]
        for s,e in ivs[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                total += (cur_e - cur_s + 1)
                cur_s, cur_e = s, e
        total += (cur_e - cur_s + 1)
        return total
    transcript_len = {tid: (transcript_exon[tid][1] - transcript_exon[tid][0] + 1) for tid in tids}
    tid2_cov = {}
    for tid in tids:
        cov_nt = union_len(per_tid_intervals[tid])
        L = transcript_len[tid]
        tid2_cov[tid] = cov_nt / L
    kept_tids = [tid for tid in tids if tid2_cov.get(tid, 0.0) >= COVERAGE_THRESHOLD]

    genome = parse_genome_file(genome_file)
    tid2_cds_seq = {}
    tid2_cds_span = {}
    tid2_pep_minmax = {}
    tids_start_is_m = {}
    tid2_reason = {}
    for tid in kept_tids:
        chrom, strand = tid2key[tid]
        exon_start, exon_end = transcript_exon[tid]
        intervals = per_tid_intervals[tid]
        if not intervals:
            tid2_reason[tid] = "无肽段或被覆盖率阈值过滤"
            continue
        pep_min = min(s for s,_,_ in intervals)
        pep_max = max(e for _,e,_ in intervals)
        tid2_pep_minmax[tid] = (pep_min, pep_max)
                           










def main():
    ms_file     = "/data/Eu/Eu_rnaseq/output_test/Eu_sp_finally.xlsx"
    gtf_file    = "/data/Eu/Eu_rnaseq/output_test/file1_nonredundant_transcripts.gtf"
    genome_file = "/data/Eu/Eu_rnaseq/genome.fa"
    output_prefix = "/data/Eu/Eu_rnaseq/output_test/filter_by_ms"
    workers = 40
    chunk_size = 600
    res = filter_and_extract_cds(
        ms_file=ms_file,
        gtf_file=gtf_file,
        genome_file=genome_file,
        workers=workers,
        chunk_size=chunk_size
    )

if __name__ == "__main__":
    main()