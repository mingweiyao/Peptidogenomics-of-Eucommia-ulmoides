import os
import re
import pandas as pd
from datetime import datetime
from tqdm import tqdm
from bisect import bisect_left
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
RE_GENE_ID       = re.compile(r'gene_id "([^"]+)"')

_PEP_BUCKET = None
_START_BUCKET = None
_TID2KEY = None
_TRANSCRIPT_EXON = None

COVERAGE_THRESHOLD = 0.10
FASTA_LINE_WRAP = 100

def _init_worker(pep_bucket, start_bucket, tid2key, transcript_exon):
    global _PEP_BUCKET, _START_BUCKET, _TID2KEY, _TRANSCRIPT_EXON
    _PEP_BUCKET  = pep_bucket
    _START_BUCKET = start_bucket
    _TID2KEY = tid2key
    _TRANSCRIPT_EXON = transcript_exon

def _worker_count_transcript_chunk(tid_chunk):
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
                    map_intervals.append((tid, pid, p['start'], p['end']))
                    seen_pid.add(pid)
            j += 1
    return map_intervals

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
        key = (str(rec['chrom']), str(rec['strand']))
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

START_CODEN = {'ATG'}
STOP_CODEN  = {'TAA','TAG','TGA'}
START_CODEN_RES = {'CAT'}
STOP_CODEN_RES = {'TTA', 'CTA', 'TCA'}
def find_cds_seq(genome, chrom, strand, pep_min, pep_max, exon_start, exon_end):
    seq = genome[chrom]
    cds_start = None
    cds_end = None
    # python splic [start, end)
    if strand == '+':
        cod = seq[pep_min-1:pep_min+2]
        if cod in START_CODEN:
            cds_start = pep_min
        else:
            p = pep_min - 3
            while p >= exon_start:
                cod = seq[p-1:p+2]
                if cod in START_CODEN:
                    cds_start = p - 1
                    break
                p -= 3
        p = pep_max
        while p < exon_end:
            cod = seq[p:p+3]
            if cod in STOP_CODEN:
                cds_end = p
                break
            p += 3
    else:
        cod = seq[pep_max-3:pep_max]
        if cod in START_CODEN_RES:
            cds_end = pep_max
        else:
            p = pep_max + 3
            while p <= exon_end:
                cod = seq[p-3: p]
                if cod in START_CODEN_RES:
                    cds_end = p
                    break
                p += 3
        p = pep_min-3
        while p >= exon_start:
            cod = seq[p-1: p+2]
            if cod in STOP_CODEN_RES:
                cds_start = p - 1
                break
            p -= 3
    return seq[cds_start:cds_end]

def filter_and_extract_cds(ms_file, gtf_file, genome_file, workers, chunk_size):
    transcript_lines, transcript_exon, tid2key = parse_gtf_transcript(gtf_file)
    peptide_df, pep_bucket, start_bucket = build_peptide_buckets(ms_file)
    tids = list(transcript_lines.keys())
    chunks = [tids[i:i+chunk_size] for i in range(0, len(tids), chunk_size)]
    per_tid_intervals = defaultdict(list)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(pep_bucket, start_bucket, tid2key, transcript_exon)) as pool:
        futs = [pool.submit(_worker_count_transcript_chunk, ch) for ch in chunks]
        for fu in tqdm(as_completed(futs), total=len(futs), desc="Mapping: "):
            intervals_part = fu.result()
            for tid, pid, s, e in intervals_part:
                per_tid_intervals[tid].append((s, e, pid))
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
    # update transcript validated by peptide in same frame
    filt_intervals = {}
    for tid, intervals in per_tid_intervals.items():
        if not intervals or tid not in transcript_exon or tid not in tid2key:
            filt_intervals[tid] = []
            continue
        exon_start, exon_end = transcript_exon[tid]
        _, strand = tid2key[tid]
        frames_info = []
        frame_cnt = Counter()
        for s, e, pid in intervals:
            off = (s - exon_start) if strand == '+' else (exon_end - e)
            if off < 0:
                continue
            fr = off % 3
            frames_info.append((fr, s, e, pid))
            frame_cnt[fr] += 1
        if not frame_cnt:
            filt_intervals[tid] = []
            continue
        max_cnt = max(frame_cnt.values())
        candidate_frames = [fr for fr, c in frame_cnt.items() if c == max_cnt]
        if len(candidate_frames) > 1:
            cover_by_frame = {}
            for fr in candidate_frames:
                iv = [(s, e, pid) for (f, s, e, pid) in frames_info if f == fr]
                cover_by_frame[fr] = union_len(iv)
            max_cover = max(cover_by_frame.values())
            candidate_frames = [fr for fr, cov in cover_by_frame.items() if cov == max_cover]
        best_frame = min(candidate_frames)
        filt_intervals[tid] = [(s, e, pid) for (fr, s, e, pid) in frames_info if fr == best_frame]
    per_tid_intervals = filt_intervals
    # filter transcript that coverage > 10%
    transcript_len = {tid: (transcript_exon[tid][1] - transcript_exon[tid][0] + 1) for tid in tids}
    tid2_cov = {}
    for tid in tids:
        cov_nt = union_len(per_tid_intervals[tid])
        L = transcript_len[tid]
        tid2_cov[tid] = cov_nt / L
    kept_tids = [tid for tid in tids if tid2_cov.get(tid, 0.0) >= COVERAGE_THRESHOLD]
    kept_tids_intervals = {}
    for tid in kept_tids:
        if tid in filt_intervals:
            kept_tids_intervals[tid] = filt_intervals[tid]
    per_tid_intervals = kept_tids_intervals
    total_counts = Counter()
    mapped_peptide_ids = set()
    mapped_pairs_all = []
    for tid, lst in per_tid_intervals.items():
        pid_set = {pid for _, _, pid in lst}
        total_counts[tid] = len(pid_set)
        mapped_peptide_ids |= pid_set
        mapped_pairs_all.extend((pid, tid) for _, _, pid in lst)
    # extract cds of transcript
    genome = parse_genome_file(genome_file)
    tid2_cds_seq = {}
    for tid in kept_tids:
        chrom, strand = tid2key[tid]
        exon_start, exon_end = transcript_exon[tid]
        intervals = per_tid_intervals[tid]
        if not intervals:
            continue
        pep_min = min(s for s,_,_ in intervals)
        pep_max = max(e for _,e,_ in intervals)
        cds_seq = find_cds_seq(genome, chrom, strand, pep_min, pep_max, exon_start, exon_end)
        tid2_cds_seq[tid] = cds_seq
    transcript_peptide_count = {tid: total_counts.get(tid, 0) for tid in tids}
    return {
        'transcript_lines': transcript_lines,
        'transcript_exon': transcript_exon,
        'tid2key': tid2key,
        'peptide_df': peptide_df,
        'mapped_pairs_all': mapped_pairs_all,
        'transcript_peptide_count': transcript_peptide_count,
        'tid2_cov': tid2_cov,
        'kept_tids': set(kept_tids),
        'tid2_cds_seq': tid2_cds_seq
    }

def generate_output(res, output_prefix):
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    # output mapped and unmapped peptide file
    peptide_df = res.get("peptide_df")
    mapped_pairs_all = res.get("mapped_pairs_all", [])
    mapped_ids = {str(pid) for pid, _ in mapped_pairs_all}
    peptide_df = peptide_df.copy()
    peptide_df["peptide_id"] = peptide_df["peptide_id"].astype(str)
    mapped_peptide_df = peptide_df[peptide_df["peptide_id"].isin(mapped_ids)].reset_index(drop=True)
    unmapped_peptides_df = peptide_df[~peptide_df["peptide_id"].isin(mapped_ids)].reset_index(drop=True)
    mapped_peptide_df.to_csv(f"{output_prefix}.mapped_peptides.csv", index=False)
    unmapped_peptides_df.to_csv(f"{output_prefix}.unmapped_peptides.csv", index=False)
    # output transcript GTF file
    transcript_lines = res.get("transcript_lines")
    tid2_cds_seq = res.get("tid2_cds_seq")
    tid2_cov = res.get("tid2_cov")
    transcript_peptide_count = res.get("transcript_peptide_count")
    tids = list(tid2_cds_seq.keys())
    output_transcript_lines = defaultdict(list)
    for tid in tids:
        output_transcript_lines[tid] = transcript_lines[tid]
    gtf_supported = f"{output_prefix}_transcripts_with_peptides.gtf"
    with open(gtf_supported, 'w') as f:
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("##source PeptideSupportedTranscripts (coverage>=threshold AND cds_extracted)\n")
        for tid, lines in output_transcript_lines.items():
            cov = tid2_cov.get(tid, 0.0)
            pc = transcript_peptide_count.get(tid, 0)
            for line in lines:
                out = line
                if (not line.startswith('#')) and ('\ttranscript\t' in line):
                    out = line.rstrip() + f' peptide_count "{pc}"; coverage "{cov:.4f}";'
                f.write(out + '\n')
    # output transcript cds and cds translate file
    gtf_cds = f"{output_prefix}_gtf_cds.fasta"
    gtf_cds_aa = f"{output_prefix}_gtf_cds_amino_acids.fasta"
    tid2key = res.get("tid2key")
    records_cds, records_cds_aa = [], []
    for tid, seq_str in tid2_cds_seq.items():
        if not seq_str:
            continue
        seq_obj = Seq(seq_str)
        strand = tid2key[tid][1]
        chrom = tid2key[tid][0]
        records_cds.append(SeqRecord(seq_obj, id=tid, description=f"{chrom}\t{strand}"))
        aa = seq_obj.translate() if strand == '+' else seq_obj.reverse_complement().translate()
        records_cds_aa.append(SeqRecord(aa, id=tid, description=f"{chrom}\t{strand}"))
    SeqIO.write(records_cds, gtf_cds, "fasta")
    SeqIO.write(records_cds_aa, gtf_cds_aa, "fasta")
    # output coverage and peptide count file
    coverage_stats = []
    for tid in tids:
        cov = tid2_cov.get(tid, 0.0)
        pc = transcript_peptide_count.get(tid, 0)
        coverage_stats.append({'transcript_id': tid, 'coverage': cov, 'peptide_count': pc})
    coverage_df = pd.DataFrame(coverage_stats)
    coverage_df.to_csv(f"{output_prefix}_coverage_stats.csv", index=False)
    peptide_count_stats = []
    for tid, count in transcript_peptide_count.items():
        peptide_count_stats.append({'transcript_id': tid, 'peptide_count': count})
    peptide_count_df = pd.DataFrame(peptide_count_stats)
    peptide_count_df.to_csv(f"{output_prefix}_peptide_count_stats.csv", index=False)

def main():
    ms_file     = "/data/Eu/Eu_rnaseq/output_test/Eu_sp_finally.xlsx"
    gtf_file    = "/data/Eu/Eu_rnaseq/output_test/file1_nonredundant_transcripts.gtf"
    genome_file = "/data/Eu/Eu_rnaseq/genome.fa"
    output_prefix = "/data/Eu/Eu_rnaseq/output_test/filter_by_ms"
    workers = 40
    chunk_size = 1000
    res = filter_and_extract_cds(
        ms_file=ms_file,
        gtf_file=gtf_file,
        genome_file=genome_file,
        workers=workers,
        chunk_size=chunk_size
    )
    generate_output(res, output_prefix)
if __name__ == "__main__":
    main()