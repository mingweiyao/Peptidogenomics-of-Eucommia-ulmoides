
import os
import re
import pandas as pd
from datetime import datetime
from tqdm import tqdm
from bisect import bisect_left
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---------------- Regex ----------------
RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
RE_GENE_ID       = re.compile(r'gene_id "([^"]+)"')

# --------------- Globals for workers (read-only) ---------------
_PEP_BUCKET = None          # {(chrom,strand): [{'start','end','peptide_id','chrom','strand'}, ...] 按start排序}
_START_BUCKET = None        # {(chrom,strand): [start,...]} 与上同步
_TRANSCRIPT_EXONS = None    # {tid: [(exon_start,exon_end), ...] 已排序}
_TID2KEY = None             # {tid: (chrom,strand)}

def _init_worker(pep_bucket, start_bucket, transcript_exons, tid2key):
    """Initializer for worker processes."""
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

def parse_gtf_transcripts(gtf_file):
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

    tid2key = {}
    for tid, info in transcript_info.items():
        key = (info['chrom'], info['strand'])
        tid2key[tid] = key

    return transcript_lines, transcript_exon_region, tid2key

def build_peptide_buckets(ms_file):
    df = pd.read_excel(ms_file, sheet_name="NCP", engine="openpyxl")
    df['start'] = df['start'].astype(int)
    df['end']   = df['end'].astype(int)
    pep_bucket = defaultdict(list)
    for rec in df[['peptide_id','chrom','start','end','strand']].to_dict(orient='records'):
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

def generate_outputs(transcript_lines, transcript_peptide_count, output_prefix):
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    gff_output_file = f"{output_prefix}_transcripts_with_peptides.gtf"
    with open(gff_output_file, 'w') as f:
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("##source PeptideSupportedTranscripts\n")
        for tid, lines in transcript_lines.items():
            pc = transcript_peptide_count.get(tid, 0)
            if pc > 0:
                for line in lines:
                    if (not line.startswith('#')) and ('\ttranscript\t' in line):
                        line = line.rstrip() + f' peptide_count "{pc}";'
                    f.write(line + '\n')
    print(f"[OK] 转录本（带肽段计数）GTF: {gff_output_file}")

    output_unmapped_tx = f"{output_prefix}_transcripts_without_peptides.gtf"
    with open(output_unmapped_tx, 'w') as f:
        f.write("##source TranscriptsNoPeptides\n")
        for tid, lines in transcript_lines.items():
            pc = transcript_peptide_count.get(tid, 0)
            if pc == 0:
                for line in lines:
                    if (not line.startswith('#')) and ('\ttranscript\t' in line):
                        line = line.rstrip() + f' peptide_count "{pc}";'
                    f.write(line + '\n')
    print(f"[OK] 无肽段支持的转录本 GTF: {output_unmapped_tx}")

    stats_output_file = f"{output_prefix}_peptide_count_distribution.csv"
    cnt = Counter(transcript_peptide_count.values())
    pd.DataFrame(
        [{'peptide_count': k, 'transcript_count': v} for k, v in sorted(cnt.items())]
    ).to_csv(stats_output_file, index=False)
    print(f"[OK] 肽段计数分布 CSV: {stats_output_file}")

def export_peptide_mapping_excels(peptide_df, mapped_pairs_all, unmapped_peptides, output_prefix):
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    map_df = pd.DataFrame(mapped_pairs_all, columns=['peptide_id', 'transcript_id'])
    if not map_df.empty:
        agg = (map_df.groupby('peptide_id')['transcript_id']
                    .agg(lambda x: ';'.join(sorted(set(x))))
                    .reset_index()
                    .rename(columns={'transcript_id': 'matched_transcripts'}))
        agg['matched_count'] = agg['matched_transcripts'].apply(lambda s: 0 if not s else s.count(';') + 1)
        mapped_detail = peptide_df.merge(agg, on='peptide_id', how='inner')
    else:
        mapped_detail = peptide_df.head(0).copy()
        mapped_detail['matched_transcripts'] = []
        mapped_detail['matched_count'] = []
    unmapped_detail = pd.DataFrame(unmapped_peptides)
    mapped_xlsx   = f"{output_prefix}_mapped_peptides.xlsx"
    unmapped_xlsx = f"{output_prefix}_unmapped_peptides.xlsx"
    with pd.ExcelWriter(mapped_xlsx, engine='openpyxl') as w:
        mapped_detail.to_excel(w, index=False, sheet_name='mapped')
    with pd.ExcelWriter(unmapped_xlsx, engine='openpyxl') as w:
        unmapped_detail.to_excel(w, index=False, sheet_name='unmapped')
    print(f"[OK] 映射肽段明细: {mapped_xlsx}")
    print(f"[OK] 未映射肽段明细: {unmapped_xlsx}")

def filter_by_ms(ms_file, gtf_file, workers, chunk_size):
    transcript_lines, transcript_exon_region, tid2key = parse_gtf_transcripts(gtf_file)
    peptide_df, pep_bucket, start_bucket = build_peptide_buckets(ms_file)
    tids = list(transcript_lines.keys())
    chunks = [tids[i:i+chunk_size] for i in range(0, len(tids), chunk_size)]

    total_counts = Counter()
    mapped_peptide_ids = set()
    mapped_pairs_all = []

    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                            initargs=(pep_bucket, start_bucket, transcript_exon_region, tid2key)) as pool:
        futs = [pool.submit(_worker_count_transcript_chunk, ch) for ch in chunks]
        for fu in tqdm(as_completed(futs), total=len(futs), desc="Counting (parallel)"):
            c_part, ids_part, pairs_part = fu.result()
            total_counts.update(c_part)
            mapped_peptide_ids |= ids_part
            mapped_pairs_all.extend(pairs_part)

    transcript_peptide_count = {tid: total_counts.get(tid, 0) for tid in tids}
    unmapped_peptides = [
        rec for rec in peptide_df[['peptide_id','chrom','start','end','strand']].to_dict(orient='records')
        if rec['peptide_id'] not in mapped_peptide_ids
    ]
    return transcript_lines, transcript_peptide_count, peptide_df, mapped_peptide_ids, mapped_pairs_all, unmapped_peptides

def main():
    ms_file = "/data/Eu/Eu_rnaseq/output_test/Eu_sp_finally.xlsx"
    gtf_file = "/data/Eu/Eu_rnaseq/output_test/file1_nonredundant_transcripts.gtf"
    output_prefix = "/data/Eu/Eu_rnaseq/output_test/filter_by_ms"
    workers = 40
    chunk_size = 600

    (transcript_lines,
     transcript_peptide_count,
     peptide_df,
     mapped_peptide_ids,
     mapped_pairs_all,
     unmapped_peptides) = filter_by_ms(ms_file, gtf_file, workers=workers, chunk_size=chunk_size)
    generate_outputs(transcript_lines, transcript_peptide_count, output_prefix)
    export_peptide_mapping_excels(peptide_df, mapped_pairs_all, unmapped_peptides, output_prefix)

if __name__ == "__main__":
    main()
