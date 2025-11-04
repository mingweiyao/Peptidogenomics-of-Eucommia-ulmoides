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
from Bio import SeqIO

RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
RE_GENE_ID       = re.compile(r'gene_id "([^"]+)"')

# ========== 全局（仅 worker 用） ==========
_PEP_BUCKET = None          # (chrom,strand) -> list of {'start','end','peptide_id', 'pep_seq'(可无)}
_START_BUCKET = None        # (chrom,strand) -> sorted list of starts
_TID2KEY = None             # tid -> (chrom,strand)
_TID2EXON = None            # tid -> (exon_start, exon_end)

# ========== 参数（可按需调整） ==========
UPSTREAM_ATG_SEARCH_BP = 3000   # 非 M 起始时，上游搜 ATG 的最大范围
FASTA_LINE_WRAP = 60            # 输出FASTA换行宽度
COVERAGE_THRESHOLD = 0.10       # 覆盖率阈值：<10% 的转录本剔除

# ----------------- Worker 初始化 -----------------
def _init_worker(pep_bucket, start_bucket, tid2key, tid2exon):
    global _PEP_BUCKET, _START_BUCKET, _TID2KEY, _TID2EXON
    _PEP_BUCKET  = pep_bucket
    _START_BUCKET = start_bucket
    _TID2KEY = tid2key
    _TID2EXON = tid2exon

# ----------------- 子进程：统计与映射 -----------------
def _worker_count_transcript_chunk(tid_chunk):
    """
    返回：
      counts: Counter[tid] = 匹配肽段个数（去重按 peptide_id）
      mapped_ids: set(peptide_id)
      map_pairs: list[(peptide_id, tid)]
      map_intervals: list[(tid, peptide_id, start, end)]
    """
    counts = Counter()
    mapped_ids = set()
    map_pairs = []
    map_intervals = []   # 供主进程计算覆盖率/最小最大区间

    for tid in tid_chunk:
        key = _TID2KEY.get(tid)
        if key is None:
            continue
        peps = _PEP_BUCKET.get(key, [])
        starts = _START_BUCKET.get(key, [])
        if not peps or not starts:
            continue

        exon_start, exon_end = _TID2EXON.get(tid, (None, None))
        if exon_start is None:
            continue

        pos = bisect_left(starts, exon_start)
        seen_pid = set()

        j = pos
        while j < len(peps) and peps[j]['start'] <= exon_end:
            p = peps[j]
            # 仅接受完全包含于外显子的肽段
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
    single_exon_tids = [tid for tid, exs in transcript_exons.items() if len(exs) == 1]
    transcript_lines_single = {tid: transcript_lines[tid] for tid in single_exon_tids if tid in transcript_lines}
    tid2exon = {tid: tuple(sorted(exs)[0]) for tid, exs in transcript_exons.items() if len(exs) == 1}
    tid2key = {}
    for tid in single_exon_tids:
        info = transcript_info.get(tid)
        if info:
            tid2key[tid] = (info['chrom'], info['strand'])
    return transcript_lines_single, tid2exon, tid2key

# ----------------- 读取肽段并分桶 -----------------
def detect_pep_seq_column(df: pd.DataFrame):
    for c in ['peptide_seq', 'sequence', 'peptide', 'Peptide', 'AASequence']:
        if c in df.columns:
            return c
    return None

def build_peptide_buckets(ms_file):
    df = pd.read_excel(ms_file, sheet_name="NCP")
    df['start'] = df['start'].astype(int)
    df['end']   = df['end'].astype(int)

    pep_seq_col = detect_pep_seq_column(df)
    pep_seq_map = {}
    if pep_seq_col:
        # 预构建 id->seq 映射，避免循环中频繁 df 查询
        tmp = df[['peptide_id', pep_seq_col]].drop_duplicates()
        pep_seq_map = dict(zip(tmp['peptide_id'].astype(str), tmp[pep_seq_col].astype(str)))

    pep_bucket = defaultdict(list)
    for rec in df[['peptide_id', 'chrom', 'start', 'end', 'strand']].to_dict(orient='records'):
        key = (str(rec['chrom']), str(rec['strand']))
        record = {
            'start': int(rec['start']),
            'end': int(rec['end']),
            'peptide_id': str(rec['peptide_id'])
        }
        if pep_seq_col:
            record['pep_seq'] = pep_seq_map.get(record['peptide_id'])
        pep_bucket[key].append(record)

    start_bucket = {}
    for key, lst in pep_bucket.items():
        lst.sort(key=lambda x: x['start'])
        start_bucket[key] = [p['start'] for p in lst]
    return df, pep_bucket, start_bucket, pep_seq_col

# ----------------- FASTA 读取 -----------------
def parse_genome_file(genome_file):
    if not genome_file:
        raise ValueError("必须提供 genome FASTA 文件路径。")
    genome_sequence = {}
    for rec in SeqIO.parse(genome_file, "fasta"):
        genome_sequence[rec.id] = str(rec.seq).upper()
    if not genome_sequence:
        raise ValueError("未能从 FASTA 读取到任何序列。请检查路径与文件格式。")
    return genome_sequence

# ----------------- 工具：序列、反向互补、翻译、寻找起止密码子 -----------------
COMPLEMENT = str.maketrans('ACGTNacgtn','TGCANtgcan')
STOP_CODONS = {'TAA','TAG','TGA'}

def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]

def get_genomic_seq(genome, chrom, start, end, strand):
    """
    1-based, inclusive coordinates.
    返回与转录本方向一致的序列：+链为原序，-链为反向互补。
    """
    seq = genome[chrom][start-1:end]
    return seq if strand == '+' else revcomp(seq)

def find_first_inframe_stop(seq, frame=0, start_offset_nt=0):
    """
    在给定读码框（0/1/2）下，从 start_offset_nt 开始，每3nt找第一个终止密码子。
    返回：相对于 seq 起点的终止密码子 **首位**的 nt 下标（0-based），若无则返回 None。
    """
    i = start_offset_nt
    while i % 3 != frame:
        i += 1
        if i >= len(seq):
            return None
    while i + 3 <= len(seq):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            return i
        i += 3
    return None

def find_nearest_upstream_atg(seq, frame=0, ref_pos_nt=0, max_up_bp=3000):
    """
    在同一读码框下，向上游（负向）搜索最近 ATG。ref_pos_nt 为参考位置（例如“最靠5′肽段”的起点偏移）。
    返回：ATG 首位的 0-based nt 下标；若未找到，返回 None。
    """
    start = max(0, ref_pos_nt - max_up_bp)
    i = ref_pos_nt
    while i % 3 != frame:
        i -= 1
        if i < 0:
            return None
    while i >= start:
        if i + 3 <= len(seq) and seq[i:i+3] == 'ATG':
            return i
        i -= 3
    return None

def filter_and_extract_cds(ms_file, gtf_file, genome_file, workers=40, chunk_size=5000):
    transcript_lines, tid2exon, tid2key = parse_gtf_transcript(gtf_file)

    # 2) 肽段桶
    peptide_df, pep_bucket, start_bucket, pep_seq_col = build_peptide_buckets(ms_file)

    # 3) 多进程映射（完全包含）
    tids = list(transcript_lines.keys())
    chunks = [tids[i:i+chunk_size] for i in range(0, len(tids), chunk_size)]

    total_counts = Counter()
    mapped_peptide_ids = set()
    mapped_pairs_all = []
    per_tid_intervals = defaultdict(list)  # tid -> list of (start,end,peptide_id)

    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(pep_bucket, start_bucket, tid2key, tid2exon)) as pool:
        futs = [pool.submit(_worker_count_transcript_chunk, ch) for ch in chunks]
        for fu in tqdm(as_completed(futs), total=len(futs), desc="Mapping (parallel)"):
            c_part, ids_part, pairs_part, intervals_part = fu.result()
            total_counts.update(c_part)
            mapped_peptide_ids |= ids_part
            mapped_pairs_all.extend(pairs_part)
            for tid, pid, s, e in intervals_part:
                per_tid_intervals[tid].append((s, e, pid))

    # 计算覆盖率（并区间长度 / 转录本长度），并过滤 <COVERAGE_THRESHOLD
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

    transcript_len = {tid: (tid2exon[tid][1] - tid2exon[tid][0] + 1) for tid in tids}
    tid2_cov = {}
    for tid in tids:
        cov_nt = union_len(per_tid_intervals[tid])
        L = transcript_len[tid]
        tid2_cov[tid] = (cov_nt / L) if L > 0 else 0.0

    kept_tids = [tid for tid in tids if tid2_cov.get(tid, 0.0) >= COVERAGE_THRESHOLD]

    # 读基因组序列
    genome = parse_genome_file(genome_file)

    # 4) 按规则推断并提取 CDS（仅对 kept_tids）
    tid2_cds_seq = {}       # tid -> CDS DNA 序列（5'->3'，与转录本方向一致）
    tid2_cds_span = {}      # tid -> (cds_genome_start, cds_genome_end) 在基因组上的范围（1-based, inclusive）
    tid2_pep_minmax = {}    # tid -> (pep_min_start, pep_max_end)（基因组坐标，小到大）
    tid2_start_is_m = {}    # tid -> bool（“最靠5′的肽段”是否以M起始）
    tid2_reason = {}        # tid -> 提示（无法提取时说明）

    for tid in kept_tids:
        chrom, strand = tid2key[tid]
        exon_start, exon_end = tid2exon[tid]
        intervals = per_tid_intervals[tid]
        if not intervals:
            tid2_reason[tid] = "无肽段或被覆盖率阈值过滤"
            continue

        pep_min = min(s for s,_,_ in intervals)
        pep_max = max(e for _,e,_ in intervals)
        tid2_pep_minmax[tid] = (pep_min, pep_max)

        # 转录本方向序列：+链为 exon_start..exon_end，-链为其反向互补
        tx_seq = get_genomic_seq(genome, chrom, exon_start, exon_end, strand)

        # 将每条肽段的“转录本起点偏移”求出，用来确定最靠5′的那条（基准肽段）
        if strand == '+':
            off_starts = [(s - exon_start, pid) for (s, e, pid) in intervals]
        else:
            # 负链：转录本 5′->3′ 对应基因组从 exon_end 向 exon_start
            # 肽段“转录本起点”应是其基因组 end 对应的偏移
            off_starts = [(exon_end - e, pid) for (s, e, pid) in intervals]

        off_starts = [x for x in off_starts if 0 <= x[0] < len(tx_seq)]
        if not off_starts:
            tid2_reason[tid] = "肽段偏移计算异常"
            continue

        # 选出转录本方向最靠 5′ 的肽段（偏移最小）
        off_first, pid_first = min(off_starts, key=lambda x: x[0])

        # 同时算出“最靠 3′ 的肽段末端偏移”，用于确保 CDS 覆盖全部肽段
        if strand == '+':
            off_max = pep_max - exon_start
        else:
            off_max = exon_end - pep_min
        if not (0 <= off_first < len(tx_seq) and 0 <= off_max < len(tx_seq)):
            tid2_reason[tid] = "肽段偏移越界（超出exon）"
            continue

        # frame 由“最靠5′肽段”的起点偏移决定
        frame = off_first % 3

        # 仅检查“最靠5′肽段”是否以 M 起始
        start_with_M = False
        if 'peptide_id' in peptide_df.columns:
            pep_seq_col = next((c for c in ['peptide_seq','sequence','peptide','Peptide','AASequence'] if c in peptide_df.columns), None)
            if pep_seq_col:
                sub = peptide_df.loc[peptide_df['peptide_id'].astype(str) == str(pid_first), pep_seq_col]
                if not sub.empty:
                    start_with_M = str(sub.iloc[0]).upper().startswith('M')

        tid2_start_is_m[tid] = start_with_M

        # 决定 CDS 起点：关键改动 —— 若非M起始且上游找不到ATG，则移除此转录本
        if start_with_M:
            cds_start_off = off_first
        else:
            atg_off = find_nearest_upstream_atg(
                tx_seq, frame=frame, ref_pos_nt=off_first, max_up_bp=UPSTREAM_ATG_SEARCH_BP
            )
            if atg_off is None:
                tid2_reason[tid] = "非M起始且上游未找到同框ATG，移除该转录本"
                continue
            cds_start_off = atg_off

        # 从 cds_start_off 起，同框向下找第一个覆盖到“最靠3′的肽段末端”的终止密码子
        first_stop_off = find_first_inframe_stop(tx_seq, frame=frame, start_offset_nt=cds_start_off)
        if first_stop_off is None or first_stop_off < off_max:
            i = max(off_max, cds_start_off)
            while i % 3 != frame:
                i += 1
                if i >= len(tx_seq):
                    break
            stop_off = None
            while i + 3 <= len(tx_seq):
                if tx_seq[i:i+3] in STOP_CODONS:
                    stop_off = i
                    break
                i += 3
            first_stop_off = stop_off

        if first_stop_off is None or first_stop_off <= cds_start_off:
            tid2_reason[tid] = "未找到同框终止密码子覆盖肽段末端"
            continue

        # CDS 序列（不含终止密码子），保证为3的倍数
        cds_seq = tx_seq[cds_start_off:first_stop_off]
        if len(cds_seq) == 0 or len(cds_seq) % 3 != 0:
            cds_seq = cds_seq[:len(cds_seq)//3*3]
        if not cds_seq:
            tid2_reason[tid] = "CDS 序列为空"
            continue

        # 映射回基因组坐标（1-based, inclusive）
        if strand == '+':
            cds_start_genome = exon_start + cds_start_off
            cds_end_genome   = exon_start + (first_stop_off - 1)
        else:
            cds_start_genome = exon_end - cds_start_off
            cds_end_genome   = exon_end - (first_stop_off - 1)
        cds_g_start = min(cds_start_genome, cds_end_genome)
        cds_g_end   = max(cds_start_genome, cds_end_genome)

        tid2_cds_seq[tid] = cds_seq
        tid2_cds_span[tid] = (cds_g_start, cds_g_end)

    # 5) 组织输出所需结构
    unmapped_peptides = [
        rec for rec in peptide_df[['peptide_id','chrom','start','end','strand']].to_dict(orient='records')
        if str(rec['peptide_id']) not in mapped_peptide_ids
    ]
    transcript_peptide_count = {tid: total_counts.get(tid, 0) for tid in tids}

    return {
        'transcript_lines': transcript_lines,
        'tid2exon': tid2exon,
        'tid2key': tid2key,
        'peptide_df': peptide_df,
        'mapped_pairs_all': mapped_pairs_all,
        'unmapped_peptides': unmapped_peptides,
        'transcript_peptide_count': transcript_peptide_count,
        'tid2_cov': tid2_cov,
        'kept_tids': set(kept_tids),          # 覆盖率达标
        'tid2_cds_seq': tid2_cds_seq,         # 真正成功提取到CDS的子集
        'tid2_cds_span': tid2_cds_span,
        'tid2_pep_minmax': tid2_pep_minmax,
        'tid2_start_is_m': tid2_start_is_m,
        'tid2_reason': tid2_reason,
    }

# ----------------- 导出：GTF/FASTA/CSV/Excel -----------------
def write_wrapped_fasta(handle, header, seq, width=60):
    handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i+width] + "\n")

def generate_outputs(res, output_prefix):
    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    transcript_lines = res['transcript_lines']
    transcript_peptide_count = res['transcript_peptide_count']
    tid2_cov = res['tid2_cov']
    kept_tids = res['kept_tids']
    cds_tids = set(res['tid2_cds_seq'].keys())  # 成功提取 CDS 的转录本

    # 支持且成功提取CDS的转录本
    gtf_supported = f"{output_prefix}_transcripts_with_peptides.gtf"
    with open(gtf_supported, 'w') as f:
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("##source PeptideSupportedTranscripts (coverage>=threshold AND cds_extracted)\n")
        for tid, lines in transcript_lines.items():
            pc = transcript_peptide_count.get(tid, 0)
            cov = tid2_cov.get(tid, 0.0)
            if (tid in kept_tids) and (tid in cds_tids) and pc > 0:
                for line in lines:
                    out = line
                    if (not line.startswith('#')) and ('\ttranscript\t' in line):
                        out = line.rstrip() + f' peptide_count "{pc}"; coverage "{cov:.4f}"; cds_extracted "1";'
                    f.write(out + '\n')
    print(f"[OK] 转录本（覆盖率达标且成功提取CDS）GTF: {gtf_supported}")

    # 因覆盖率<阈值被剔除的转录本
    gtf_unsupported_cov = f"{output_prefix}_transcripts_removed_by_coverage.gtf"
    with open(gtf_unsupported_cov, 'w') as f:
        f.write("##source TranscriptsFilteredByCoverage (coverage<threshold)\n")
        for tid, lines in transcript_lines.items():
            pc = transcript_peptide_count.get(tid, 0)
            cov = tid2_cov.get(tid, 0.0)
            if tid not in kept_tids:
                for line in lines:
                    out = line
                    if (not line.startswith('#')) and ('\ttranscript\t' in line):
                        out = line.rstrip() + f' peptide_count "{pc}"; coverage "{cov:.4f}";'
                    f.write(out + '\n')
    print(f"[OK] 因覆盖率<阈值被剔除的转录本 GTF: {gtf_unsupported_cov}")

    # 覆盖率达标但未能提取CDS（例如非M且无ATG、或无终止密码子等）
    gtf_removed_no_cds = f"{output_prefix}_transcripts_removed_no_cds.gtf"
    with open(gtf_removed_no_cds, 'w') as f:
        f.write("##source TranscriptsRemovedNoCDS (coverage>=threshold BUT cds_extracted=0)\n")
        for tid, lines in transcript_lines.items():
            if (tid in kept_tids) and (tid not in cds_tids):
                pc = transcript_peptide_count.get(tid, 0)
                cov = tid2_cov.get(tid, 0.0)
                for line in lines:
                    out = line
                    if (not line.startswith('#')) and ('\ttranscript\t' in line):
                        reason = res['tid2_reason'].get(tid, '')
                        out = line.rstrip() + f' peptide_count "{pc}"; coverage "{cov:.4f}"; cds_extracted "0"; note "{reason}";'
                    f.write(out + '\n')
    print(f"[OK] 覆盖率达标但未能提取CDS的转录本 GTF: {gtf_removed_no_cds}")

    # CDS FASTA（仅对成功提取的 kept_tids 输出）
    fasta_cds = f"{output_prefix}_cds.fasta"
    with open(fasta_cds, 'w') as f:
        for tid in sorted(cds_tids):
            seq = res['tid2_cds_seq'].get(tid)
            if not seq:
                continue
            s,e = res['tid2_cds_span'][tid]
            start_is_m = res['tid2_start_is_m'].get(tid, False)
            header = f"{tid} cds:{s}-{e} start_is_M:{int(start_is_m)}"
            write_wrapped_fasta(f, header, seq, width=FASTA_LINE_WRAP)
    print(f"[OK] 提取到的 CDS FASTA: {fasta_cds}")

    # 覆盖率与计数分布（同时标记是否成功提取CDS）
    stats_cov = f"{output_prefix}_transcript_coverage.csv"
    df_cov = pd.DataFrame(
        [{'transcript_id': tid,
          'peptide_count': transcript_peptide_count.get(tid, 0),
          'coverage': tid2_cov.get(tid, 0.0),
          'kept_by_cov': tid in kept_tids,
          'cds_extracted': tid in cds_tids}
         for tid in transcript_lines.keys()]
    )
    df_cov.to_csv(stats_cov, index=False)
    print(f"[OK] 转录本覆盖率/肽段计数: {stats_cov}")

    # 已映射/未映射肽段明细
    map_df = pd.DataFrame(res['mapped_pairs_all'], columns=['peptide_id', 'transcript_id'])
    if not map_df.empty:
        agg = (map_df.groupby('peptide_id')['transcript_id']
                   .agg(lambda x: ';'.join(sorted(set(x))))
                   .reset_index()
                   .rename(columns={'transcript_id': 'matched_transcripts'}))
        agg['matched_count'] = agg['matched_transcripts'].apply(lambda s: 0 if not s else s.count(';') + 1)
        mapped_detail = res['peptide_df'].merge(agg, on='peptide_id', how='inner')
    else:
        mapped_detail = res['peptide_df'].head(0).copy()
        mapped_detail['matched_transcripts'] = []
        mapped_detail['matched_count'] = []
    unmapped_detail = pd.DataFrame(res['unmapped_peptides'])

    mapped_xlsx   = f"{output_prefix}_mapped_peptides.xlsx"
    unmapped_xlsx = f"{output_prefix}_unmapped_peptides.xlsx"
    with pd.ExcelWriter(mapped_xlsx, engine='openpyxl') as w:
        mapped_detail.to_excel(w, index=False, sheet_name='mapped')
    with pd.ExcelWriter(unmapped_xlsx, engine='openpyxl') as w:
        unmapped_detail.to_excel(w, index=False, sheet_name='unmapped')
    print(f"[OK] 映射肽段明细: {mapped_xlsx}")
    print(f"[OK] 未映射肽段明细: {unmapped_xlsx}")

    # CDS 提取备注（未能提取的原因）
    reason_tsv = f"{output_prefix}_cds_extraction_notes.tsv"
    rows = []
    for tid in sorted(transcript_lines.keys()):
        rows.append({
            'transcript_id': tid,
            'reason_or_empty': res['tid2_reason'].get(tid, '')
        })
    pd.DataFrame(rows).to_csv(reason_tsv, sep='\t', index=False)
    print(f"[OK] CDS 提取备注: {reason_tsv}")

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
    generate_outputs(res, output_prefix)

if __name__ == "__main__":
    main()
