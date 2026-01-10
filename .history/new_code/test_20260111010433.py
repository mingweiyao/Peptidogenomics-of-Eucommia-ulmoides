import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio.Seq import Seq

_genome_dict = None
_max_scan_nt = None

# -----------------------------
# 1) 常量与全局变量
# -----------------------------
START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
STOP_CODONS  = {"TAA", "TAG", "TGA"}
# 负链：沿用你原脚本的“在正向基因组序列上匹配的三联体集合”
MINUS_START_CODONS = {"CAT", "CAG", "CAC", "CAA", "CGT"}
MINUS_STOP_CODONS  = {"TTA", "CTA", "TCA"}

# 1) 汇总 accession 覆盖信息
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
        if accession not in database_dict: continue
        seq = database_dict[accession]
        length_aa = len(seq)
        min_start, max_end = merge_segments(segments)
        cov_length_aa = (max_end - min_start + 1) / 3.0
        coverage = cov_length_aa / length_aa if length_aa > 0 else None
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
# 2) 多进程扫描并枚举候选
def init_worker(genome_file, max_scan_nt):
    global _genome_dict, _max_scan_nt
    _max_scan_nt = max_scan_nt
    _genome_dict = {}
    for rec in SeqIO.parse(genome_file, "fasta"):
        _genome_dict[rec.id] = rec.seq
def _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        s = min_start - 3 * i
        if s - 1 < 0 or s + 2 > L: continue
        triplet = seq_str[s - 1:s + 2]
        if triplet in START_CODONS:
            yield s, triplet
def _find_stop_plus(seq_str, max_end, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for j in range(max_steps):
        e = max_end + 3 * j
        if e + 3 > L: break
        triplet = seq_str[e:e + 3]
        if triplet in STOP_CODONS: return e
    return None
def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank=6):
    genome_seq = Seq(str(genome_seq).upper())
    if strand == '+': ctx = genome_seq[phy_start - 1 - flank: phy_start + flank + 2]
    else: ctx = genome_seq[phy_end - 3 - flank: phy_end + flank].reverse_complement()
    return {"context": str(ctx)}
def enumerate_orf_candidates_plus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
    seq_str = str(Seq(str(genome_seq)).upper())
    candidates = []
    for phy_start, triplet in _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
        phy_end = _find_stop_plus(seq_str, max_end, max_scan_nt)
        if phy_end is None: continue
        if not (phy_start <= min_start and phy_end >= max_end): continue
        kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='+', flank=flank)
        candidates.append({
            "phy_start": phy_start,
            "phy_end": phy_end,
            "prior": triplet,
            "kozak_seq": kozak["context"],
            "start_to_peptide_nt": int(min_start - phy_start)
        })
    return candidates
def _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        left = max_end + 3 * (i - 1)
        right = max_end + 3 * i
        if right <= 0 or left >= L: continue
        triplet = seq_str[left:right]
        if len(triplet) != 3: continue
        if triplet in MINUS_START_CODONS:
            phy_end = max_end + 3 * i
            yield phy_end, triplet
def _find_stop_minus(seq_str, min_start, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        left = (min_start - 1) - 3 * (i + 1)
        right = (min_start - 1) - 3 * i
        if right <= 0 or left >= L: continue
        triplet = seq_str[left:right]
        if len(triplet) != 3: continue
        phy_start = min_start - 3 * i
        if triplet in MINUS_STOP_CODONS: return phy_start
    return None
def enumerate_orf_candidates_minus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
    seq_str = str(Seq(str(genome_seq)).upper())
    candidates = []
    for phy_end, triplet_raw in _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
        phy_start = _find_stop_minus(seq_str, min_start, max_scan_nt)
        if phy_start is None: continue
        if not (phy_start <= min_start and phy_end >= max_end): continue
        kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='-', flank=flank)
        triplet_rc = str(Seq(triplet_raw).reverse_complement())
        candidates.append({
            "phy_start": phy_start,
            "phy_end": phy_end,
            "prior": triplet_rc,
            "kozak_seq": kozak["context"],
            "start_to_peptide_nt": int(phy_end - max_end)
        })
    return candidates
def run_scan_and_output_for_item(item):
    global _genome_dict, _max_scan_nt
    chrom = item['chrom']
    strand = item['strand']
    min_start = int(item['min_start'])
    max_end = int(item['max_end'])
    if chrom not in _genome_dict:
        item['note'] = 'chrom_not_found'
        item['candidates'] = []
        return item
    gseq = _genome_dict[chrom]
    if strand == '+':
        candidates = enumerate_orf_candidates_plus(gseq, min_start, max_end, _max_scan_nt, flank=6)
    else:
        candidates = enumerate_orf_candidates_minus(gseq, min_start, max_end, _max_scan_nt, flank=6)
    item['candidates'] = candidates
    return item
def run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=20):
    with Pool(processes=nproc, initializer=init_worker, initargs=(genome_file, max_scan_nt)) as pool:
        stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
    return stats_update
def main():
    peptide_file  = "finally_expressed_sp_info.xlsx"
    database_file = "Eu_peptide_database_customized_5.fa"
    genome_file   = "Eu_genome.fasta"
    out_cand = "output_candidates.csv"
    # 1) 汇总 accession 覆盖信息
    stats = filter_peptide_seq_cal_cov(peptide_file, database_file)
    # 2) 多进程扫描并枚举候选
    stats_update = run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=100)
    # 3) 输出 candidates 表
    cand_rows = []
    for it in stats_update:
        base = {k: it[k] for k in it.keys() if k not in ("candidates")}
        for rank, cand in enumerate(it.get("candidates", []), start=1):
            cand_rows.append({**base, "rank": rank, **cand})
    df_cand = pd.DataFrame(cand_rows)
    df_cand.to_csv(out_cand, index=False)
if __name__ == "__main__":
    main()