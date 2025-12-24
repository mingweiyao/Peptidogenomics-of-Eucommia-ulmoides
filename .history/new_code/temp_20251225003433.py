import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
from multiprocessing import Pool
import warnings
warnings.filterwarnings("ignore")
# -----------------------------
# 1) 常量与全局变量
# -----------------------------
START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
STOP_CODONS  = {"TAA", "TAG", "TGA"}
# 负链：沿用你原脚本的“在正向基因组序列上匹配的三联体集合”
MINUS_START_CODONS = {"CAT", "CAG", "CAC", "CAA", "CGT"}
MINUS_STOP_CODONS  = {"TTA", "CTA", "TCA"}
CODON_PRIOR = {
    "ATG": 0.0,
    "CTG": -1.0, "GTG": -1.0, "TTG": -1.0,
    "ACG": -1.5, "AUA": -1.5, "AUU": -1.5, "AUC": -1.5,
    "AAG": -2.0, "AGG": -2.0, "CGU": -2.0, "CGC": -2.0, "CGG": -2.0, "CAG": -2.0
}
MINUS_CODON_PRIOR = {
    "CAT": 0.0,
    "CAG": -1.0, "CAC": -1.0, "CAA": -1.0,
    "CGT": -1.5, "TAT": -1.5, "AAT": -1.5, "GAT": -1.5,
    "CTT": -2.0, "CCT": -2.0, "ACG": -2.0, "GCG": -2.0, "CCG": -2.0, "CTG": -2.0
}
_genome_dict = None
_max_scan_nt = None
def init_worker(genome_file, max_scan_nt):
    global _genome_dict, _max_scan_nt
    _max_scan_nt = max_scan_nt
    _genome_dict = {}
    for rec in SeqIO.parse(genome_file, "fasta"):
        _genome_dict[rec.id] = rec.seq
# -----------------------------
# 2) 读取 peptide 表、按 accession 合并区间、算 coverage
# -----------------------------
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
# -----------------------------
# 3) Kozak 评分
# -----------------------------
def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank=6):
    genome_seq = Seq(str(genome_seq).upper())
    if strand == '+':
        ctx = genome_seq[phy_start - 1 - flank: phy_start + flank + 2]
    else:
        ctx = genome_seq[phy_end - 3 - flank: phy_end + flank].reverse_complement()
    codon_pos = flank
    if len(ctx) < codon_pos + 4:
        return {"core": None, "extended": None, "total": None, "context": None}
    core = 0.0
    if codon_pos - 3 >= 0 and ctx[codon_pos - 3] in ('A', 'G'):
        core += 2.0
    if codon_pos + 3 < len(ctx) and ctx[codon_pos + 3] == 'G':
        core += 2.0
    extended = 0.0
    for rel, ref in {-6: 'G', -5: 'C', -4: 'C', -2: 'C', +4: 'G'}.items():
        pos = codon_pos + rel
        if 0 <= pos < len(ctx) and ctx[pos] == ref:
            extended += 0.5
    return {"core": core, "extended": extended, "total": core + extended, "context": str(ctx)}
# -----------------------------
# 4) 枚举候选 ORF（+）
# -----------------------------
def _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        s = min_start - 3 * i
        if s - 1 < 0 or s + 2 > L:
            continue
        triplet = seq_str[s - 1:s + 2]
        if triplet in START_CODONS:
            yield s, triplet, CODON_PRIOR.get(triplet, -10.0)
def _find_stop_plus(seq_str, start_pos, max_end, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for j in range(max_steps):
        e = max_end + 3 * j
        if e + 3 > L:
            break
        if (e - start_pos) % 3 != 0:
            continue
        triplet = seq_str[e:e + 3]
        if triplet in STOP_CODONS:
            return e
    return None
def enumerate_orf_candidates_plus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
    seq_str = str(Seq(str(genome_seq)).upper())
    candidates = []
    for phy_start, triplet, prior_value in _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
        phy_end = _find_stop_plus(seq_str, phy_start, max_end, max_scan_nt)
        if phy_end is None:
            continue
        if not (phy_start <= min_start and phy_end >= max_end):
            continue
        kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='+', flank=flank)
        if kozak["total"] is None:
            continue
        total_score = prior_value + kozak["total"]
        candidates.append({
            "phy_start": phy_start,
            "phy_end": phy_end,
            "prior": triplet,
            "prior_value": prior_value,
            "kozak_total": kozak["total"],
            "kozak_seq": kozak["context"],
            "total_score": total_score,
            "start_to_peptide_nt": int(min_start - phy_start)
        })
    candidates.sort(key=lambda d: d["total_score"], reverse=True)
    return candidates
# -----------------------------
# 5) 枚举候选 ORF（-）——沿用你原坐标/切片定义
# -----------------------------
def _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        left = max_end + 3 * (i - 1)
        right = max_end + 3 * i
        if right <= 0 or left >= L:
            continue
        triplet = seq_str[left:right]
        if len(triplet) != 3:
            continue
        if triplet in MINUS_START_CODONS:
            phy_end = max_end + 3 * i
            prior_value = MINUS_CODON_PRIOR.get(triplet, -10.0)
            yield phy_end, triplet, prior_value
def _find_stop_minus(seq_str, phy_end, min_start, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        left = (min_start - 1) - 3 * (i + 1)
        right = (min_start - 1) - 3 * i
        if right <= 0 or left >= L:
            continue
        triplet = seq_str[left:right]
        if len(triplet) != 3:
            continue
        phy_start = min_start - 3 * i
        if (phy_end - phy_start) % 3 != 0:
            continue
        if triplet in MINUS_STOP_CODONS:
            return phy_start
    return None
def enumerate_orf_candidates_minus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
    seq_str = str(Seq(str(genome_seq)).upper())
    candidates = []
    for phy_end, triplet_raw, prior_value in _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
        phy_start = _find_stop_minus(seq_str, phy_end, min_start, max_scan_nt)
        if phy_start is None:
            continue
        if not (phy_start <= min_start and phy_end >= max_end):
            continue
        kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='-', flank=flank)
        if kozak["total"] is None:
            continue
        total_score = prior_value + kozak["total"]
        triplet_rc = str(Seq(triplet_raw).reverse_complement())
        candidates.append({
            "phy_start": phy_start,
            "phy_end": phy_end,
            "prior": triplet_rc,
            "prior_value": prior_value,
            "kozak_total": kozak["total"],
            "kozak_seq": kozak["context"],
            "total_score": total_score,
            "start_to_peptide_nt": int(phy_end - max_end)
        })
    candidates.sort(key=lambda d: d["total_score"], reverse=True)
    return candidates
# -----------------------------
# 6) worker：返回 best + candidates
# -----------------------------
def run_scan_and_output_for_item(item):
    global _genome_dict, _max_scan_nt
    chrom = item['chrom']
    strand = item['strand']
    min_start = int(item['min_start'])
    max_end = int(item['max_end'])
    if chrom not in _genome_dict:
        item['note'] = 'chrom_not_found'
        item['best'] = None
        item['candidates'] = []
        return item
    gseq = _genome_dict[chrom]
    if strand == '+':
        candidates = enumerate_orf_candidates_plus(gseq, min_start, max_end, _max_scan_nt, flank=6)
    else:
        candidates = enumerate_orf_candidates_minus(gseq, min_start, max_end, _max_scan_nt, flank=6)
    item['best'] = candidates[0] if candidates else None
    item['candidates'] = candidates
    return item
def run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=20):
    with Pool(processes=nproc, initializer=init_worker, initargs=(genome_file, max_scan_nt)) as pool:
        stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
    return stats_update
# -----------------------------
# 7) main：输出两张表
# -----------------------------
def main():
    peptide_file  = r"D:\Desktop\peptidemicro\00file\01figure\finally_expressed_sp_info.xlsx"
    database_file = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_peptide_database_customized_5.fa"
    genome_file   = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_genome.fasta"
    out_best = r"D:\Desktop\peptidemicro\00file\01figure\figure4\output_best.csv"
    out_cand = r"D:\Desktop\peptidemicro\00file\01figure\figure4\output_candidates.csv"
    # 1) 汇总 accession 覆盖信息
    stats = filter_peptide_seq_cal_cov(peptide_file, database_file)
    # 2) 多进程扫描并枚举候选
    stats_update = run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=100)
    # 3) 拆成 best 表与 candidates 表
    best_rows = []
    cand_rows = []
    for it in stats_update:
        base = {k: it[k] for k in it.keys() if k not in ("best", "candidates")}
        # best：每个 accession 一行
        if it.get("best") is not None:
            best_rows.append({**base, **it["best"]})
        else:
            best_rows.append({
                **base,
                "phy_start": None, "phy_end": None, "prior": None,
                "prior_value": None, "kozak_total": None, "kozak_seq": None,
                "total_score": None, "start_to_peptide_nt": None
            })
        # candidates：每个 accession 多行（rank 从 1 开始）
        for rank, cand in enumerate(it.get("candidates", []), start=1):
            cand_rows.append({**base, "rank": rank, **cand})
    df_best = pd.DataFrame(best_rows)
    df_cand = pd.DataFrame(cand_rows)
    df_best.to_csv(out_best, index=False)
    df_cand.to_csv(out_cand, index=False)
    print(f"✅ 已输出 best：{out_best}（每个 accession 最优 ORF）")
    print(f"✅ 已输出 candidates：{out_cand}（每个 accession 所有候选 ORF，含 rank）")
if __name__ == "__main__":
    main()
