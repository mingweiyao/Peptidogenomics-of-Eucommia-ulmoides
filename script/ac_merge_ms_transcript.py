import pandas as pd
from Bio import SeqIO
from datetime import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor

# ===== 全局变量（在子进程中只读）=====
_CANDIDATE_INDEX = None  # {(chrom,strand): [(transcript_id, [(cds_start,cds_end),...])]}
def _init_worker(candidate_index):
    # 每个子进程启动时挂载一次，后续任务直接复用，避免反复序列化巨大对象
    global _CANDIDATE_INDEX
    _CANDIDATE_INDEX = candidate_index

def _worker_match_one_peptide(peptide_tuple):
    """
    子进程执行：对一个肽段，返回匹配到的 transcript_id 列表（位于任一CDS且相位对齐）
    peptide_tuple: (chrom, start, end, strand)
    """
    chrom, pstart, pend, strand = peptide_tuple
    hit_ids = []
    # 候选仅限同染色体+同链，减少对比规模
    for tid, cds_regions in _CANDIDATE_INDEX.get((chrom, strand), []):
        for cds_start, cds_end in cds_regions:
            # 与你原逻辑完全一致：肽段全落在某个CDS内，且与CDS起点相位对齐
            if pstart >= cds_start and pend <= cds_end and ((pstart - cds_start) % 3 == 0):
                hit_ids.append(tid)
                break
    return hit_ids

def filter_by_ms(ms_file, gff_file, pep_file, workers=None):
    # ---------- 读 GFF3，收集每个转录本的整段行 + mRNA基本信息 + CDS区间 ----------
    transcript_lines = {}
    transcript_mrna_info = {}
    transcript_cds_regions = {}

    current_gene_lines = []
    current_transcript_id = None

    with open(gff_file, 'r') as f:
        for line in tqdm(f, desc="Loading GFF3"):
            if not line or line.startswith('#'):
                continue
            line = line.rstrip('\n')
            fields = line.split('\t')
            chrom, _, ftype, start, end, _, strand, _, attrs = fields
            if ftype == 'gene':
                if current_transcript_id and current_gene_lines:
                    transcript_lines[current_transcript_id] = current_gene_lines.copy()
                current_gene_lines = [line]
                current_transcript_id = None
            elif ftype == 'mRNA':
                tid = attrs.split(';')[1].split('=')[-1]
                current_transcript_id = tid
                current_gene_lines.append(line)
                transcript_mrna_info[tid] = {
                    'chrom': chrom,
                    'strand': strand,
                    'start': int(start),
                    'end': int(end)
                }
            else:
                if current_gene_lines:
                    current_gene_lines.append(line)
            if current_transcript_id:
                transcript_lines[current_transcript_id] = current_gene_lines
        if current_transcript_id and current_gene_lines:
            transcript_lines[current_transcript_id] = current_gene_lines
            
    for tid, lines in tqdm(transcript_lines.items(), desc="Parsing CDS regions"):
        cds = []
        for line in lines:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if f[2] == 'CDS':
                cds.append((int(f[3]), int(f[4])))
        transcript_cds_regions[tid] = sorted(cds)

    candidate_index = defaultdict(list)
    for tid, info in transcript_mrna_info.items():
        cds_regions = transcript_cds_regions.get(tid, [])
        if not cds_regions:
            continue
        key = (info['chrom'], info['strand'])
        candidate_index[key].append((tid, cds_regions))

    df = pd.read_excel(ms_file, sheet_name="NCP", engine="openpyxl")
    peptides = list(df[['chrom', 'start', 'end', 'strand']].itertuples(index=False, name=None))

    counter = Counter()
    remain_id = set()

    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(candidate_index,)) as ex:
        for hit_list in tqdm(ex.map(_worker_match_one_peptide, peptides), total=len(peptides), desc="Matching peptides (parallel)"):
            if hit_list:
                counter.update(hit_list)
                remain_id.update(hit_list)

    # ---------- 读取 PEP 序列（仅取命中的 ID） ----------
    pep_sequences = {}
    if remain_id:
        for record in SeqIO.parse(pep_file, "fasta"):
            if record.id in remain_id:
                pep_sequences[record.id] = str(record.seq)

    # 将 Counter 转为普通 dict 以保持你原调用习惯
    transcript_peptide_count = {tid: counter.get(tid, 0) for tid in transcript_lines.keys()}

    return transcript_lines, transcript_peptide_count, pep_sequences

def generate_outputs(transcript_lines, transcript_peptide_count, pep_sequences, output_prefix):
    # 1) GFF3（仅输出有肽段支持的转录本）
    gff_output_file = f"{output_prefix}_transcripts_with_peptides.gff3"
    with open(gff_output_file, 'w') as f:
        f.write("##gff-version 3\n")
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("##source PeptideSupportedTranscripts\n")
        for tid, lines in transcript_lines.items():
            pc = transcript_peptide_count.get(tid, 0)
            if pc > 0:
                for line in lines:
                    if (not line.startswith('#')) and ('\tmRNA\t' in line):
                        # 安全追加属性
                        if line.endswith(';') or ';' in line:
                            line = line.rstrip() + f"peptide_count={pc}"
                        else:
                            line = line.rstrip() + f";peptide_count={pc}"
                    f.write(line + '\n')
                f.write("###\n")

    # 2) 支持的转录本蛋白序列
    fasta_output_file = f"{output_prefix}_transcripts_with_peptides.fasta"
    if pep_sequences:
        recs = [SeqRecord(Seq(seq), id=tid, description="") for tid, seq in pep_sequences.items()]
        with open(fasta_output_file, "w") as fh:
            SeqIO.write(recs, fh, "fasta")

    # 3) 统计
    stats_output_file = f"{output_prefix}_peptide_count_distribution.csv"
    cnt = Counter(transcript_peptide_count.values())
    pd.DataFrame(
        [{'peptide_count': k, 'transcript_count': v} for k, v in sorted(cnt.items())]
    ).to_csv(stats_output_file, index=False)

def main():
    ms_file = "/media/wanglab/caca/Eu_peptido/20251018imeta/file/00raw/Eu_sp_finally.xlsx"
    gff_file = "/media/wanglab/caca/Eu_peptido/20251018imeta/file/00raw/output/filter_file1_nonredundant_coding_transcripts.gff3"
    pep_file = "/media/wanglab/caca/Eu_peptido/20251018imeta/file/00raw/output/filter_file2_nonredundant_coding_transcript_pep.fasta"
    output_prefix = "/media/wanglab/caca/Eu_peptido/20251018imeta/file/01new_gene/analysis_results/filter_by_ms"
    transcript_lines, transcript_peptide_count, pep_sequences = filter_by_ms(
        ms_file, gff_file, pep_file, workers=100
    )
    generate_outputs(transcript_lines, transcript_peptide_count, pep_sequences, output_prefix)

if __name__ == "__main__":
    main()