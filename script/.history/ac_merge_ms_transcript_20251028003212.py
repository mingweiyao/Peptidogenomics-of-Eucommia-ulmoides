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
        for (exon_start, exon_end) in exon_regions:
            if pstart >= exon_start and pend <= exon_end:
                hit_ids.append(tid)     
                break                  
    return hit_ids

def write_unmapped_peptides_gtf(unmapped_peptides, output_prefix):
    out = f"{output_prefix}_unmapped_peptides.gtf"
    with open(out, 'w') as f:
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("##source UnmappedPeptides\n")
        for i, rec in enumerate(unmapped_peptides, 1):
            peptide_id = str(rec.get('peptide_id'))
            chrom  = str(rec.get('chrom'))
            start  = int(rec.get('start'))
            end    = int(rec.get('end'))
            strand = str(rec.get('strand', '.')) if str(rec.get('strand', '.')) in ['+','-','.'] else '.'
            attrs = [f'transcript_id "{peptide_id}";', f'gene_id "{peptide_id}";']
            attr_str = ' '.join(attrs)
            transcript_line = (
                f"{chrom}\tUnmappedPeptides\ttranscript\t"
                f"{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n"
            )
            f.write(transcript_line)
            exon_line = (
                f"{chrom}\tUnmappedPeptides\texon\t"
                f"{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n"
            )
            f.write(exon_line)
    print(f"未映射肽段的GTF已保存: {out}")

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
    df['start'] = df['start'].astype(int)
    df['end']   = df['end'].astype(int)
    peptides = list(df[['peptide_id', 'chrom', 'start', 'end', 'strand']].itertuples(index=False, name=None))
    counter = Counter()
    unmapped_peptides = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(candidate_index,)) as ex:
        for rec, hit_list in tqdm(ex.map(_worker_match_one_peptide, peptides), total=len(peptides), desc="Matching peptides (parallel)"):
            if hit_list:
                counter.update(hit_list)
            else:
                unmapped_peptides.append(rec)
    transcript_peptide_count = {tid: counter.get(tid, 0) for tid in transcript_lines.keys()}
    return transcript_lines, transcript_peptide_count, unmapped_peptides

def generate_outputs(transcript_lines, transcript_peptide_count, unmapped_peptides, output_prefix):
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

    if unmapped_peptides:
        write_unmapped_peptides_gtf(unmapped_peptides, output_prefix)
    else:
        print("所有肽段都有命中（未生成 unmapped GTF）。")    

def main():
    ms_file = "/media/wanglab/caca/Eu_peptido/20251018imeta/file_no_transdecoder/00raw/Eu_sp_finally.xlsx"
    gtf_file = "/media/wanglab/caca/Eu_peptido/20251018imeta/file_no_transdecoder/00raw/output/file1_nonredundant_coding_transcripts.gtf"
    output_prefix = "/media/wanglab/caca/Eu_peptido/20251018imeta/file_no_transdecoder/01new_gene/analysis_results/filter_by_ms"
    transcript_lines, transcript_peptide_count, unmapped_peptides = filter_by_ms(
        ms_file, gtf_file, workers=70
    )
    generate_outputs(transcript_lines, transcript_peptide_count, unmapped_peptides, output_prefix)

if __name__ == "__main__":
    main()