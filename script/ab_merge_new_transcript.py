import os
import re
from dataclasses import dataclass
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---------- 可选：尽量用 polars，极快；没有则回退到 pandas ----------
try:
    import polars as pl
    _HAS_POLARS = True
except Exception:
    _HAS_POLARS = False
import pandas as pd

# 预编译正则（快）
RE_TRANSCRIPT_ID = re.compile(r'(?:^|;)\s*transcript_id\s+"([^"]+)"')
RE_GENE_ID       = re.compile(r'(?:^|;)\s*gene_id\s+"([^"]+)"')
RE_REF_ID        = re.compile(r'(?:^|;)\s*reference_id\s+"([^"]+)"')
RE_REF_GENE_ID   = re.compile(r'(?:^|;)\s*ref_gene_id\s+"([^"]+)"')
RE_FASTA_ID      = re.compile(r'>(\S+)')

@dataclass
class Transcript:
    gene_id: str
    transcript_id: str
    reference_id: str
    ref_gene_id: str
    chrom: str
    strand: str
    exons: list  # [(start, end), ...] 已排序
    source_file: str
    transcript_len: int = 0
    def __post_init__(self):
        if self.exons:
            self.transcript_len = sum(e - s + 1 for s, e in self.exons)
        else:
            self.transcript_len = 0
    def structure_key(self):
        return (self.chrom, self.strand, tuple(self.exons))

def load_coding_ids(cpc2_file: str, plek_file: str):
    coding = set()
    try:
        df_cpc2 = pd.read_csv(cpc2_file, sep='\t', header=0, dtype='string', engine='c')
    except Exception:
        df_cpc2 = pd.read_csv(cpc2_file, sep='\t', header=0, dtype='string')
    coding.update(
        df_cpc2.loc[df_cpc2.iloc[:, -1].eq('coding'), df_cpc2.columns[0]]
        .dropna().tolist()
    )
    try:
        df_plek = pd.read_csv(plek_file, sep='\t', header=None, usecols=[0, 2], dtype='string', engine='c')
    except Exception:
        df_plek = pd.read_csv(plek_file, sep='\t', header=None, usecols=[0, 2], dtype='string')
    ids = df_plek.loc[df_plek[0].eq('Coding'), 2].str.extract(RE_FASTA_ID, expand=False)
    coding.update(ids.dropna().tolist())
    return coding

def parse_gtf_fast(gtf_file: str, coding_set: set, source_name: str):
    if _HAS_POLARS:
        return _parse_gtf_polars(gtf_file, coding_set, source_name)
    else:
        return _parse_gtf_pandas(gtf_file, coding_set, source_name)

def _parse_gtf_polars(gtf_file: str, coding_set: set, source_name: str):
    df = pl.read_csv(
        gtf_file,
        separator="\t",
        comment_char="#",
        has_header=False,
        new_columns=['chrom','source','feature','start','end','score','strand','frame','attrs'],
        infer_schema_length=0,
        ignore_errors=True
    )
    exons = df.filter(pl.col('feature') == 'exon')
    trans = df.filter(pl.col('feature') == 'transcript')
    def rex(pat): return pl.col('attrs').str.extract(pat, 1)
    exons = exons.with_columns([
        rex(RE_TRANSCRIPT_ID).alias('transcript_id'),
        rex(RE_GENE_ID).alias('gene_id'),
    ])
    trans = trans.with_columns([
        rex(RE_TRANSCRIPT_ID).alias('transcript_id'),
        rex(RE_GENE_ID).alias('gene_id'),
        rex(RE_REF_ID).alias('reference_id'),
        rex(RE_REF_GENE_ID).alias('ref_gene_id'),
    ])
    exons = exons.filter(pl.col('transcript_id').is_in(list(coding_set)))
    trans = trans.filter(pl.col('transcript_id').is_in(list(coding_set)))
    trans_info = (
        trans
        .unique(subset=['transcript_id'], keep='first')
        .select(['transcript_id','gene_id','reference_id','ref_gene_id','chrom','strand'])
        .with_columns([
            pl.col('strand').fill_null('.')
        ])
        .to_dict(as_series=False)
    )
    tdict = {tid: (gene, refid, refgene, chrom, strand)
             for tid, gene, refid, refgene, chrom, strand in zip(
                 trans_info['transcript_id'], trans_info['gene_id'],
                 trans_info['reference_id'], trans_info['ref_gene_id'],
                 trans_info['chrom'], trans_info['strand']
             )}
    exons2 = (
        exons
        .with_columns([
            pl.struct(['start','end']).alias('se')
        ])
        .group_by('transcript_id')
        .agg([
            pl.col('chrom').first().alias('chrom'),
            pl.col('strand').first().alias('strand'),
            pl.col('se').sort_by('start').alias('se_list'),
        ])
    )
    transcripts = []
    for row in exons2.iter_rows(named=True):
        tid = row['transcript_id']
        chrom = row['chrom']
        strand = row['strand'] if row['strand'] is not None else '.'
        ex_list = [(se['start'], se['end']) for se in row['se_list']]

        gene, refid, refgene, tchrom, tstrand = (None, None, None, chrom, strand)
        if tid in tdict:
            gene, refid, refgene, tchrom, tstrand = tdict[tid]
            if tchrom is not None:
                chrom = tchrom
            if tstrand is not None:
                strand = tstrand

        transcripts.append(Transcript(
            gene_id=gene, transcript_id=tid, reference_id=refid, ref_gene_id=refgene,
            chrom=chrom, strand=strand, exons=ex_list, source_file=source_name
        ))
    return transcripts

def _parse_gtf_pandas(gtf_file: str, coding_set: set, source_name: str):
    try:
        df = pd.read_csv(
            gtf_file, sep='\t', comment='#', header=None,
            names=['chrom','source','feature','start','end','score','strand','frame','attrs'],
            dtype={'chrom':'string','source':'string','feature':'string','start':'int64','end':'int64','score':'string','strand':'string','frame':'string','attrs':'string'},
            engine='c'
        )
    except Exception:
        df = pd.read_csv(
            gtf_file, sep='\t', comment='#', header=None,
            names=['chrom','source','feature','start','end','score','strand','frame','attrs'],
            dtype='string'
        )
        df['start'] = df['start'].astype('int64')
        df['end'] = df['end'].astype('int64')
    exons = df[df['feature'] == 'exon'].copy()
    trans = df[df['feature'] == 'transcript'].copy()
    def extract(series, pat):
        return series.str.extract(pat, expand=False)
    exons['transcript_id'] = extract(exons['attrs'], RE_TRANSCRIPT_ID)
    exons['gene_id']       = extract(exons['attrs'], RE_GENE_ID)
    trans['transcript_id'] = extract(trans['attrs'], RE_TRANSCRIPT_ID)
    trans['gene_id']       = extract(trans['attrs'], RE_GENE_ID)
    trans['reference_id']  = extract(trans['attrs'], RE_REF_ID)
    trans['ref_gene_id']   = extract(trans['attrs'], RE_REF_GENE_ID)
    exons = exons[exons['transcript_id'].isin(coding_set)]
    trans = trans[trans['transcript_id'].isin(coding_set)]
    trans_info = (trans.drop_duplicates('transcript_id')
                      .set_index('transcript_id')[['gene_id','reference_id','ref_gene_id','chrom','strand']])
    transcripts = []
    for tid, g in exons.groupby('transcript_id', sort=False):
        g_sorted = g.sort_values('start', kind='mergesort')
        ex_list = list(zip(g_sorted['start'].tolist(), g_sorted['end'].tolist()))
        gene = refid = refgene = None
        chrom = g_sorted['chrom'].iloc[0]
        strand = g_sorted['strand'].iloc[0] if pd.notna(g_sorted['strand'].iloc[0]) else '.'
        if tid in trans_info.index:
            trow = trans_info.loc[tid]
            gene, refid, refgene = trow['gene_id'], trow['reference_id'], trow['ref_gene_id']
            chrom = trow['chrom'] if pd.notna(trow['chrom']) else chrom
            strand = trow['strand'] if pd.notna(trow['strand']) else strand
        transcripts.append(Transcript(
            gene_id=gene, transcript_id=tid, reference_id=refid, ref_gene_id=refgene,
            chrom=chrom, strand=strand, exons=ex_list, source_file=source_name
        ))
    return transcripts

def _process_one_pair(pair):
    base = pair['base_name']
    coding = load_coding_ids(pair['cpc2'], pair['plek'])
    transcripts = parse_gtf_fast(pair['gtf'], coding, base)
    local = {}
    for t in transcripts:
        key = t.structure_key()
        if key not in local:
            local[key] = {'count': 1, 'representative': t, 'sources': {base}}
        else:
            rep = local[key]['representative']
            if t.transcript_len == rep.transcript_len:
                local[key]['count'] += 1
            else:
                local[key]['count'] += 1
                if t.transcript_len > rep.transcript_len:
                    local[key]['representative'] = t
    return base, local

class FastTranscriptProcessor:
    def __init__(self):
        self.fingerprint_dict = {}
        self.all_transcripts = []
        self.processed_files = 0
        self.saturation_data = []

    @staticmethod
    def find_matching_file_pairs(gtf_dir, cpc2_dir, plek_dir):
        def index_dir(dirp, exts):
            m = {}
            for f in os.listdir(dirp):
                if any(f.endswith(e) for e in exts):
                    base = os.path.splitext(f)[0].split('_')[0]
                    m[base] = os.path.join(dirp, f)
            return m
        gtf = index_dir(gtf_dir, ['.gtf'])
        cpc2 = index_dir(cpc2_dir, ['.txt', '.tsv'])
        plek = index_dir(plek_dir, ['.txt', '.tsv'])
        pairs = []
        for b in sorted(gtf.keys()):
            if b in cpc2 and b in plek:
                pairs.append({'gtf': gtf[b], 'cpc2': cpc2[b], 'plek': plek[b], 'base_name': b})
        return pairs

    def process_file_pairs(self, file_pairs):
        print(f"开始处理 {len(file_pairs)} 对文件（并行）...")
        results = {}
        with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
            futs = [ex.submit(_process_one_pair, p) for p in file_pairs]
            for fu in as_completed(futs):
                base, local = fu.result()
                results[base] = local
        cumulative_unique = 0
        for i, pair in enumerate(file_pairs):
            base = pair['base_name']
            local = results[base]
            new_unique = 0
            for k, v in local.items():
                if k not in self.fingerprint_dict:
                    self.fingerprint_dict[k] = {
                        'count': v['count'],
                        'representative': v['representative'],
                        'source_files': set(v['sources']),
                    }
                    self.all_transcripts.append(v['representative'])
                    new_unique += 1
                else:
                    g = self.fingerprint_dict[k]
                    g['count'] += v['count']
                    g['source_files'] |= v['sources']
                    cur = g['representative']
                    cand = v['representative']
                    if cand.transcript_len > cur.transcript_len:
                        g['representative'] = cand
            cumulative_unique = len(self.fingerprint_dict)
            self.processed_files += 1
            self.saturation_data.append({
                'file_name': base,
                'files_processed': i + 1,
                'new_unique_transcripts': new_unique,
                'cumulative_unique_transcripts': cumulative_unique,
                'total_transcripts_so_far': len(self.all_transcripts)
            })
            print(f"[{i+1}/{len(file_pairs)}] {base}: 新增 {new_unique}，累计唯一 {cumulative_unique}")
        print(f"\n处理完成，共 {self.processed_files} 个文件对。")

    def generate_outputs(self, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        self._write_nonredundant_gtf(output_dir)
        self._write_saturation(output_dir)
        self._write_summary(output_dir)

    def _write_nonredundant_gtf(self, output_dir):
        out = os.path.join(output_dir, "file1_nonredundant_coding_transcripts.gtf")
        lines = []
        for k, info in sorted(self.fingerprint_dict.items(), key=lambda x: x[1]['count'], reverse=True):
            t = info['representative']
            if t.exons:
                t_start = t.exons[0][0]
                t_end   = t.exons[-1][1]
            else:
                t_start = t_end = 1
            lines.append(
                f"{t.chrom}\t.\ttranscript\t{t_start}\t{t_end}\t.\t{t.strand}\t.\t"
                f'gene_id "{t.gene_id}"; transcript_id "{t.transcript_id}"; reference_id "{t.reference_id}"; ref_gene_id "{t.ref_gene_id}";\n'
            )
            for i, (s, e) in enumerate(t.exons, 1):
                lines.append(
                    f"{t.chrom}\t.\texon\t{s}\t{e}\t.\t{t.strand}\t.\t"
                    f'gene_id "{t.gene_id}"; transcript_id "{t.transcript_id}"; exon_number "{i}"; reference_id "{t.reference_id}"; ref_gene_id "{t.ref_gene_id}";\n'
                )
        with open(out, 'w') as f:
            f.writelines(lines)
        print(f"非冗余GTF已保存: {out}")

    def _write_saturation(self, output_dir):
        out = os.path.join(output_dir, "file2_saturation_curve_data.tsv")
        lines = ["File_Order\tFile_Name\tFiles_Processed\tNew_Unique_Transcripts\tCumulative_Unique_Transcripts\tTotal_Transcripts\n"]
        for i, d in enumerate(self.saturation_data, 1):
            lines.append(f"{i}\t{d['file_name']}\t{d['files_processed']}\t{d['new_unique_transcripts']}\t{d['cumulative_unique_transcripts']}\t{d['total_transcripts_so_far']}\n")
        with open(out, 'w') as f:
            f.writelines(lines)
        print(f"饱和曲线数据已保存: {out}")

    def _write_summary(self, output_dir):
        out = os.path.join(output_dir, "file3_coding_transcripts_processing_summary.txt")
        lines = []
        lines.append("=== 编码转录本处理摘要报告 ===\n\n")
        lines.append(f"处理文件对数: {self.processed_files}\n")
        lines.append(f"总编码转录本数(代表作): {len(self.all_transcripts)}\n")
        lines.append(f"唯一编码转录本结构数: {len(self.fingerprint_dict)}\n")
        cnt = defaultdict(int)
        for info in self.fingerprint_dict.values():
            cnt[len(info['source_files'])] += 1
        lines.append("支持文件数分布:\n")
        total_unique = len(self.fingerprint_dict) if self.fingerprint_dict else 1
        for c in sorted(cnt.keys(), reverse=True):
            freq = cnt[c]
            pct = (freq / total_unique) * 100
            lines.append(f"  支持 {c:2d} 个文件: {freq:4d} 个唯一编码转录本 ({pct:.1f}%)\n")
        chrom_cnt = defaultdict(int)
        for info in self.fingerprint_dict.values():
            chrom_cnt[info['representative'].chrom] += 1
        lines.append(f"\n染色体分布 (共 {len(chrom_cnt)} 条染色体):\n")
        for chrom in sorted(chrom_cnt.keys()):
            lines.append(f"  {chrom}: {chrom_cnt[chrom]} 个唯一编码转录本\n")
        with open(out, 'w') as f:
            f.writelines(lines)
        print(f"处理摘要报告已保存: {out}")

def main():
    gtf_dir = "/media/wanglab/caca/Eu_peptido/20251018imeta/file_no_transdecoder/00raw/predict_gtf"
    cpc2_dir = "/media/wanglab/caca/Eu_peptido/20251018imeta/file_no_transdecoder/00raw/predict_cpc2"
    plek_dir = "/media/wanglab/caca/Eu_peptido/20251018imeta/file_no_transdecoder/00raw/predict_plek"
    out_dir = "/media/wanglab/caca/Eu_peptido/20251018imeta/file_no_transdecoder/00raw/output"
    proc = FastTranscriptProcessor()
    pairs = proc.find_matching_file_pairs(gtf_dir, cpc2_dir, plek_dir)
    print(f"找到 {len(pairs)} 对匹配的文件")
    proc.process_file_pairs(pairs)
    proc.generate_outputs(out_dir)
    print("\n=== 处理完成 ===")
    print(f"输出文件保存在: {out_dir}")
    print("生成的文件:")
    print("  1. file1_nonredundant_coding_transcripts.gtf")
    print("  2. file2_saturation_curve_data.tsv")
    print("  3. file3_coding_transcripts_processing_summary.txt")

if __name__ == "__main__":
    main()
