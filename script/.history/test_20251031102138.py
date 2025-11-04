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

# 外显子长度统计
import pandas as pd
from collections import defaultdict
import re

RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
RE_GENE_ID       = re.compile(r'gene_id "([^"]+)"')

def parse_gtf_exons(gtf_file):
    exon_lengths = []
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            if feature == 'exon':
                start = int(start)
                end = int(end)
                exon_length = end - start + 1
                exon_lengths.append(exon_length)
                if exon_length == 1:
                    print(f"{attrs}")
    return exon_lengths

def plot_exon_length_distribution(exon_lengths, output_file):
    length_distribution = defaultdict(int)
    for length in exon_lengths:
        length_distribution[length] += 1
    df = pd.DataFrame(length_distribution.items(), columns=['Exon Length', 'Frequency'])
    df = df.sort_values(by='Exon Length')  # 按长度排序
    df.to_csv(output_file, index=False)
    print(f"Exon length distribution saved to {output_file}")

def main():
    gtf_file = '/Volumes/caca/Eu_peptido/20251018imeta/new_file_no_predict/01new_gene/analysis_results/filter_by_ms_transcripts_with_peptides.gtf'  # 替换为你的 GTF 文件路径
    output_file = '/Volumes/caca/Eu_peptido/20251018imeta/new_file_no_predict/01new_gene/analysis_results/mapped_transcript_exon_length_distribution.csv'  # 输出的 CSV 文件路径
    exon_lengths = parse_gtf_exons(gtf_file)
    plot_exon_length_distribution(exon_lengths, output_file)

if __name__ == "__main__":
    main()
