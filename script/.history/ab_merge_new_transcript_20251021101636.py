import os
import re
from collections import defaultdict
import gffutils
from Bio import SeqIO

class Transcript:
    """转录本对象，存储完整的结构信息"""
    def __init__(self, transcript_id, chrom, strand, exons, cds_regions, utr5, utr3, source_file):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.exons = sorted(exons, key=lambda x: x[0])
        self.cds_regions = sorted(cds_regions, key=lambda x: x[0])
        self.utr5 = sorted(utr5, key=lambda x: x[0]) if utr5 else []
        self.utr3 = sorted(utr3, key=lambda x: x[0]) if utr3 else []
        self.source_file = source_file
    
    def get_structure_fingerprint(self):
        """生成结构指纹：染色体_链方向_外显子边界_CDS边界_UTR边界"""
        fingerprint_parts = [f"{self.chrom}_{self.strand}"]        
        # 添加外显子边界
        exon_boundaries = []
        for exon in self.exons:
            exon_boundaries.extend([str(exon[0]), str(exon[1])])
        fingerprint_parts.append("-".join(exon_boundaries) if exon_boundaries else "None")
        # 添加CDS边界
        cds_boundaries = []
        for cds in self.cds_regions:
            cds_boundaries.extend([str(cds[0]), str(cds[1])])
        fingerprint_parts.append("-".join(cds_boundaries) if cds_boundaries else "None")        
        # 添加5'UTR边界
        utr5_boundaries = []
        for utr in self.utr5:
            utr5_boundaries.extend([str(utr[0]), str(utr[1])])
        fingerprint_parts.append("-".join(utr5_boundaries) if utr5_boundaries else "None")        
        # 添加3'UTR边界
        utr3_boundaries = []
        for utr in self.utr3:
            utr3_boundaries.extend([str(utr[0]), str(utr[1])])
        fingerprint_parts.append("-".join(utr3_boundaries) if utr3_boundaries else "None")        
        return "|".join(fingerprint_parts)
    
    def get_genomic_span(self):
        """获取基因组跨度"""
        if not self.exons:
            return 0
        starts = [exon[0] for exon in self.exons]
        ends = [exon[1] for exon in self.exons]
        return max(ends) - min(starts) + 1


class TranscriptProcessor:
    def __init__(self):
        self.fingerprint_dict = {}     # 指纹 -> 聚合信息
        self.all_transcripts = []      # 所有转录本对象
        self.processed_files = 0       # 已处理文件计数
    
    def parse_fasta_headers(self, fasta_file):
        """解析单个FASTA文件头，构建转录本ID到染色体的映射"""
        transcript_to_chrom = {}
        pattern = r'(\S+)\s+loc:([^|]+)\|(\d+)-(\d+)\|\.\s+exons:(\d+)-(\d+)\s+segs:(\d+)-(\d+)'
        try:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                header = rec.description
                match = re.search(pattern, header)
                if match:
                    transcript_id = match.group(1)
                    chrom = match.group(2)
                    transcript_to_chrom[transcript_id] = chrom
        except Exception as e:
            print(f"解析FASTA文件 {fasta_file} 时出错: {e}")
        print(f"  从FASTA文件中解析出 {len(transcript_to_chrom)} 个转录本-染色体映射")
        return transcript_to_chrom
    
    def find_matching_file_pairs(self, gff3_dir, fasta_dir):
        """查找GFF3和FASTA文件对"""
        file_pairs = []
        # 获取所有GFF3文件
        gff3_files = {}
        for f in os.listdir(gff3_dir):
            if f.endswith('.gff3') or f.endswith('.gff'):
                base_name = os.path.splitext(f)[0].split('.')[0]
                gff3_files[base_name] = os.path.join(gff3_dir, f)
        # 获取所有FASTA文件并匹配
        for f in os.listdir(fasta_dir):
            if f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fna'):
                base_name = os.path.splitext(f)[0].split('.')[0]
                if base_name in gff3_files:
                    file_pairs.append({
                        'gff3': gff3_files[base_name],
                        'fasta': os.path.join(fasta_dir, f),
                        'base_name': base_name
                    })
        return file_pairs
    
    def parse_gff3_file(self, gff3_file, transcript_to_chrom, source_name):
        """解析单个GFF3文件，提取转录本结构信息"""
        transcripts = []
        try:
            # 使用gffutils解析GFF3文件
            db = gffutils.create_db(gff3_file, dbfn=":memory:", force=True, 
                                  keep_order=True, merge_strategy='merge',
                                  sort_attribute_values=True)
            for mrna in db.features_of_type('mRNA'):
                transcript_id = mrna.id
                seq_id = mrna.seqid
                chrom = transcript_to_chrom.get(seq_id)
                if not chrom:
                    # 如果直接匹配失败，尝试其他ID格式
                    for key in transcript_to_chrom.keys():
                        if seq_id in key or key in seq_id:
                            chrom = transcript_to_chrom[key]
                            break
                strand = mrna.strand if mrna.strand else '+'
                # 获取外显子
                exons = []
                for exon in db.children(mrna.id, featuretype='exon'):
                    exons.append((exon.start, exon.end))
                # 获取CDS区域
                cds_regions = []
                for cds in db.children(mrna.id, featuretype='CDS'):
                    cds_regions.append((cds.start, cds.end))
                # 获取UTR区域
                utr5 = []
                utr3 = []
                for utr in db.children(mrna.id, featuretype=['five_prime_UTR', 'three_prime_UTR']):
                    if utr.featuretype == 'five_prime_UTR':
                        utr5.append((utr.start, utr.end))
                    else:
                        utr3.append((utr.start, utr.end))                
                # 创建转录本对象
                transcript_obj = Transcript(
                    transcript_id=transcript_id,
                    chrom=chrom,
                    strand=strand,
                    exons=exons,
                    cds_regions=cds_regions,
                    utr5=utr5,
                    utr3=utr3,
                    source_file=source_name
                )
                transcripts.append(transcript_obj)                
        except Exception as e:
            print(f"解析GFF3文件 {gff3_file} 时出错: {e}")
        return transcripts
    
    def process_file_pairs(self, file_pairs):
        """处理所有文件对"""
        print(f"开始处理 {len(file_pairs)} 对文件...")
        for i, pair in enumerate(file_pairs):
            self.processed_files += 1
            print(f"处理文件对 {self.processed_files}/{len(file_pairs)}: {pair['base_name']}")
            # 步骤1: 解析FASTA文件获取染色体映射
            transcript_to_chrom = self.parse_fasta_headers(pair['fasta'])
            # 步骤2: 解析GFF3文件
            transcripts = self.parse_gff3_file(pair['gff3'], transcript_to_chrom, pair['base_name'])
            print(f"  从 {pair['base_name']} 中找到 {len(transcripts)} 个转录本")
            # 步骤3: 更新指纹字典
            for transcript in transcripts:
                fingerprint = transcript.get_structure_fingerprint()
                if fingerprint not in self.fingerprint_dict:
                    self.fingerprint_dict[fingerprint] = {
                        'count': 1,
                        'representative': transcript,
                        'source_files': [pair['base_name']],
                        'all_transcripts': [transcript]
                    }
                else:
                    self.fingerprint_dict[fingerprint]['count'] += 1
                    self.fingerprint_dict[fingerprint]['source_files'].append(pair['base_name'])
                    self.fingerprint_dict[fingerprint]['all_transcripts'].append(transcript)
            self.all_transcripts.extend(transcripts)
        print(f"\n处理完成。共处理 {self.processed_files} 个文件对")
        print(f"总转录本数: {len(self.all_transcripts)}")
        print(f"唯一结构数: {len(self.fingerprint_dict)}")
    
    def generate_output_files(self, output_dir):
        """生成所有输出文件"""
        os.makedirs(output_dir, exist_ok=True)
        # 1. 生成统计表格
        self._generate_statistics_table(output_dir)
        # 2. 生成非冗余GFF3文件
        self._generate_nonredundant_gff3(output_dir)
        # 3. 生成支持文件详细列表
        self._generate_support_details(output_dir)
        # 4. 生成处理摘要
        self._generate_summary_report(output_dir)
    
    def _generate_statistics_table(self, output_dir):
        """生成统计表格"""
        output_file = os.path.join(output_dir, "unique_transcripts_statistics.tsv")
        with open(output_file, 'w') as f:
            f.write("Fingerprint\tChromosome\tStrand\tSupport_Count\tRepresentative_Transcript\tExon_Count\tCDS_Count\t5UTR_Count\t3UTR_Count\tGenomic_Span\n")
            for fingerprint, info in sorted(self.fingerprint_dict.items(), 
                                          key=lambda x: x[1]['count'], reverse=True):
                transcript = info['representative']
                f.write(f"{fingerprint}\t{transcript.chrom}\t{transcript.strand}\t")
                f.write(f"{info['count']}\t{transcript.transcript_id}\t")
                f.write(f"{len(transcript.exons)}\t{len(transcript.cds_regions)}\t")
                f.write(f"{len(transcript.utr5)}\t{len(transcript.utr3)}\t")
                f.write(f"{transcript.get_genomic_span()}\n")
        print(f"统计表格已保存至: {output_file}")
    
    def _generate_nonredundant_gff3(self, output_dir):
        """生成非冗余GFF3文件"""
        output_file = os.path.join(output_dir, "nonredundant_transcripts.gff3")
        with open(output_file, 'w') as f:
            f.write("##gff-version 3\n")
            for fingerprint, info in sorted(self.fingerprint_dict.items(), 
                                          key=lambda x: x[1]['count'], reverse=True):
                transcript = info['representative']
                # 写入转录本行
                transcript_start = min([e[0] for e in transcript.exons]) if transcript.exons else 1
                transcript_end = max([e[1] for e in transcript.exons]) if transcript.exons else 1
                f.write(f"{transcript.chrom}\t.\ttranscript\t{transcript_start}\t")
                f.write(f"{transcript_end}\t.\t{transcript.strand}\t.\t")
                f.write(f"ID={transcript.transcript_id};Support_Count={info['count']};Source_Files={','.join(info['source_files'][:5])}{'...' if len(info['source_files']) > 5 else ''}\n")
                # 写入外显子
                for i, exon in enumerate(transcript.exons, 1):
                    f.write(f"{transcript.chrom}\t.\texon\t{exon[0]}\t{exon[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tParent={transcript.transcript_id};exon_id={transcript.transcript_id}.exon{i}\n")          
                # 写入CDS区域
                for i, cds in enumerate(transcript.cds_regions, 1):
                    f.write(f"{transcript.chrom}\t.\tCDS\t{cds[0]}\t{cds[1]}\t.\t")
                    f.write(f"{transcript.strand}\t{i}\tParent={transcript.transcript_id}\n")                
                # 写入5'UTR
                for utr in transcript.utr5:
                    f.write(f"{transcript.chrom}\t.\tfive_prime_UTR\t{utr[0]}\t{utr[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tParent={transcript.transcript_id}\n")                
                # 写入3'UTR
                for utr in transcript.utr3:
                    f.write(f"{transcript.chrom}\t.\tthree_prime_UTR\t{utr[0]}\t{utr[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tParent={transcript.transcript_id}\n")        
        print(f"非冗余GFF3文件已保存至: {output_file}")
    
    def _generate_support_details(self, output_dir):
        """生成支持文件详细列表"""
        output_file = os.path.join(output_dir, "transcript_support_details.tsv")        
        with open(output_file, 'w') as f:
            f.write("Fingerprint\tRepresentative_Transcript\tSupport_Count\tsource_files\n")            
            for fingerprint, info in sorted(self.fingerprint_dict.items(), 
                                          key=lambda x: x[1]['count'], reverse=True):
                transcript = info['representative']
                source_files = ",".join(info['source_files'])
                f.write(f"{fingerprint}\t{transcript.transcript_id}\t{info['count']}\t{source_files}\n")
        print(f"支持文件详细列表已保存至: {output_file}")
    
    def _generate_summary_report(self, output_dir):
        """生成处理摘要报告"""
        output_file = os.path.join(output_dir, "processing_summary.txt")        
        with open(output_file, 'w') as f:
            f.write("=== 转录本处理摘要报告 ===\n\n")
            f.write(f"处理文件对数: {self.processed_files}\n")
            f.write(f"总转录本数: {len(self.all_transcripts)}\n")
            f.write(f"唯一结构数: {len(self.fingerprint_dict)}\n")
            f.write(f"冗余度: {len(self.all_transcripts) / len(self.fingerprint_dict):.2f}\n\n")            
            # 支持文件数分布
            count_distribution = defaultdict(int)
            for info in self.fingerprint_dict.values():
                count_distribution[info['count']] += 1            
            f.write("支持文件数分布:\n")
            for count in sorted(count_distribution.keys(), reverse=True):
                freq = count_distribution[count]
                percentage = (freq / len(self.fingerprint_dict)) * 100
                f.write(f"  支持 {count:2d} 个文件: {freq:4d} 个唯一转录本 ({percentage:.1f}%)\n")            
            # 染色体分布
            chrom_distribution = defaultdict(int)
            for info in self.fingerprint_dict.values():
                chrom = info['representative'].chrom
                chrom_distribution[chrom] += 1            
            f.write(f"\n染色体分布 (共 {len(chrom_distribution)} 条染色体):\n")
            for chrom in sorted(chrom_distribution.keys()):
                f.write(f"  {chrom}: {chrom_distribution[chrom]} 个唯一转录本\n")        
        print(f"处理摘要报告已保存至: {output_file}")

def main():
    # 配置路径
    gff3_directory = "/Users/lemon/Desktop/GFF3"     # GFF3文件目录
    fasta_directory = "/Users/lemon/Desktop/fasta"   # FASTA文件目录
    output_directory = "/Users/lemon/Desktop/output"            # 输出目录
    # 创建处理器
    processor = TranscriptProcessor()
    # 查找匹配的文件对
    file_pairs = processor.find_matching_file_pairs(gff3_directory, fasta_directory)
    # 处理所有文件对
    processor.process_file_pairs(file_pairs)
    # 生成输出文件
    processor.generate_output_files(output_directory)
if __name__ == "__main__":
    main()