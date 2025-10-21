import os
import re
from collections import defaultdict
import gffutils
from Bio import SeqIO
import pandas as pd

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
        self.file_deduplication_stats = defaultdict(dict)  # 单文件去重统计
        self.coding_potential_data = {}  # 编码潜力预测数据
    
    def load_coding_potential_predictions(self, cpc2_dir, plek_dir):
        """加载CPC2和PLEK编码潜力预测结果"""
        print("加载编码潜力预测数据...")
        coding_transcripts = set()        
        # 加载CPC2预测结果
        if cpc2_dir and os.path.exists(cpc2_dir):
            print("加载CPC2预测结果...")
            for f in os.listdir(cpc2_dir):
                if f.endswith('.txt') or f.endswith('.tsv') or f.endswith('.csv'):
                    file_path = os.path.join(cpc2_dir, f)
                    try:
                        # CPC2结果格式：ID\tlength\tpeptide\tcoding_probability\tlabel
                        df = pd.read_csv(file_path, sep='\t', header=0)
                        for _, row in df.iterrows():
                            transcript_id = row.iloc[0]
                            # 检查是否为编码（coding）
                            if len(row) >= 5 and str(row.iloc[4]).lower() in ['coding', '1']:
                                coding_transcripts.add(transcript_id)
                    except Exception as e:
                        print(f"  解析CPC2文件 {f} 时出错: {e}")
        # 加载PLEK预测结果
        if plek_dir and os.path.exists(plek_dir):
            print("加载PLEK预测结果...")
            for f in os.listdir(plek_dir):
                if f.endswith('.txt') or f.endswith('.tsv') or f.endswith('.csv'):
                    file_path = os.path.join(plek_dir, f)
                    try:
                        # PLEK结果格式可能不同，根据实际格式调整
                        df = pd.read_csv(file_path, sep='\t', header=0)
                        for _, row in df.iterrows():
                            transcript_id = row.iloc[0]  # 第一列为转录本ID
                            # 检查是否为编码（coding）
                            if len(row) >= 2 and str(row.iloc[1]).lower() in ['coding', '1', 'yes']:
                                coding_transcripts.add(transcript_id)
                    except Exception as e:
                        print(f"  解析PLEK文件 {f} 时出错: {e}")
        self.coding_potential_data = coding_transcripts
        print(f"  总共找到 {len(coding_transcripts)} 个具有编码潜力的转录本")
        return coding_transcripts
    
    def parse_gff3_file_with_coding_filter(self, gff3_file, transcript_to_chrom, source_name):
        """解析单个GFF3文件，先根据编码潜力筛选，再进行文件内部去重"""
        transcripts = []
        file_fingerprint_dict = {}  # 单文件内的指纹字典
        
        try:
            db = gffutils.create_db(gff3_file, dbfn=":memory:", force=True, 
                                  keep_order=True, merge_strategy='merge',
                                  sort_attribute_values=True)
            
            # 统计原始转录本数量
            original_transcripts = list(db.features_of_type('mRNA'))
            original_count = len(original_transcripts)
            
            # 第一步：根据编码潜力筛选
            coding_transcripts_in_file = []
            non_coding_count = 0
            
            for mrna in original_transcripts:
                transcript_id = mrna.id
                
                # 检查是否为编码转录本
                if transcript_id in self.coding_potential_data:
                    coding_transcripts_in_file.append(mrna)
                else:
                    non_coding_count += 1
            
            print(f"  编码潜力筛选: {original_count} → {len(coding_transcripts_in_file)} (移除了 {non_coding_count} 个非编码转录本)")
            
            # 第二步：对编码转录本进行结构解析和去重
            for mrna in coding_transcripts_in_file:
                transcript_id = mrna.id
                seq_id = mrna.seqid
                chrom = transcript_to_chrom.get(seq_id)
                
                if not chrom:
                    # 如果直接匹配失败，尝试其他ID格式
                    for key in transcript_to_chrom.keys():
                        if seq_id in key or key in seq_id:
                            chrom = transcript_to_chrom[key]
                            break
                
                if not chrom:
                    continue  # 跳过没有染色体映射的转录本
                
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
                
                # 单文件内去重
                fingerprint = transcript_obj.get_structure_fingerprint()
                if fingerprint not in file_fingerprint_dict:
                    file_fingerprint_dict[fingerprint] = transcript_obj
                    transcripts.append(transcript_obj)
                else:
                    # 记录重复信息
                    existing_transcript = file_fingerprint_dict[fingerprint]
                    print(f"  文件内重复: {transcript_id} 与 {existing_transcript.transcript_id}")
            
            # 记录单文件去重统计
            coding_count = len(coding_transcripts_in_file)
            unique_count = len(transcripts)
            duplicates_removed = coding_count - unique_count
            
            self.file_deduplication_stats[source_name] = {
                'original': original_count,
                'coding': coding_count,
                'unique': unique_count,
                'non_coding_removed': non_coding_count,
                'duplicates_removed': duplicates_removed
            }
            
            # 报告单文件处理结果
            print(f"  文件内去重: {coding_count} → {unique_count} (去除了 {duplicates_removed} 个结构重复)")
                        
        except Exception as e:
            print(f"解析GFF3文件 {gff3_file} 时出错: {e}")
        
        return transcripts
    
    def parse_fasta_headers(self, fasta_file):
        """解析单个FASTA文件头，构建转录本ID到染色体的映射"""
        transcript_to_chrom = {}
        pattern = r'(\S+)\s+loc:([^|]+)'
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
    
    def find_matching_file_pairs(self, gff3_dir, fasta_dir, cpc2_dir, plek_dir):
        """查找GFF3和FASTA文件对，并加载编码潜力预测"""
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
        # 加载编码潜力预测数据
        if cpc2_dir or plek_dir:
            self.load_coding_potential_predictions(cpc2_dir, plek_dir)        
        return file_pairs
    
    def process_file_pairs(self, file_pairs):
        """处理所有文件对 - 先编码筛选，再去重合并"""
        print(f"开始处理 {len(file_pairs)} 对文件...")
        
        total_original = 0
        total_coding = 0
        total_after_dedup = 0
        
        for i, pair in enumerate(file_pairs):
            self.processed_files += 1
            print(f"处理文件对 {self.processed_files}/{len(file_pairs)}: {pair['base_name']}")
            
            # 步骤1: 解析FASTA文件获取染色体映射
            transcript_to_chrom = self.parse_fasta_headers(pair['fasta'])
            
            if not transcript_to_chrom:
                print(f"  警告: 无法从FASTA文件中解析染色体信息: {pair['fasta']}")
                continue
            
            # 步骤2: 解析GFF3文件（先编码筛选，再去重）
            transcripts = self.parse_gff3_file_with_coding_filter(pair['gff3'], transcript_to_chrom, pair['base_name'])
            
            if not transcripts:
                print(f"  警告: 未从GFF3文件中找到有效编码转录本: {pair['gff3']}")
                continue
            
            # 从file_stats获取准确的统计信息
            file_stats = self.file_deduplication_stats.get(pair['base_name'], {})
            original_count = file_stats.get('original', 0)
            coding_count = file_stats.get('coding', 0)
            unique_count = file_stats.get('unique', len(transcripts))
            
            total_original += original_count
            total_coding += coding_count
            total_after_dedup += unique_count
            
            print(f"  处理结果: 原始{original_count} → 编码{coding_count} → 唯一{unique_count}")
            
            # 步骤3: 更新全局指纹字典
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
                    # 检查是否重复添加
                    if pair['base_name'] not in self.fingerprint_dict[fingerprint]['source_files']:
                        self.fingerprint_dict[fingerprint]['count'] += 1
                        self.fingerprint_dict[fingerprint]['source_files'].append(pair['base_name'])
                        self.fingerprint_dict[fingerprint]['all_transcripts'].append(transcript)
                    else:
                        print(f"  警告: 指纹 {fingerprint} 已包含文件 {pair['base_name']}")
                
            self.all_transcripts.extend(transcripts)
        
        # 总体统计
        print(f"\n处理完成。共处理 {self.processed_files} 个文件对")
        print(f"编码潜力筛选: {total_original} → {total_coding} (移除了 {total_original - total_coding} 个非编码转录本)")
        print(f"单文件去重: {total_coding} → {total_after_dedup} (去除了 {total_coding - total_after_dedup} 个文件内重复)")
        print(f"全局去重后: {len(self.all_transcripts)} → {len(self.fingerprint_dict)} 个唯一编码转录本结构")
        
        # 输出单文件去重摘要
        self._print_file_processing_summary()
    
    def _print_file_processing_summary(self):
        """输出单文件处理统计摘要"""
        print("\n=== 单文件处理统计摘要 ===")
        files_with_coding = []
        
        for file_name, stats in self.file_deduplication_stats.items():
            if stats['coding'] > 0:
                files_with_coding.append((file_name, stats))
        
        if files_with_coding:
            print(f"发现 {len(files_with_coding)} 个文件包含编码转录本:")
            for file_name, stats in files_with_coding[:10]:
                print(f"  {file_name}: 原始{stats['original']} → 编码{stats['coding']} → 唯一{stats['unique']}")
            if len(files_with_coding) > 10:
                print(f"  ... 还有 {len(files_with_coding) - 10} 个文件")
        else:
            print("所有文件均无编码转录本")
    
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
        # 5. 生成单文件处理详细报告
        self._generate_processing_report(output_dir)
    
    def _generate_statistics_table(self, output_dir):
        """生成统计表格"""
        output_file = os.path.join(output_dir, "unique_coding_transcripts_statistics.tsv")
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
        output_file = os.path.join(output_dir, "nonredundant_coding_transcripts.gff3")
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
        output_file = os.path.join(output_dir, "coding_transcript_support_details.tsv")        
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
        output_file = os.path.join(output_dir, "coding_transcripts_processing_summary.txt")        
        with open(output_file, 'w') as f:
            f.write("=== 编码转录本处理摘要报告 ===\n\n")
            f.write(f"处理文件对数: {self.processed_files}\n")
            f.write(f"总编码转录本数: {len(self.all_transcripts)}\n")
            f.write(f"唯一编码转录本结构数: {len(self.fingerprint_dict)}\n")
            f.write(f"冗余度: {len(self.all_transcripts) / len(self.fingerprint_dict):.2f}\n\n")            
            # 支持文件数分布
            count_distribution = defaultdict(int)
            for info in self.fingerprint_dict.values():
                count_distribution[info['count']] += 1            
            f.write("支持文件数分布:\n")
            for count in sorted(count_distribution.keys(), reverse=True):
                freq = count_distribution[count]
                percentage = (freq / len(self.fingerprint_dict)) * 100
                f.write(f"  支持 {count:2d} 个文件: {freq:4d} 个唯一编码转录本 ({percentage:.1f}%)\n")            
            # 染色体分布
            chrom_distribution = defaultdict(int)
            for info in self.fingerprint_dict.values():
                chrom = info['representative'].chrom
                chrom_distribution[chrom] += 1            
            f.write(f"\n染色体分布 (共 {len(chrom_distribution)} 条染色体):\n")
            for chrom in sorted(chrom_distribution.keys()):
                f.write(f"  {chrom}: {chrom_distribution[chrom]} 个唯一编码转录本\n")        
        print(f"处理摘要报告已保存至: {output_file}")
    
    def _generate_processing_report(self, output_dir):
        """生成单文件处理详细报告"""
        output_file = os.path.join(output_dir, "file_processing_report.tsv")
        with open(output_file, 'w') as f:
            f.write("File_Name\tOriginal_Count\tCoding_Count\tUnique_Count\tNonCoding_Removed\tDuplicates_Removed\tCoding_Rate\tDeduplication_Rate\n")
            for file_name, stats in sorted(self.file_deduplication_stats.items()):
                coding_rate = (stats['coding'] / stats['original']) * 100 if stats['original'] > 0 else 0
                dedup_rate = (stats['duplicates_removed'] / stats['coding']) * 100 if stats['coding'] > 0 else 0
                f.write(f"{file_name}\t{stats['original']}\t{stats['coding']}\t{stats['unique']}\t{stats['non_coding_removed']}\t{stats['duplicates_removed']}\t{coding_rate:.2f}%\t{dedup_rate:.2f}%\n")
        print(f"单文件处理报告已保存至: {output_file}")


def main():
    # 配置路径
    gff3_directory = r"G:\Eu_peptido\20251018 imeta\file\00raw\GFF3"     # GFF3文件目录
    fasta_directory = r"G:\Eu_peptido\20251018 imeta\file\00raw\fasta"   # FASTA文件目录
    cpc2_directory = r"G:\Eu_peptido\20251018 imeta\file\00raw\CPC2"     # CPC2预测结果目录
    plek_directory = r"G:\Eu_peptido\20251018 imeta\file\00raw\PLEK"     # PLEK预测结果目录
    output_directory = r"G:\Eu_peptido\20251018 imeta\file\00raw\output" # 输出目录
    processor = TranscriptProcessor()
    file_pairs = processor.find_matching_file_pairs(
        gff3_directory, 
        fasta_directory, 
        cpc2_directory, 
        plek_directory
    )
    processor.process_file_pairs(file_pairs)
    processor.generate_output_files(output_directory)

if __name__ == "__main__":
    main()