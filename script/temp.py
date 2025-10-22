import os
import re
from collections import defaultdict
import gffutils
from Bio import SeqIO
import pandas as pd

class Transcript:
    def __init__(self, transcript_id, transcript_seqid, chrom, strand, exons, cds_regions, utr5, utr3, source_file):
        self.transcript_id = transcript_id
        self.transcript_seqid = transcript_seqid
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
        exon_boundaries = []
        for exon in self.exons:
            exon_boundaries.extend([str(exon[0]), str(exon[1])])
        fingerprint_parts.append("-".join(exon_boundaries) if exon_boundaries else "None")
        cds_boundaries = []
        for cds in self.cds_regions:
            cds_boundaries.extend([str(cds[0]), str(cds[1])])
        fingerprint_parts.append("-".join(cds_boundaries) if cds_boundaries else "None")        
        utr5_boundaries = []
        for utr in self.utr5:
            utr5_boundaries.extend([str(utr[0]), str(utr[1])])
        fingerprint_parts.append("-".join(utr5_boundaries) if utr5_boundaries else "None")        
        utr3_boundaries = []
        for utr in self.utr3:
            utr3_boundaries.extend([str(utr[0]), str(utr[1])])
        fingerprint_parts.append("-".join(utr3_boundaries) if utr3_boundaries else "None")        
        return "|".join(fingerprint_parts)

class TranscriptProcessor:
    def __init__(self):
        self.fingerprint_dict = {}
        self.all_transcripts = []
        self.processed_files = 0
        self.file_duplication_stats = defaultdict(dict)
        self.coding_potential_data = {}

    def parse_fasta_header(self, fasta_file):
        transcript_to_chrom = {}
        pattern = r'(\S+)\s+loc:([^|]+)'
        for rec in SeqIO.parse(fasta_file, "fasta"):
            header = rec.description
            match = re.search(pattern, header)
            transcript_id = match.group(1)
            chrom = match.group(2)
            transcript_to_chrom[transcript_id] = chrom
        print(f"  从FASTA文件中解析出 {len(transcript_to_chrom)} 个转录本-染色体映射")
        return transcript_to_chrom

    def find_matching_file_pairs(self, gff3_dir, fasta_dir, cpc2_dir, plek_dir):
        file_pairs = []
        gff3_files = {}
        for f in os.listdir(gff3_dir):
            if f.endswith('.gff3') or f.endswith('.gff'):
                base_name = os.path.splitext(f)[0].split('_')[0]
                gff3_files[base_name] = os.path.join(gff3_dir, f)
        fasta_files = {}
        for f in os.listdir(fasta_dir):
            if f.endswith('.fasta') or f.endswith('.fa') or f.endswith('.fna'):
                base_name = os.path.splitext(f)[0].split('_')[0]
                fasta_files[base_name] = os.path.join(fasta_dir, f)
        cpc2_files = {}
        for f in os.listdir(cpc2_dir):
            if f.endswith('.txt') or f.endswith('tsv'):
                base_name = os.path.splitext(f)[0].split('_')[0]
                cpc2_files[base_name] = os.path.join(cpc2_dir, f)
        for f in os.listdir(plek_dir):
            if f.endswith('.txt') or f.endswith('tsv'):
                base_name = os.path.splitext(f)[0].split('_')[0]
                if (base_name in gff3_files and 
                    base_name in fasta_files and 
                    base_name in cpc2_files):
                    file_pairs.append({
                        "gff3": gff3_files[base_name],
                        "fasta": fasta_files[base_name],
                        "cpc2": cpc2_files[base_name],
                        "plek": os.path.join(plek_dir, f),
                        'base_name': base_name
                    })
        return file_pairs
        
    def loading_coding_potential_predictions(self, cpc2_file, plek_file):
        coding_transcripts = set()
        df_cpc2 = pd.read_csv(cpc2_file, sep='\t', header=0)
        df_cpc2_coding = df_cpc2[df_cpc2.iloc[:, -1] == 'coding']
        coding_transcripts.update(df_cpc2_coding.iloc[:, 0].tolist())
        df_plek = pd.read_csv(plek_file, sep='\t', header=None)
        df_plek_coding = df_plek[df_plek.iloc[:, 0] == 'Coding']
        plek_ids = df_plek_coding.iloc[:, 2].apply(lambda x: str(x).split('>')[-1].split(' ')[0])
        coding_transcripts.update(plek_ids.tolist())
        self.coding_potential_data = coding_transcripts
        return coding_transcripts

    def parse_gff3_file_with_coding_filter(self, gff3_file, transcript_to_chrom, coding_transcripts, source_name):
        transcripts = []
        file_fingerprint_dict = {}
        db = gffutils.create_db(gff3_file, dbfn=":memory:", force=True,     
                                  keep_order=True, merge_strategy='merge',
                                  sort_attribute_values=True)
        original_transcripts = list(db.features_of_type('mRNA'))
        original_count = len(original_transcripts)
        for mrna in original_transcripts:
            if mrna.seqid in coding_transcripts:
                seq_id = mrna.seqid
                chrom = transcript_to_chrom.get(seq_id)
                strand = mrna.strand if mrna.strand else '+'
                exons = []
                for exon in db.children(mrna.id, featuretype='exon'):
                    exons.append((exon.start, exon.end))
                cds_regions = []
                for cds in db.children(mrna.id, featuretype='CDS'):
                    cds_regions.append((cds.start, cds.end))
                utr5 = []
                utr3 = []
                for utr in db.children(mrna.id, featuretype=['five_prime_UTR', 'three_prime_UTR']):
                    if utr.featuretype == 'five_prime_UTR':
                        utr5.append((utr.start, utr.end))
                    else:
                        utr3.append((utr.start, utr.end))    
                transcript_obj = Transcript(
                    transcript_id=mrna.id,
                    transcript_seqid = mrna.seqid,
                    chrom=chrom,
                    strand=strand,
                    exons=exons,
                    cds_regions=cds_regions,
                    utr5=utr5,
                    utr3=utr3,
                    source_file=source_name
                )
                fingerprint = transcript_obj.get_structure_fingerprint()
                if fingerprint not in self.fingerprint_dict:
                    file_fingerprint_dict[fingerprint] = transcript_obj
                    transcripts.append(transcript_obj)
        return transcripts

    def process_file_pairs(self, file_pairs):
        print(f"开始处理 {len(file_pairs)} 对文件...")
        for i, pair in enumerate(file_pairs):
            self.processed_files += 1
            print(f"处理文件对 {self.processed_files}/{len(file_pairs)}: {pair['base_name']}")
            # 步骤1:解析FASTA文件获取染色体映射
            transcript_to_chrom = self.parse_fasta_header(pair['fasta'])
            # 步骤2:加载编码潜力转录本ID
            coding_transcripts = self.loading_coding_potential_predictions(pair['cpc2'], pair['plek'])
            # 步骤3:解析GFF3文件并过滤 编码转录本
            transcripts = self.parse_gff3_file_with_coding_filter(pair['gff3'], transcript_to_chrom, coding_transcripts, pair['base_name'])
            # 步骤4:更新全局指纹字典
            unique_transcript = []
            for transcript in transcripts:
                fingerprint = transcript.get_structure_fingerprint()
                if fingerprint not in self.fingerprint_dict:
                    self.fingerprint_dict[fingerprint] = {
                        'count': 1,
                        'source_files': [pair['base_name']],
                        'all_transcripts': [transcript],
                        "representative": transcript
                    }
                    unique_transcript.append(transcript)
                else:
                    self.fingerprint_dict[fingerprint]['count'] += 1
                    self.fingerprint_dict[fingerprint]['source_files'].append(pair['base_name'])
                    self.fingerprint_dict[fingerprint]['all_transcripts'].append(transcript)
            self.all_transcripts.extend(unique_transcript)
        print(f"\n处理完成。共处理 {self.processed_files} 个文件对")
    
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
        output_file = os.path.join(output_dir, "file1_unique_coding_transcripts_statistics.tsv")
        with open(output_file, 'w') as f:
            f.write("Fingerprint\tChromosome\tStrand\tSupport_Count\tExon_Count\tCDS_Count\t5UTR_Count\t3UTR_Count\n")
            for fingerprint, info in sorted(self.fingerprint_dict.items(), 
                                          key=lambda x: x[1]['count'], reverse=True):
                transcript = info['representative']
                f.write(f"{fingerprint}\t{transcript.chrom}\t{transcript.strand}\t")
                f.write(f"{info['count']}\t{transcript.transcript_id}\t")
                f.write(f"{len(transcript.exons)}\t{len(transcript.cds_regions)}\t")
                f.write(f"{len(transcript.utr5)}\t{len(transcript.utr3)}\n")
        print(f"统计表格已保存至: {output_file}")
    
    def _generate_nonredundant_gff3(self, output_dir):
        """生成非冗余GFF3文件"""
        output_file = os.path.join(output_dir, "file2_nonredundant_coding_transcripts.gff3")
        with open(output_file, 'w') as f:
            f.write("##gff-version 3\n")
            for fingerprint, info in sorted(self.fingerprint_dict.items(), 
                                          key=lambda x: x[1]['count'], reverse=True):
                transcript = info['representative']
                # 写入基因行
                transcript_start = min([e[0] for e in transcript.exons]) if transcript.exons else 1
                transcript_end = max([e[1] for e in transcript.exons]) if transcript.exons else 1
                f.write(f"{transcript.chrom}\t.\tgene\t{transcript_start}\t")
                f.write(f"{transcript_end}\t.\t{transcript.strand}\t.\t")
                f.write(f"ID={transcript.transcript_seqid};Name={fingerprint};Support_Count={info['count']}\n")
                # 写入mRNA
                f.write(f"{transcript.chrom}\t.\tmRNA\t{transcript_start}\t")
                f.write(f"{transcript_end}\t.\t{transcript.strand}\t.\t")
                f.write(f"ID={transcript.transcript_id};Parent={transcript.transcript_seqid};Support_Count={info['count']}\n")              
                # 写入5'UTR
                for utr in transcript.utr5:
                    f.write(f"{transcript.chrom}\t.\tfive_prime_UTR\t{utr[0]}\t{utr[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tParent={transcript.transcript_seqid}\n")                
                # 写入外显子
                for i, exon in enumerate(transcript.exons, 1):
                    f.write(f"{transcript.chrom}\t.\texon\t{exon[0]}\t{exon[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tParent={transcript.transcript_seqid};exon_id={transcript.transcript_id}.exon{i}\n")          
                # 写入CDS区域
                for i, cds in enumerate(transcript.cds_regions, 1):
                    f.write(f"{transcript.chrom}\t.\tCDS\t{cds[0]}\t{cds[1]}\t.\t")
                    f.write(f"{transcript.strand}\t{i}\tParent={transcript.transcript_seqid}\n")  
                # 写入3'UTR
                for utr in transcript.utr3:
                    f.write(f"{transcript.chrom}\t.\tthree_prime_UTR\t{utr[0]}\t{utr[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tParent={transcript.transcript_seqid}\n")        
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
    gff3_directory = "G:/Eu_peptido/20251018 imeta/file/00raw/GFF3"     # GFF3文件目录
    fasta_directory = "G:/Eu_peptido/20251018 imeta/file/00raw/fasta"   # FASTA文件目录
    cpc2_directory = "G:/Eu_peptido/20251018 imeta/file/00raw/cpc2"     # CPC2预测结果目录
    plek_directory = "G:/Eu_peptido/20251018 imeta/file/00raw/plek"     # PLEK预测结果目录
    output_directory = "G:/Eu_peptido/20251018 imeta/file/00raw/output" # 输出目录
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