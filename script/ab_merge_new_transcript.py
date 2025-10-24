import os
import re
from collections import defaultdict
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
        self.fingerprint = None
    
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
        self.sequence_dict = {}
        self.saturation_data = [] 

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

    def find_matching_file_pairs(self, gff3_dir, fasta_dir, cpc2_dir, plek_dir, pep_dir):
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
        pep_files = {}
        for f in os.listdir(pep_dir):
            if f.endswith('.pep'):
                base_name = os.path.splitext(f)[0].split('_')[0]
                pep_files[base_name] = os.path.join(pep_dir, f)
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
                        "pep": pep_files[base_name],
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

    def extract_fasta(self, transcripts, pep_file):
        current_file_sequences = {}
        transcript_ids = [t.transcript_id for t in transcripts]
        for rec in SeqIO.parse(pep_file, "fasta"):
            if rec.id in transcript_ids:
                current_file_sequences[rec.id] = str(rec.seq)
        self.sequence_dict.update(current_file_sequences)

    def parse_gff3_file_with_coding_filter(self, gff3_file, transcript_to_chrom, 
                                           coding_transcripts, source_name, pep_file):
        transcripts = []
        file_fingerprint_dict = {}
        db = gffutils.create_db(gff3_file, dbfn=":memory:", force=True,     
                                  keep_order=True, merge_strategy='merge',
                                  sort_attribute_values=True)
        original_transcripts = list(db.features_of_type('mRNA'))
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
                    transcript_id=mrna.id, # have .p1
                    transcript_seqid = mrna.seqid, # do not have .p1
                    chrom=chrom,
                    strand=strand,
                    exons=exons,
                    cds_regions=cds_regions,
                    utr5=utr5,
                    utr3=utr3,
                    source_file=source_name
                )
                fingerprint = transcript_obj.get_structure_fingerprint()
                transcript_obj.fingerprint = fingerprint
                if fingerprint not in file_fingerprint_dict:
                    file_fingerprint_dict[fingerprint] = transcript_obj
                    transcripts.append(transcript_obj)
        self.extract_fasta(transcripts, pep_file)
        return transcripts
    
    def collapse_by_mrna_fingerprint_keep_longest_cds(self):
        def mrna_key_from_fingerprint(fp: str):
            # fingerprint 形如: "chrom_strand|exon_bounds|cds_bounds|utr5|utr3"
            parts = fp.split('|')
            return parts[0]
        def is_redundant(start1, end1, start2, end2):
            return start1 >= start2 and end1 <= end2
        groups = defaultdict(list)
        for fp, info in self.fingerprint_dict.items():
            rep = info["representative"]
            mrna_key = mrna_key_from_fingerprint(fp)
            groups[mrna_key].append((fp, rep))
        keep_fps = set()
        for mrna_key, items in groups.items():
            if len(items) == 1:
                keep_fps.add(items[0][0])
            else:
                items.sort(key=lambda kv: (kv[1].cds_regions[0][0], -kv[1].cds_regions[-1][1]))
                best_fp = None
                for i, (fp1, rep1) in enumerate(items):
                    start1, end1 = rep1.cds_regions[0][0], rep1.cds_regions[-1][1]
                    is_valid = True
                    if best_fp is not None:
                        best_rep = self.fingerprint_dict[best_fp]["representative"]
                        best_start, best_end = best_rep.cds_regions[0][0], best_rep.cds_regions[-1][1]
                        if is_redundant(start1, end1, best_start, best_end):
                            is_valid = False 
                    if is_valid:
                        keep_fps.add(fp1)
                        best_fp = fp1
        self.fingerprint_dict = {fp: info for fp, info in self.fingerprint_dict.items() if fp in keep_fps}
        self.all_transcripts = [info["representative"] for info in self.fingerprint_dict.values()]

    def process_file_pairs(self, file_pairs):
        print(f"开始处理 {len(file_pairs)} 对文件...")
        cumulative_unique = 0
        file_order = []
        for i, pair in enumerate(file_pairs):
            self.processed_files += 1
            print(f"处理文件对 {self.processed_files}/{len(file_pairs)}: {pair['base_name']}")
            # 步骤1:解析FASTA文件获取染色体映射
            transcript_to_chrom = self.parse_fasta_header(pair['fasta'])
            # 步骤2:加载编码潜力转录本ID
            coding_transcripts = self.loading_coding_potential_predictions(pair['cpc2'], pair['plek'])
            # 步骤3:解析GFF3文件并过滤 编码转录本
            transcripts = self.parse_gff3_file_with_coding_filter(pair['gff3'], transcript_to_chrom, coding_transcripts, pair['base_name'], pair['pep'])
            # 步骤4:更新全局指纹字典
            unique_transcript = []
            new_unique_count = 0
            for transcript in transcripts:
                fingerprint = transcript.fingerprint
                if fingerprint not in self.fingerprint_dict:
                    self.fingerprint_dict[fingerprint] = {
                        'count': 1,
                        'source_files': [pair['base_name']],
                        'all_transcripts': [transcript],
                        "representative": transcript
                    }
                    new_unique_count += 1
                    unique_transcript.append(transcript)
                else:
                    self.fingerprint_dict[fingerprint]['count'] += 1
                    self.fingerprint_dict[fingerprint]['source_files'].append(pair['base_name'])
                    self.fingerprint_dict[fingerprint]['all_transcripts'].append(transcript)
            self.extract_fasta(unique_transcript, pair['pep'])
            self.all_transcripts.extend(unique_transcript)
            cumulative_unique = len(self.fingerprint_dict)
            file_order.append(pair['base_name'])
            self.saturation_data.append({
                'file_name': pair['base_name'],
                'files_processed': i + 1,
                'new_unique_transcripts': new_unique_count,
                'cumulative_unique_transcripts': cumulative_unique,
                'total_transcripts_so_far': len(self.all_transcripts)
            })            
        print(f"\n处理完成。共处理 {self.processed_files} 个文件对")
    
    def generate_output_files(self, output_dir):
        """生成所有输出文件"""
        os.makedirs(output_dir, exist_ok=True)
        # 1. 生成去重合并后的GFF3文件
        self._generate_nonredundant_gff3(output_dir)
        # 2. 生成去重合并后的FASTA文件
        self._generate_nonredundant_fasta(output_dir)        
        # 3. 生成饱和曲线数据
        self._generate_saturation_curve_data(output_dir)        
        # 4. 生成处理摘要
        self._generate_summary_report(output_dir)
    
    def _generate_nonredundant_gff3(self, output_dir):
        """生成非冗余GFF3文件"""
        output_file = os.path.join(output_dir, "filter_file1_nonredundant_coding_transcripts.gff3")
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
                f.write(f"ID={transcript.transcript_id};Name={fingerprint}\n")
                # 写入mRNA
                f.write(f"{transcript.chrom}\t.\tmRNA\t{transcript_start}\t")
                f.write(f"{transcript_end}\t.\t{transcript.strand}\t.\t")
                f.write(f"ID={transcript.transcript_id}.1;Parent={transcript.transcript_id}\n")              
                # 写入5'UTR
                for utr in transcript.utr5:
                    f.write(f"{transcript.chrom}\t.\tfive_prime_UTR\t{utr[0]}\t{utr[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tID={transcript.transcript_id}.1.utr5;Parent={transcript.transcript_id}.1\n")                
                # 写入外显子
                for i, exon in enumerate(transcript.exons, 1):
                    f.write(f"{transcript.chrom}\t.\texon\t{exon[0]}\t{exon[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tID={transcript.transcript_id}.1.exon;Parent={transcript.transcript_id}.1\n")          
                # 写入CDS区域
                for i, cds in enumerate(transcript.cds_regions, 1):
                    f.write(f"{transcript.chrom}\t.\tCDS\t{cds[0]}\t{cds[1]}\t.\t")
                    f.write(f"{transcript.strand}\t0\tID=cds.{transcript.transcript_id}.1;Parent={transcript.transcript_id}.1\n")  
                # 写入3'UTR
                for utr in transcript.utr3:
                    f.write(f"{transcript.chrom}\t.\tthree_prime_UTR\t{utr[0]}\t{utr[1]}\t.\t")
                    f.write(f"{transcript.strand}\t.\tID={transcript.transcript_id}.1.utr3;Parent={transcript.transcript_id}.1\n")        
        print(f"非冗余GFF3文件已保存至: {output_file}")

    def _generate_nonredundant_fasta(self, output_dir):
        output_file = os.path.join(output_dir, "filter_file2_nonredundant_coding_transcript_pep.fasta")
        records = []
        for fp, info in self.fingerprint_dict.items():
            rep_id = info["representative"].transcript_id
            seq = self.sequence_dict.get(rep_id)
            if not seq:
                continue  # 代表转录本在 pep 里没找到的话就跳过或记录
            records.append(SeqRecord(Seq(seq), id=rep_id))
        with open(output_file, 'w') as f:
            SeqIO.write(records, f, "fasta")
                
    def _generate_saturation_curve_data(self, output_dir):
        output_file = os.path.join(output_dir, "filter_file3_saturation_curve_data.tsv")
        with open(output_file, 'w') as f:
            f.write("File_Order\tFile_Name\tFiles_Processed\tNew_Unique_Transcripts\tCumulative_Unique_Transcripts\tTotal_Transcripts\n")
            for i, data in enumerate(self.saturation_data, 1):
                f.write(f"{i}\t{data['file_name']}\t{data['files_processed']}\t")
                f.write(f"{data['new_unique_transcripts']}\t{data['cumulative_unique_transcripts']}\t")
                f.write(f"{data['total_transcripts_so_far']}\n")    
        print(f"饱和曲线数据已保存至: {output_file}")        
    
    def _generate_summary_report(self, output_dir):
        """生成处理摘要报告"""
        output_file = os.path.join(output_dir, "filter_file4_coding_transcripts_processing_summary.txt")        
        with open(output_file, 'w') as f:
            f.write("=== 编码转录本处理摘要报告 ===\n\n")
            f.write(f"处理文件对数: {self.processed_files}\n")
            f.write(f"总编码转录本数: {len(self.all_transcripts)}\n")
            f.write(f"唯一编码转录本结构数: {len(self.fingerprint_dict)}\n")         
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

def main():
    # 配置路径
    gff3_directory = r"G:\Eu_peptido\20251018imeta\file\00raw\gff3"     # GFF3文件目录
    fasta_directory = r"G:\Eu_peptido\20251018imeta\file\00raw\fasta"   # FASTA文件目录
    cpc2_directory = r"G:\Eu_peptido\20251018imeta\file\00raw\cpc2"     # CPC2预测结果目录
    plek_directory = r"G:\Eu_peptido\20251018imeta\file\00raw\plek"     # PLEK预测结果目录
    pep_directory = r"G:\Eu_peptido\20251018imeta\file\00raw\pep"
    output_directory = r"G:\Eu_peptido\20251018imeta\file\00raw\output" # 输出目录
    processor = TranscriptProcessor()
    file_pairs = processor.find_matching_file_pairs(
        gff3_directory, 
        fasta_directory, 
        cpc2_directory, 
        plek_directory,
        pep_directory
    )
    print(f"找到 {len(file_pairs)} 对匹配的文件")
    processor.process_file_pairs(file_pairs)
    processor.collapse_by_mrna_fingerprint_keep_longest_cds()
    processor.generate_output_files(output_directory)
    print("\n=== 处理完成 ===")
    print(f"输出文件保存在: {output_directory}")
    print("生成的文件:")
    print("  1. nonredundant_coding_transcripts.gff3 - 去重合并的GFF3文件")
    print("  2. nonredundant_coding_transcripts.fasta - 去重合并的FASTA文件") 
    print("  3. saturation_curve_data.tsv - 饱和曲线数据")
    print("  4. processing_summary.txt - 处理摘要报告")

if __name__ == "__main__":
    main()
