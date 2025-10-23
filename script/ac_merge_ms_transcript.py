import pandas as pd
from Bio import SeqIO
from datetime import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def filter_by_ms(ms_file, gff_file, pep_file):
    transcript_lines = {}
    current_gene_lines = []
    current_transcript_id = None
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            fields = line.split('\t')                
            feature_type = fields[2]
            attributes = fields[8]
            if feature_type == 'gene':
                if current_transcript_id and current_gene_lines:
                    transcript_lines[current_transcript_id] = current_gene_lines.copy()
                current_gene_lines = [line]
                current_transcript_id = None
            elif feature_type == 'mRNA':
                current_transcript_id = attributes.split(';')[0].split('=')[-1]
                if current_transcript_id:
                    current_gene_lines.append(line)
            elif current_gene_lines:
                current_gene_lines.append(line)
        if current_transcript_id and current_gene_lines:
            transcript_lines[current_transcript_id] = current_gene_lines
    
    df = pd.read_excel(ms_file, sheet_name="NCP")
    transcript_peptide_count = {}
    transcript_cds_regions = {}
    for transcript_id, lines in transcript_lines.items():
        cds_regions = []
        for line in lines:
            if not line.startswith('#') and 'CDS' in line:
                fields = line.split('\t')
                start = int(fields[3])
                end = int(fields[4])
                cds_regions.append((start, end))
        transcript_cds_regions[transcript_id] = sorted(cds_regions)
        transcript_peptide_count[transcript_id] = 0

    for _, peptide_row in df.iterrows():
        peptide_chrom = peptide_row['chrom']
        peptide_start = peptide_row['start']
        peptide_end = peptide_row['end']
        peptide_strand = peptide_row['strand']
        remain_id = []        
        for transcript_id, lines in transcript_lines.items():
            # 从mRNA行提取位置信息
            mrna_line = None
            for line in lines:
                if not line.startswith('#') and 'mRNA' in line:
                    mrna_line = line
                    break            
            if not mrna_line:
                continue                
            fields = mrna_line.split('\t')
            chrom = fields[0]
            strand = fields[6]            
            if chrom != peptide_chrom or strand != peptide_strand:
                continue            
            # 精确检查：肽段是否在CDS区域内
            cds_regions = transcript_cds_regions[transcript_id]
            peptide_in_cds = False            
            for cds_start, cds_end in cds_regions:
                # 检查肽段是否完全在当前CDS区域内
                if (peptide_start >= cds_start and peptide_end <= cds_end) and ((peptide_start - cds_start) % 3 == 0):
                    peptide_in_cds = True
                    remain_id.append(transcript_id)
                    break            
            if peptide_in_cds:
                transcript_peptide_count[transcript_id] += 1    
    # 加载PEP文件序列
    pep_sequences = {}
    for record in SeqIO.parse(pep_file, "fasta"):
        if record.id in remain_id:
            pep_sequences[record.id] = str(record.seq)    
    return transcript_lines, transcript_peptide_count, pep_sequences

def generate_outputs(transcript_lines, transcript_peptide_count, pep_sequences, output_prefix):
    """生成所有要求的输出"""    
    # 1. 生成GFF3格式输出（只输出有肽段支持的转录本）
    gff_output_file = f"{output_prefix}_transcripts_with_peptides.gff3"
    with open(gff_output_file, 'w') as f:
        f.write("##gff-version 3\n")        
        # 遍历所有转录本，只输出有肽段支持的
        for transcript_id, lines in transcript_lines.items():
            if transcript_peptide_count.get(transcript_id, 0) > 0:
                # 更新mRNA行，添加peptide_count属性
                updated_lines = []
                for line in lines:
                    if 'mRNA' in line and not line.startswith('#'):
                        if ';' in line:
                            line = line.rstrip() + f";peptide_count={transcript_peptide_count[transcript_id]}"
                        else:
                            line = line.rstrip() + f"peptide_count={transcript_peptide_count[transcript_id]}"
                    updated_lines.append(line)
                for line in updated_lines:
                    f.write(line + '\n')
                f.write("###\n")
    
    # 2. 生成FASTA文件（有肽段支持的转录本的蛋白质序列）
    fasta_output_file = f"{output_prefix}_transcripts_with_peptides.fasta"
    fasta_records = []
    for transcript_id, seq in pep_sequences.items():
        record = SeqRecord(Seq(seq), id=transcript_id, description="")
        fasta_records.append(record)
    with open(fasta_output_file, "w") as f:
        SeqIO.write(fasta_records, f, "fasta")

    # 3. 统计不同peptide_count的转录本数量
    print("\n=== 不同peptide_count的转录本数量统计 ===")
    count_distribution = {}
    for count in transcript_peptide_count.values():
        count_distribution[count] = count_distribution.get(count, 0) + 1
    
    # 4. 保存统计分布到CSV
    stats_output_file = f"{output_prefix}_peptide_count_distribution.csv"
    df_stats = pd.DataFrame([
        {'peptide_count': count, 'transcript_count': num} 
        for count, num in count_distribution.items()
    ])
    df_stats.to_csv(stats_output_file, index=False)

def main():
    ms_file = r"D:\Desktop\Eu_sp_finally.xlsx"
    gff_file = r"D:\Desktop\output\file1_nonredundant_coding_transcripts.gff3"
    pep_file = r"D:\Desktop\output\file2_nonredundant_coding_transcript_pep.fasta"
    output_prefix = r"D:\Desktop\analysis_results"
    transcript_lines, transcript_peptide_count, pep_sequences = filter_by_ms(ms_file, gff_file, pep_file)
    generate_outputs(transcript_lines, transcript_peptide_count, pep_sequences, output_prefix)

if __name__ == "__main__":
    main()