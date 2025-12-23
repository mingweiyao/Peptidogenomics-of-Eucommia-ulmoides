from Bio import SeqIO
import pandas as pd
def file_analysis(genome_file, gff_file, new_genome_file, new_gff_file, excel_info, new_excel_info, tss_file, new_tss_file):
    chrom_ids = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        chrom_ids[record.id] = record.description.split("\t")[1].split("=")[-1]
    new_records = []
    for record in SeqIO.parse(genome_file, "fasta"):
        if record.id in chrom_ids:
            new_id = chrom_ids[record.id]
            record.id = new_id
            record.description = f"{new_id}\tLen={len(record.seq)}"
            new_records.append(record)
    SeqIO.write(new_records, new_genome_file, "fasta")
    print(f"已创建新的基因组文件: {new_genome_file}")
    with open(gff_file, 'r') as infile, open(new_gff_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                seq_id = fields[0]
                if seq_id in chrom_ids:
                    fields[0] = chrom_ids[seq_id]
                    new_line = '\t'.join(fields) + '\n'
                    outfile.write(new_line)
                else:
                    outfile.write(line)
            else:
                outfile.write(line)
    print(f"已创建新的GFF文件: {new_gff_file}")
    df = pd.read_excel(excel_info, sheet_name="unique")
    df['chrom'] = df['chrom'].map(chrom_ids)
    df.to_excel(new_excel_info, sheet_name="unique", index=False)
    tss_df = pd.read_excel(tss_file, sheet_name="TSS")
    tss_df['chrom'] = tss_df['chrom'].map(chrom_ids)
    tss_df.to_excel(new_tss_file, sheet_name="TSS", index=False)

def main():
    genome_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\Eu_genome.fasta"
    gff_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\GWHBISF00000000.gff"
    excel_info = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\finally_expressed_sp_info.xlsx"
    tss_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\figure.xlsx"
    new_genome_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_genome.fasta"
    new_gff_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_genome.gff"
    new_excel_info = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_sp_info.xlsx"
    new_tss_file = r"D:\Desktop\peptidemicro\00file\01figure\Eu_genome_modified\new_figure.xlsx"
    file_analysis(genome_file, gff_file, new_genome_file, new_gff_file, excel_info, new_excel_info, tss_file, new_tss_file)

if __name__ == "__main__":
    main()
