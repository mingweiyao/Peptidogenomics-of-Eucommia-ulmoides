from datetime import datetime
import pandas as pd
def excel_to_gff3(df, output_gff):
    required_columns = ['ID', 'start', 'end', 'strand', 'chrom']
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        print(f"错误: 缺少必要的列: {', '.join(missing_cols)}")
        return
    gff_header = f"""##gff-version 3
##date {datetime.now().strftime('%Y-%m-%d')}
##source {df}
##genome-build v1.0
"""
    gff_lines = []
    for idx, row in df.iterrows():
        seqid = row['chrom']
        source = "EuNCP"
        start = int(row['start'])
        end = int(row['end'])
        strand = row['strand']
        gene_id = row['ID']
        gene_line = (
            f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={gene_id};Name={gene_id}"
        )
        gff_lines.append(gene_line)
        mrna_id = f"{gene_id}.t1"
        mrna_line = (
            f"{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={mrna_id};Parent={gene_id};product=predicted protein"
        )
        gff_lines.append(mrna_line)
        exon_id = f"{mrna_id}.exon1"
        exon_line = (
            f"{seqid}\t{source}\texon\t{start}\t{end}\t.\t{strand}\t.\t"
            f"ID={exon_id};Parent={mrna_id}"
        )
        gff_lines.append(exon_line)
        cds_id = f"{mrna_id}.cds"
        cds_line = (
            f"{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
            f"ID={cds_id};Parent={mrna_id}"
        )
        gff_lines.append(cds_line)
    try:
        with open(output_gff, 'w') as f:
            f.write(gff_header)
            f.write("\n".join(gff_lines))
        print(f"成功生成GFF3文件: {output_gff}")
        print(f"转换了 {len(df)} 条基因记录，共生成 {len(gff_lines)} 行GFF记录")
    except Exception as e:
        print(f"写入GFF文件失败: {e}")
excel_file = r"F:\CLE\New_CLE_analysis\理化性质\基因组坐标.xlsx"
output_gff = r"F:\CLE\New_CLE_analysis\染色体定位\EuCLE_partial.gff"
df = pd.read_excel(excel_file, sheet_name="Sheet2")
excel_to_gff3(df, output_gff)