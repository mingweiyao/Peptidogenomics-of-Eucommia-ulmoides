import pandas as pd
from Bio import SeqIO
import gffutils
from datetime import datetime

def filter_by_ms(ms_file, gff_file, pep_file):
    db = gffutils.create_db(gff_file, dbfn=":memory:", force=True)
    transcript_cds = {}
    transcript_features = {}
    for transcript in db.features_of_type('mRNA'):
        cds_regions = []
        cds_features = []
        exon_features = []
        utr_features = []
        for cds in db.children(transcript.id, featuretype='CDS'):
            cds_regions.append((cds.start, cds.end))
            cds_features.append(cds)
        for exon in db.children(transcript.id, featuretype='exon'):
            exon_features.append(exon) 
        for utr in db.children(transcript.id, featuretype=['five_prime_UTR', 'three_prime_UTR']):
            utr_features.append(utr)
        transcript_cds[transcript.id] = {
            'chrom': transcript.seqid,
            'strand': transcript.strand,
            'cds_regions': sorted(cds_regions),
            'peptide_count': 0,
            'transcript_feature': transcript,
            'cds_features': cds_features,
            'exon_features': exon_features,
            'utr_features': utr_features
        }
    
    df = pd.read_excel(ms_file, sheet_name="Sheet1")
    for _, peptide_row in df.iterrows():
        peptide_chrom = peptide_row['chrom']
        peptide_start = peptide_row['start']
        peptide_end = peptide_row['end']
        peptide_strand = peptide_row['strand']        
        for transcript_id, transcript_info in transcript_cds.items():
            if (transcript_info['chrom'] != peptide_chrom) or (transcript_info['strand'] != peptide_strand):
                continue                
            for cds_start, cds_end in transcript_info['cds_regions']:
                if (peptide_start >= cds_start and 
                    peptide_end <= cds_end and 
                    (peptide_start - cds_start) % 3 == 0):
                    transcript_info['peptide_count'] += 1
                    break    
    pep_sequences = {}
    for record in SeqIO.parse(pep_file, "fasta"):
        pep_sequences[record.id] = str(record.seq)    
    return transcript_cds, pep_sequences

def generate_gff_output(transcript_cds, output_file, db):    
    with open(output_file, 'w') as f:
        # 写入GFF3文件头
        f.write("##gff-version 3\n")
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        # 按染色体和位置排序输出
        transcripts_with_peptides = [
            (tid, info) for tid, info in transcript_cds.items() 
            if info['peptide_count'] > 0
        ]
        transcripts_with_peptides.sort(key=lambda x: (x[1]['chrom'], x[1]['transcript_feature'].start))    
        for transcript_id, info in transcripts_with_peptides:
            transcript = info['transcript_feature']            
            # 写入基因行（如果存在）
            gene_parents = list(db.parents(transcript_id, featuretype='gene'))
            if gene_parents:
                gene = gene_parents[0]
                gene_attrs = f"ID={gene.id};Name={gene.id};peptide_count={info['peptide_count']}"
                f.write(f"{gene.seqid}\t{gene.source}\t{gene.featuretype}\t{gene.start}\t{gene.end}\t{gene.score}\t{gene.strand}\t{gene.frame}\t{gene_attrs}\n")
            # 写入mRNA行，添加peptide_count属性
            transcript_attrs = transcript.attributes.copy()
            transcript_attrs['peptide_count'] = [str(info['peptide_count'])]
            attrs_str = ';'.join([f"{k}={','.join(v)}" for k, v in transcript_attrs.items()])
            f.write(f"{transcript.seqid}\t{transcript.source}\t{transcript.featuretype}\t{transcript.start}\t{transcript.end}\t{transcript.score}\t{transcript.strand}\t{transcript.frame}\t{attrs_str}\n")
            # 写入外显子
            for exon in info['exon_features']:
                exon_attrs = f"Parent={transcript_id}"
                if 'ID' in exon.attributes:
                    exon_attrs = f"ID={exon.attributes['ID'][0]};{exon_attrs}"
                f.write(f"{exon.seqid}\t{exon.source}\t{exon.featuretype}\t{exon.start}\t{exon.end}\t{exon.score}\t{exon.strand}\t{exon.frame}\t{exon_attrs}\n")
            # 写入CDS
            for cds in info['cds_features']:
                cds_attrs = f"Parent={transcript_id}"
                if 'ID' in cds.attributes:
                    cds_attrs = f"ID={cds.attributes['ID'][0]};{cds_attrs}"
                f.write(f"{cds.seqid}\t{cds.source}\t{cds.featuretype}\t{cds.start}\t{cds.end}\t{cds.score}\t{cds.strand}\t{cds.frame}\t{cds_attrs}\n")
            # 写入UTR
            for utr in info['utr_features']:
                utr_attrs = f"Parent={transcript_id}"
                if 'ID' in utr.attributes:
                    utr_attrs = f"ID={utr.attributes['ID'][0]};{utr_attrs}"
                f.write(f"{utr.seqid}\t{utr.source}\t{utr.featuretype}\t{utr.start}\t{utr.end}\t{utr.score}\t{utr.strand}\t{utr.frame}\t{utr_attrs}\n")
            f.write("###\n")

def generate_outputs(transcript_cds, pep_sequences, db, output_prefix):
    # 1. 生成GFF3格式输出
    gff_output_file = f"{output_prefix}_transcripts_with_peptides.gff3"
    generate_gff_output(transcript_cds, gff_output_file, db)
    # 2. 统计不同peptide_count的转录本数量
    print("\n=== 不同peptide_count的转录本数量统计 ===")
    count_distribution = {}
    for transcript_id, info in transcript_cds.items():
        count = info['peptide_count']
        if count not in count_distribution:
            count_distribution[count] = 0
        count_distribution[count] += 1    
    # 3. 输出有肽段支持的转录本及其PEP序列
    print("\n=== 有肽段支持的转录本及其PEP序列 ===")
    non_zero_transcripts = {tid: info for tid, info in transcript_cds.items() if info['peptide_count'] > 0}
    df_stats = pd.DataFrame([
        {'peptide_count': count, 'transcript_count': num} 
        for count, num in count_distribution.items()
    ])
    df_stats.to_csv(f"{output_prefix}_peptide_count_distribution.csv", index=False)

def main():
    ms_file = "path/to/your/ms_data.xlsx"
    gff_file = "path/to/your/annotation.gff3"
    pep_file = "path/to/your/peptide.fasta"
    output_prefix = "analysis_results"
    transcript_cds, pep_sequences, db = filter_by_ms(ms_file, gff_file, pep_file)
    generate_outputs(transcript_cds, pep_sequences, db, output_prefix)

if __name__ == "__main__":
    main()