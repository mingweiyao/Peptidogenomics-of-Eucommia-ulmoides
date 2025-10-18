# # 多外显子gff文件
# import pandas as pd
# from datetime import datetime
# def excel_to_gff3(input_excel, output_gff, sheet_name="CPs_intron"):
#     try:
#         df = pd.read_excel(input_excel, sheet_name=sheet_name)
#         print(f"成功读取Excel文件: {input_excel}")
#         print(f"找到 {len(df)} 条基因记录")
#     except Exception as e:
#         print(f"读取Excel文件失败: {e}")
#         return
#     required_columns = ['ID', 'accessions', 'start', 'end', 'strand', 'chrom']
#     missing_cols = [col for col in required_columns if col not in df.columns]
#     if missing_cols:
#         print(f"错误: 缺少必要的列: {', '.join(missing_cols)}")
#         return
#     gff_header = f"""##gff-version 3
# ##date {datetime.now().strftime('%Y-%m-%d')}
# ##source {input_excel}
# ##genome-build v1.0
# """
#     gff_lines = []
#     for idx, row in df.iterrows():
#         seqid = row['chrom']
#         source = "EuNCP"
#         strand = row['strand']
#         accession = row['accessions']
#         gene_id = row['ID']
#         try:
#             exon1_start, exon1_end = map(int, str(row['start']).split('/'))
#             exon2_start, exon2_end = map(int, str(row['end']).split('/'))
#         except:
#             print(f"警告: 行 {idx+2} 的外显子格式无效，跳过")
#             continue
#         gene_start = min(exon1_start, exon2_start)
#         gene_end = max(exon1_end, exon2_end)
#         gene_line = (
#             f"{seqid}\t{source}\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t"
#             f"ID={gene_id};Name={gene_id}"
#         )
#         gff_lines.append(gene_line)
#         mrna_id = f"{gene_id}.t1"
#         mrna_line = (
#             f"{seqid}\t{source}\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t"
#             f"ID={mrna_id};Parent={gene_id};product=predicted protein"
#         )
#         gff_lines.append(mrna_line)
#         exons = [
#             (exon1_start, exon1_end, 1),
#             (exon2_start, exon2_end, 2)
#         ]
#         if strand == '-':
#             exons.reverse()
#         for exon_num, (exon_start, exon_end, orig_num) in enumerate(exons, 1):
#             exon_id = f"{mrna_id}.exon{exon_num}"
#             exon_line = (
#                 f"{seqid}\t{source}\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t"
#                 f"ID={exon_id};Parent={mrna_id}"
#             )
#             gff_lines.append(exon_line)
#             cds_id = f"{mrna_id}.cds{exon_num}"
#             phase = "." if exon_num != len(exons) else "0"
#             cds_line = (
#                 f"{seqid}\t{source}\tCDS\t{exon_start}\t{exon_end}\t.\t{strand}\t{phase}\t"
#                 f"ID={cds_id};Parent={mrna_id}"
#             )
#             gff_lines.append(cds_line)
#     try:
#         with open(output_gff, 'w') as f:
#             f.write(gff_header)
#             f.write("\n".join(gff_lines))
#         print(f"成功生成GFF3文件: {output_gff}")
#         print(f"转换了 {len(df)} 条基因记录，共生成 {len(gff_lines)} 行GFF记录")
#     except Exception as e:
#         print(f"写入GFF文件失败: {e}")
# if __name__ == "__main__":
#     input_excel = "G:/peptidegenomics/sp_data/output/Eu_sp_finally.xlsx"
#     output_gff = "G:/peptidegenomics/sp_data/output/output_exon.gff"
#     excel_to_gff3(input_excel, output_gff)

import pandas as pd
from collections import defaultdict
from tqdm import tqdm
from datetime import datetime
from intervaltree import Interval, IntervalTree

def filter_overlapping_regions(df):
    filtered_dfs = []
    for chrom, group in df.groupby('chrom'):
        group = group.sort_values('start').reset_index(drop=True)
        filtered_indices = set(range(len(group)))
        tree = IntervalTree()
        for idx, row in group.iterrows():
            tree.addi(row['start'], row['end'], idx)        
        to_remove = set()        
        for idx in range(len(group)):
            if idx in to_remove:
                continue                
            current = group.loc[idx]
            start, end = current['start'], current['end']            
            overlaps = tree.overlap(start, end)
            overlapping_indices = {iv.data for iv in overlaps}            
            if len(overlapping_indices) == 1:
                continue                
            overlapping_indices = sorted(overlapping_indices)
            if len(overlapping_indices) == 2:
                other_idx = overlapping_indices[1] if overlapping_indices[0] == idx else overlapping_indices[0]
                other_row = group.loc[other_idx]                
                if (other_row['end'] - other_row['start']) > (end - start):
                    to_remove.add(idx)
                else:
                    to_remove.add(other_idx)
            else:
                best_remove = None
                max_non_overlap = 0                
                for candidate in overlapping_indices:
                    temp_remove = to_remove | {candidate}
                    remaining = [i for i in overlapping_indices if i not in temp_remove]
                    has_overlap = False
                    for i in range(len(remaining)):
                        for j in range(i+1, len(remaining)):
                            row1 = group.loc[remaining[i]]
                            row2 = group.loc[remaining[j]]
                            if row1['end'] > row2['start']:
                                has_overlap = True
                                break
                        if has_overlap:
                            break                    
                    if not has_overlap:
                        non_overlap_count = len(remaining)
                        if non_overlap_count > max_non_overlap:
                            max_non_overlap = non_overlap_count
                            best_remove = candidate                
                if best_remove is not None:
                    to_remove.add(best_remove)
                else:
                    longest_idx = overlapping_indices[0]
                    max_length = group.loc[longest_idx]['end'] - group.loc[longest_idx]['start']
                    for i in overlapping_indices[1:]:
                        current_length = group.loc[i]['end'] - group.loc[i]['start']
                        if current_length > max_length:
                            longest_idx = i
                            max_length = current_length
                    for i in overlapping_indices:
                        if i != longest_idx:
                            to_remove.add(i)        
        filtered_dfs.append(group[~group.index.isin(to_remove)])    
    return pd.concat(filtered_dfs, ignore_index=True)

def remove_fully_contained(df):
    df = df.copy()
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    filtered_rows = []
    grouped = df.groupby(['chrom', 'strand'], group_keys=False)    
    for (chrom, strand), group in grouped:     
        tree = IntervalTree()
        for idx, row in group.iterrows():
            start, end = row['start'], row['end']
            if start > end:
                start, end = end, start
            if start < end:
                tree[start:end] = idx        
        to_remove = set()
        for interval in tree:
            containing = tree.envelop(interval.begin, interval.end)
            for other in containing:
                if other.data != interval.data:
                    to_remove.add(interval.data)
        filtered_group = group.drop(index=to_remove)
        filtered_rows.append(filtered_group)    
    return pd.concat(filtered_rows, ignore_index=True)

def process_junctions(df, gff_df):
    gff_dict = defaultdict(list)
    for _, row in gff_df.iterrows():
        key = (row['chrom'], row['strand'])
        if int(row['start']) < int(row['end']):
            gff_dict[key].append((
                int(row['start']),
                int(row['end']),
                row['type']
            ))
    cds_trees = defaultdict(IntervalTree)
    for key in gff_dict:
        for start, end, type_ in gff_dict[key]:
            if type_ == 'CDS':
                if start < end:
                    cds_trees[key].addi(start, end)
    to_drop = []
    junction_mask = df['type'] == 'junction'    
    for idx, row in tqdm(df[junction_mask].iterrows(), total=len(df[junction_mask]), desc="Processing junctions"):
        chrom = row['chrom']
        strand = row['strand']
        j_start = int(row['start'])
        j_end = int(row['end'])
        if j_start >= j_end:
            continue
        key = (chrom, strand)
        if key in cds_trees and cds_trees[key].overlaps(j_start, j_end):
            to_drop.append(idx)
    df = df.drop(index=to_drop)
    return df

def excel_to_gff3(df, output_gff):
    required_columns = ['ID', 'accessions', 'start', 'end', 'strand', 'chrom']
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
        accession = row['accessions']
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

def main():
    NCP_file = "/Volumes/caca/test_fractionation/00raw/sp_loc/Eu_sp_finally.xlsx"
    output_gff = "/Volumes/caca/test_fractionation/00raw/sp_loc/Eu_NCPs.gff"
    gff_file = "/Volumes/caca/test_fractionation/00raw/sp_loc/peptide_analysis_results.xlsx"
    try:
        print("步骤1: 读取原始数据...")
        df_NCP = pd.read_excel(NCP_file, sheet_name="GFF")
        df_gff = pd.read_excel(gff_file, sheet_name="Genomic_Features")
        print(f"成功读取 {len(df_NCP)} 行数据")
        
        print("\n步骤2: 移除完全包含的行...")
        df_remove = remove_fully_contained(df_NCP)
        print(f"过滤后剩余 {len(df_remove)} 行数据")
        
        print("\n步骤3: 处理junction区域（去除CDS重叠）...")
        df_r_CDS = process_junctions(df_remove, df_gff)
        
        print("\n步骤4: 处理重叠区域...")
        df_filtered = filter_overlapping_regions(df_r_CDS)
        print(f"重叠处理后剩余 {len(df_filtered)} 行数据")
        
        print("\n步骤5: 生成gff文件")
        excel_to_gff3(df_filtered, output_gff)

        df_filtered.to_excel("/Volumes/caca/test_fractionation/00raw/sp_loc/gff_results_test.xlsx", index=False)
    
    except Exception as e:
        print(f"\n处理过程中出错: {e}")

if __name__ == "__main__":
    main()