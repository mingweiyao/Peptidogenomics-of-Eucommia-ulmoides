"""
这个脚本包含了用于肽段基因组定位、蛋白质理化性质分析、基因组位置完全不重叠的肽段GFF文件的生成、转录组定量结果合并部分的代码，
满足了图形数据分析时数据的准备
"""

# # 肽段基因组位置定位
# from Bio import SeqIO
# import pandas as pd
# import logging
# import re
# from tqdm import tqdm
# from intervaltree import IntervalTree
# from Bio.SeqUtils.ProtParam import ProteinAnalysis
# import os
# from collections.abc import Mapping
# from openpyxl.utils.dataframe import dataframe_to_rows
# from collections import defaultdict
# #$ log文件设置
# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s - %(levelname)s - %(message)s',
#     handlers=[
#         logging.FileHandler('peptide_analysis.log'),
#         logging.StreamHandler()
#     ]
# )
# logger = logging.getLogger(__name__)
# # 2.加载数据库信息
# def loading_peptide_db_information(peptide_db):
#     peptide_db_dict = {}
#     # >P1 GWHBISF00000001 1 |+| 1/110 1/330
#     try:
#         for rec in SeqIO.parse(peptide_db, "fasta"):
#             db_id = rec.id
#             parts = rec.description.split(" ")
#             chrom = parts[1]
#             strand = parts[3].split("|")[1]
#             frame = parts[2]
#             if strand == "-":
#                 frame = str(int(frame) + 3)
#             start = parts[-1].split("/")[0]
#             end = parts[-1].split("/")[1]
#             peptide_db_dict[db_id] = (chrom, strand, frame, start, end, rec.seq)
#     except Exception as e:
#         logger.error(f"加载数据库失败: {str(e)}")
#     return  peptide_db_dict
# def loading_protein_db_information(protein_db):
#     protein_db_dict = {}
#     '''
#     >GWHPBISF000003	mRNA=GWHTBISF000003	Gene=GWHGBISF000003	Position=GWHBISF00000485: 221065-221100, 231722-231962, 232040-232285, 
#                                     233270-233326, 233451-233711, 234204-234256, 234891-235005, 235094-235348, 238238-239319: +	
#                                     Frame=0	OriID=evm.model.Chr1.10	OriTrascriptID=evm.model.Chr1.10	transl_table=1	OriGeneID=evm.TU.Chr1.10	OriSeqID=Chr1
#     '''
#     chrom_aa_dict = defaultdict(lambda: defaultdict(list))
#     try:
#         for rec in SeqIO.parse(protein_db, "fasta"):
#             positions = []
#             position_list = []
#             mapped_positions = []
#             db_id = rec.id
#             parts = rec.description.split("\t")
#             chrom_part = parts[3].split("=")[1].split(":")
#             chrom = chrom_part[0]
#             strand = chrom_part[-1].split(' ')[-1]
#             position = chrom_part[1].split(", ")
#             for r in position:
#                 start, end = map(int, r.split('-'))
#                 positions.append((start, end))
#             for start, end in positions:
#                 position_list.extend(range(start, end + 1))
#             if strand == "+":
#                 for i, char in enumerate(rec.seq):
#                     mapped_positions.append((char, position_list[i * 3:i * 3 + 3]))
#                     codon_pos = position_list[i * 3:i * 3 + 3]
#                     chrom_aa_dict[chrom][char].append({"pos": codon_pos, "protein_id": db_id})
#             else:
#                 position_list = position_list[:-3][::-1]
#                 for i, char in enumerate(rec.seq):
#                     mapped_positions.append((char, position_list[i * 3:i * 3 + 3]))
#                     codon_pos = position_list[i * 3:i * 3 + 3]
#                     chrom_aa_dict[chrom][char].append({"pos": codon_pos, "protein_id": db_id})
#             protein_db_dict[db_id] = (chrom, strand, mapped_positions, rec.seq)
#     except Exception as e:
#         logger.error(f"加载数据库失败: {str(e)}")
#     chrom_aa_dict_final = {
#         chrom: dict(aa_pos) 
#         for chrom, aa_pos in chrom_aa_dict.items()
#     }
#     return protein_db_dict, chrom_aa_dict_final
# # 3.染色体长度统计
# def genome_chrom_length(genome_file):
#     # >GWHBISF00000001	OriSeqID=F000542F	Len=330
#     chrom_length_dict = {}
#     try:
#         for rec in SeqIO.parse(genome_file, "fasta"):
#             match = re.search(r'Len=(\d+)', rec.description)
#             if match:
#                 chrom_length_dict[rec.id] = int(match.group(1))
#         logger.info(f"获取{len(chrom_length_dict)}条染色体长度信息")
#     except Exception as e:
#         logger.error(f"获取染色体长度失败: {str(e)}")
#         raise
#     return chrom_length_dict
# # 4. 计算肽段坐标
# def calculate_peptide_coordinates(peptides, peptide_db_info, chrom_length_dict):
#     peptides = peptides.copy()
#     peptides["start"] = None
#     peptides["end"] = None
#     peptides["strand"] = None
#     peptides["chrom"] = None
#     peptides["frame"] = None
#     for index, row in tqdm(peptides.iterrows(), total=len(peptides), desc="计算肽段坐标"):
#         try:
#             id = row['accessions']
#             seq = row["sequence"]
#             if id not in peptide_db_info:
#                 continue
#             db_chrom, db_strand, db_frame, db_start, db_end, db_seq = peptide_db_info[id]
#             chrom_length = chrom_length_dict[db_chrom]
#             pos = db_seq.find(seq)
#             if pos == -1:
#                 continue
#             if db_strand == "+":
#                 start = int(db_start) + pos * 3
#                 end = start + len(seq) * 3 - 1
#             else:
#                 start_strand = int(db_start) + pos * 3
#                 end_strand = start_strand + len(seq) * 3 - 1
#                 end = chrom_length - start_strand + 1
#                 start = chrom_length - end_strand + 1
#             peptides.at[index, 'start'] = start
#             peptides.at[index, 'end'] = end
#             peptides.at[index, 'strand'] = db_strand
#             peptides.at[index, 'chrom'] = db_chrom
#             peptides.at[index, 'frame'] = db_frame
#         except Exception as e:
#             logger.warning(f"行 {index} 处理失败: {str(e)}")
#             continue
#     return peptides
# def calculate_protein_coordinates(proteins, protein_db_info):
#     proteins = proteins.copy()
#     proteins["start"] = None
#     proteins["end"] = None
#     proteins["strand"] = None
#     proteins["chrom"] = None
#     proteins["frame"] = 0
#     for index, row in tqdm(proteins.iterrows(), total=len(proteins), desc="计算蛋白坐标"):
#         try:
#             protein_id = row['accessions']
#             seq = row["sequence"]
#             if protein_id not in protein_db_info:
#                 continue
#             db_chrom, db_strand, db_mapped_positions, db_seq = protein_db_info[protein_id]
#             seq_len = len(seq)
#             pos = db_seq.find(seq)
#             if pos == -1 or (pos + seq_len) > len(db_mapped_positions):
#                 continue
#             first_pos = db_mapped_positions[pos][1]
#             last_pos = db_mapped_positions[pos + seq_len - 1][1]
#             if db_strand == "+":
#                 start = first_pos[0]
#                 end = last_pos[-1]
#             else:
#                 start = last_pos[-1]
#                 end = first_pos[0]
#             proteins.at[index, 'start'] = start
#             proteins.at[index, 'end'] = end
#             proteins.at[index, 'strand'] = db_strand
#             proteins.at[index, 'chrom'] = db_chrom   
#         except Exception as e:
#             logger.warning(f"行 {index} 处理失败: {str(e)}")
#             continue       
#     return proteins
# # 5. 解析gff文件
# def region_calculate(mrnas, exons):
#     required_cols = ['chrom', 'type', 'strand', 'start', 'end']
#     if not all(col in mrnas.columns for col in required_cols):
#         raise ValueError("mRNAs DataFrame缺少必要列")
#     if not all(col in exons.columns for col in required_cols):
#         raise ValueError("Exons DataFrame缺少必要列")
#     results = []
#     mrnas_group = mrnas.groupby(['chrom', 'strand'])
#     exons_group = exons.groupby(['chrom', 'strand'])
#     for (chrom, strand), mrna_group in mrnas_group:
#         tree = IntervalTree()    
#         for _, row in mrna_group.iterrows():
#             tree.addi(row["start"], row["end"] + 1)
#         if (chrom, strand) in exons_group.groups:
#             for _, exon_row in exons_group.get_group((chrom, strand)).iterrows():
#                 tree.chop(exon_row["start"], exon_row["end"] + 1)
#         for interval in sorted(tree):
#             results.append({
#                 "chrom": chrom,
#                 "type": "intron",
#                 "start": interval.begin,
#                 "end": interval.end - 1,
#                 "strand": strand
#             })   
#     return pd.DataFrame(results)
# def intergenic_calculate(chrom_length_dict, genes, strand):
#     if not isinstance(chrom_length_dict, dict):
#         raise TypeError("chrom_length_dict应为字典")
#     required_cols = ['chrom', 'start', 'end', 'strand']
#     if not all(col in genes.columns for col in required_cols):
#         raise ValueError(f"genes DataFrame缺少必要列，需要: {required_cols}")
#     if strand not in ('+', '-'):
#         raise ValueError("strand参数必须是'+'或'-'")
#     intergenic_regions = []
#     for chrom, length in chrom_length_dict.items():
#         try:
#             chrom_length = int(length)
#             if chrom_length <= 0:
#                 raise ValueError(f"染色体{chrom}长度必须为正整数")
#         except (ValueError, TypeError):
#             logger.warning(f"忽略无效的染色体长度: {chrom}={length}")
#             continue
#         tree = IntervalTree()
#         tree.addi(1, chrom_length + 1)
#         chrom_genes = genes[(genes['chrom'] == chrom)]
#         if not chrom_genes.empty:
#             for _, gene in chrom_genes.iterrows():
#                 try:
#                     tree.chop(gene["start"], gene["end"] + 1)
#                 except ValueError as e:
#                     logger.warning(f"忽略无效基因区间 {chrom}:{gene['start']}-{gene['end']}: {str(e)}")
#                     continue
#         tree.merge_overlaps(strict=False)
#         for interval in sorted(tree):
#             intergenic_regions.append({
#                 "chrom": chrom,
#                 "type": "intergenic",
#                 "start": interval.begin,
#                 "end": interval.end - 1,
#                 "strand": strand
#             })
#     return pd.DataFrame(intergenic_regions)
# def parse_gff_file(gff_file, chrom_length_dict):
#     try:
#         gff_data = []
#         with open(gff_file, "r") as f:
#             for line_num, line in enumerate(f, 1):
#                 line = line.strip()
#                 if not line or line.startswith("#"):
#                     continue
#                 parts = line.split("\t")
#                 if len(parts) < 9:
#                     logger.warning(f"忽略不完整的行(列数不足): 第{line_num}行")
#                     continue
#                 try:
#                     record = {
#                         'chrom': parts[0],
#                         'type': parts[2],
#                         'start': int(parts[3]),
#                         'end': int(parts[4]),
#                         'strand': parts[6]
#                     }
#                     if record['start'] > record['end']:
#                         logger.warning(f"忽略无效坐标(start>end): 第{line_num}行")
#                         continue
#                     gff_data.append(record)
#                 except (ValueError, IndexError) as e:
#                     logger.warning(f"忽略格式错误的行: 第{line_num}行: {str(e)}")
#                     continue
#         gff_df = pd.DataFrame(gff_data)
#         if gff_df.empty:
#             logger.error("GFF文件无有效数据")
#             return {}
#         def filter_features(df, feature_type, strand=None):
#             mask = df['type'] == feature_type
#             if strand is not None:
#                 mask &= df['strand'] == strand
#             return df[mask].copy()
#         mrnas_pos = filter_features(gff_df, 'mRNA', '+')
#         mrnas_neg = filter_features(gff_df, 'mRNA', '-')
#         exons_pos = filter_features(gff_df, 'exon', '+')
#         exons_neg = filter_features(gff_df, 'exon', '-')
#         cds_pos = filter_features(gff_df, 'CDS', '+')
#         cds_neg = filter_features(gff_df, 'CDS', '-')
#         utr3_pos = filter_features(gff_df, 'three_prime_UTR', '+')
#         utr3_neg = filter_features(gff_df, 'three_prime_UTR', '-')
#         utr5_pos = filter_features(gff_df, 'five_prime_UTR', '+')
#         utr5_neg = filter_features(gff_df, 'five_prime_UTR', '-')
#         intron_pos = region_calculate(mrnas_pos, exons_pos)
#         intron_neg = region_calculate(mrnas_neg, exons_neg)
#         intergenic_pos = intergenic_calculate(chrom_length_dict, mrnas_pos, "+")
#         intergenic_neg = intergenic_calculate(chrom_length_dict, mrnas_neg, "-")
#         features = pd.concat([
#             cds_pos, cds_neg,
#             utr3_pos, utr3_neg,
#             utr5_pos, utr5_neg,
#             intron_pos, intron_neg,
#             intergenic_pos, intergenic_neg
#         ], ignore_index=True)
#         feature_dict = {}
#         for _, row in features.iterrows():
#             key = (row['chrom'], row['strand'])
#             feature_dict.setdefault(key, []).append(
#                 (row['start'], row['end'], row['type'])
#             )
#         logger.info(
#             f"成功解析GFF文件: 共处理{len(gff_df)}条记录, "
#             f"得到{len(feature_dict)}个染色体/链组合的特征"
#         )
#         return feature_dict
#     except Exception as e:
#         logger.error(f"解析GFF文件失败: {str(e)}")
#         raise
# # 6. 过滤N区域
# def detect_n_regions(genome_file):
#     n_regions = []
#     try:
#         for rec in SeqIO.parse(genome_file, "fasta"):
#             chrom = rec.description.split()[0]
#             seq = str(rec.seq).upper()
#             n_start = None
#             for i, base in enumerate(seq, 1):
#                 if base == 'N':
#                     if n_start is None:
#                         n_start = i
#                 else:
#                     if n_start is not None:
#                         n_regions.append([chrom, n_start, i-1])
#                         n_start = None
#             if n_start is not None:
#                 n_regions.append([chrom, n_start, len(seq)])
#         n_df = pd.DataFrame(n_regions, columns=['chrom', 'start', 'end'])
#         logger.info(f"发现 {len(n_df)} 个N区域")
#         return n_df
#     except Exception as e:
#         logger.error(f"N区域检测失败: {str(e)}")
#         raise
# def filter_n_overlaps(peptide_df, n_regions):
#     required_cols = ['chrom', 'start', 'end']
#     if not all(col in peptide_df.columns for col in required_cols):
#         raise ValueError("肽段DataFrame缺少必要列")
#     if not all(col in n_regions.columns for col in required_cols):
#         raise ValueError("N区域DataFrame缺少必要列")
#     from intervaltree import IntervalTree
#     chrom_trees = {}
#     for _, row in n_regions.iterrows():
#         chrom = row['chrom']
#         if chrom not in chrom_trees:
#             chrom_trees[chrom] = IntervalTree()
#         chrom_trees[chrom].addi(row['start'], row['end']+1)
#     peptide_df = peptide_df.copy()
#     peptide_df['overlaps_n'] = False
#     for idx, row in tqdm(peptide_df.iterrows(), total=len(peptide_df), desc="过滤N区域"):
#         chrom = row['chrom']
#         start = row['start']
#         end = row['end']
#         if pd.isna(start) or pd.isna(end) or chrom not in chrom_trees:
#             continue
#         if chrom_trees[chrom].overlaps(start, end+1):
#             peptide_df.at[idx, 'overlaps_n'] = True
#     filtered = peptide_df[~peptide_df['overlaps_n']].copy()
#     removed = peptide_df[peptide_df['overlaps_n']].copy()
#     filtered.drop(columns=['overlaps_n'], inplace=True)
#     removed.drop(columns=['overlaps_n'], inplace=True)
#     logger.info(f"过滤结果: 保留 {len(filtered)} 条, 移除 {len(removed)} 条")
#     return filtered, removed
# # 7. 注释肽段基因组位置
# def annotate_peptide_loc(gff_info, peptide_file_coord):
#     required_cols = ['chrom', 'start', 'end', 'strand']
#     if not all(col in peptide_file_coord.columns for col in required_cols):
#         raise ValueError("肽段坐标DataFrame缺少必要列")
#     type_list = []
#     for _, row in tqdm(peptide_file_coord.iterrows(), total=len(peptide_file_coord), desc="注释肽段基因组位置"):
#         chrom = row['chrom']
#         start = row['start']
#         end = row['end']
#         strand = row['strand']
#         search_key = (chrom, strand)
#         current_type = None
#         junction_found = False
#         if search_key in gff_info:
#             features = gff_info[search_key]
#             for feat_start, feat_end, feat_type in features:
#                 if feat_start <= start < end <= feat_end:
#                     current_type = feat_type
#                     break
#                 elif (start < feat_start < end) or (start < feat_end < end):
#                     junction_found = True
#             if current_type is None and junction_found:
#                 current_type = "junction"
#         type_list.append(current_type)
#     peptide_file_coord = peptide_file_coord.copy()
#     peptide_file_coord['type'] = type_list
#     return peptide_file_coord
# def annotate_protein_loc(protein_file_coord, protein_db):
#     position_dict = {}
#     for record in SeqIO.parse(protein_db, "fasta"):
#         pos_info = record.description.split("\t")[3].split(": ")
#         position_info = pos_info[-2]
#         strand = pos_info[-1]
#         ranges = position_info.split(',')
#         boundaries = []
#         for r in ranges:
#             start, end = r.strip().split('-')
#             start = int(start)
#             end = int(end)
#             boundaries.append((start, end))
#         position_dict[record.id] = {
#             'strand': strand,
#             'boundaries': sorted(boundaries)
#         }
#     type_list = []
#     first_boundary = []
#     second_boundary = []
#     for _, row in tqdm(protein_file_coord.iterrows(), total=len(protein_file_coord), desc="注释肽段基因组位置"):
#         start = row['start']
#         end = row['end']
#         accession = row['accessions']
#         chrom_data = position_dict[accession]
#         boundaries = chrom_data['boundaries']
#         if start != None:
#             type_list.append("CPs")
#             flag = 0
#             for i in range(len(boundaries)-1):
#                 if (start <= boundaries[i][1]) and (end >= boundaries[i+1][0]):
#                     first_boundary.append(boundaries[i][1])
#                     second_boundary.append(boundaries[i+1][0])
#                     flag = 1
#                     break
#             if flag==0:
#                 first_boundary.append(None)
#                 second_boundary.append(None)
#         else:
#             type_list.append("error")
#             first_boundary.append(None)
#             second_boundary.append(None)
#     protein_file_coord['type'] = type_list
#     protein_file_coord['first_boundary'] = first_boundary
#     protein_file_coord['second_boundary'] = second_boundary
#     return protein_file_coord
# # 8. 分析CDS肽段
# def cds_peptide_analysis(annotated_peptides, protein_db_info, chrom_aa_dict):
#     required_cols = ['type', 'accessions', 'sequence', 'chrom', 'start', 'end', 'strand']
#     if not all(col in annotated_peptides.columns for col in required_cols):
#         raise ValueError(f"输入DataFrame缺少必要列: {required_cols}")
#     cds_types = pd.Series(index=annotated_peptides.index, dtype=object)
#     for idx, row in tqdm(annotated_peptides.iterrows(),
#                         total=len(annotated_peptides),
#                         desc="CDS肽段验证"):
#         if row['type'] != 'CDS':
#             cds_types[idx] = None
#             continue
#         peptide_id = row['accessions']
#         seq = row['sequence']
#         chrom = row['chrom']
#         start = row['start']
#         end = row['end']
#         strand = row['strand']
#         if peptide_id in protein_db_info:
#             db_chrom, db_strand, mapped_position, db_seq = protein_db_info[peptide_id]
#             if strand != db_strand:
#                 cds_types[idx] = "CDS(strand mismatch)"
#                 continue
#             index = db_seq.find(seq)
#             if index == -1:
#                 cds_types[idx] = "CDS(out of frame)"
#                 continue
#             if strand == "+":
#                 seq_start = mapped_position[index][1][0]
#                 seq_end = mapped_position[index + len(seq) - 1][1][-1]
#             else:
#                 seq_start = mapped_position[index + len(seq) - 1][1][-1]
#                 seq_end = mapped_position[index][1][0]    
#             if seq_start == start and seq_end == end:
#                 cds_types[idx] = "CDS(verified by protein_db)"
#             else:
#                 cds_types[idx] = "CDS(out of frame)"
#         elif chrom in chrom_aa_dict:
#             first_aa = seq[0]
#             matched = False
#             if first_aa in chrom_aa_dict[chrom]:
#                 for pos_info in chrom_aa_dict[chrom][first_aa]:
#                     pos_group = pos_info["pos"]
#                     ref_pos = pos_group[0] if strand == '+' else pos_group[-1]
#                     if ref_pos == start:
#                         matched_protein_id = pos_info["protein_id"] 
#                         if matched_protein_id in protein_db_info:
#                             db_chrom, db_strand, mapped_position, db_seq = protein_db_info[matched_protein_id]
#                             if strand != db_strand:
#                                 cds_types[idx] = "CDS(strand mismatch)"
#                                 matched = True
#                                 break
#                             index = db_seq.find(seq)
#                             if index == -1:
#                                 cds_types[idx] = "CDS(out of frame)"
#                                 matched = True
#                                 break
#                             if strand == "+":
#                                 seq_start = mapped_position[index][1][0]
#                                 seq_end = mapped_position[index + len(seq) - 1][1][-1]
#                             else:
#                                 seq_start = mapped_position[index + len(seq) - 1][1][-1]
#                                 seq_end = mapped_position[index][1][0]                                  
#                             if seq_start == start and seq_end == end:
#                                 cds_types[idx] = "CDS(verified by protein_db)"
#                             else:
#                                 cds_types[idx] = "CDS(out of frame)"
#                             matched = True
#                             break                
#                 if not matched:
#                     cds_types[idx] = "CDS(out of frame)"
#             else:
#                 cds_types[idx] = "CDS(out of frame)"
#         else:
#             cds_types[idx] = "CDS(out of frame)"
#     result_df = annotated_peptides.copy()
#     result_df['cds_type'] = cds_types
#     return result_df
# # 9. 分析肽段性质
# def analyze_peptide_properties(sequence):
#     if not sequence or not isinstance(sequence, str):
#         return {
#             "length": 0,
#             "molecular_weight": None,
#             "isoelectric_point": None,
#             "gravy": None,
#             "aromaticity": None,
#             "instability_index": None
#         }
#     try:
#         analysis = ProteinAnalysis(sequence)
#         properties = {
#             "length": len(sequence),
#             "molecular_weight": analysis.molecular_weight(),
#             "isoelectric_point": analysis.isoelectric_point(),
#             "gravy": analysis.gravy(),
#         }
#         try:
#             properties.update({
#                 "instability_index": analysis.instability_index()
#             })
#         except Exception as e:
#             logger.warning(f"部分性质计算失败: {str(e)}")
#         return properties  
#     except Exception as e:
#         logger.error(f"肽段分析失败: {str(e)}")
#         return {
#             "length": len(sequence),
#             "molecular_weight": None,
#             "isoelectric_point": None,
#             "gravy": None,
#             "aromaticity": None,
#             "instability_index": None
#         }
# def batch_analyze_peptides(annotated_peptides_finally):
#     results = []
#     for _, row in tqdm(annotated_peptides_finally.iterrows(), 
#                           total=len(annotated_peptides_finally),
#                           desc="分析肽段性质"):
#         seq = row.get('sequence', '')
#         results.append(analyze_peptide_properties(seq))
#     result_df = pd.DataFrame(results)
#     output_df = pd.concat([annotated_peptides_finally, result_df], axis=1)
#     return output_df
# # 10. 输出文件
# def save_to_excel(data_dict, output_file):
#     try:
#         os.makedirs(os.path.dirname(output_file), exist_ok=True)
#         with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
#             for sheet_name, data in data_dict.items():
#                 if isinstance(data, Mapping):
#                     if all(isinstance(v, (list, tuple)) for v in data.values()):
#                         rows = []
#                         for key, values in data.items():
#                             chrom, strand = key
#                             for (start, end, feat_type) in values:
#                                 rows.append([chrom, strand, start, end, feat_type])
#                         df = pd.DataFrame(rows, columns=["chrom", "strand", "start", "end", "type"])
#                     else:
#                         df = pd.DataFrame.from_dict(data, orient='index')
#                 else:
#                     df = data if isinstance(data, pd.DataFrame) else pd.DataFrame()
#                 if df.empty:
#                     df = pd.DataFrame({"Info": [f"No data in {sheet_name}"]})
#                 safe_name = re.sub(r'[\\/*?\[\]:]', '_', str(sheet_name))[:31]
#                 df.to_excel(
#                     writer,
#                     sheet_name=safe_name,
#                     index=False,
#                     header=not df.empty,
#                     freeze_panes=(1, 0)
#                 )
#                 if not df.empty:
#                     worksheet = writer.sheets[safe_name]
#                     for col in worksheet.columns:
#                         max_len = max((
#                             len(str(cell.value)) 
#                             for cell in col
#                         ), default=0) + 2
#                         worksheet.column_dimensions[col[0].column_letter].width = min(max_len, 50)
#         logger.info(f"成功保存结果到: {output_file}")
#         return True
#     except Exception as e:
#         logger.error(f"保存失败: {str(e)}")
#         if os.path.exists(output_file):
#             os.remove(output_file)
#         return False
# # 整体工作流程
# def run_workflow(peptide_db, protein_db, peptide_file, protein_file, gff_file, genome_file, output_dir):
#     logger.info("=== 开始肽段分析流程 ===")
#     try:
#         logger.info("1. 加载肽段数据...")
#         peptides = pd.read_excel(peptide_file)
#         peptides.columns = peptides.columns.str.strip().str.lower()
#         proteins = pd.read_excel(protein_file)
#         proteins.columns = proteins.columns.str.strip().str.lower()
#         logger.info("2. 加载数据库数据...")
#         peptide_db_info = loading_peptide_db_information(peptide_db)
#         protein_db_info, chrom_aa_dict = loading_protein_db_information(protein_db)
#         logger.info("3. 染色体长度统计")
#         chrom_length_dict = genome_chrom_length(genome_file)
#         logger.info("4. 计算肽段坐标")
#         peptide_file_coord = calculate_peptide_coordinates(peptides, peptide_db_info, chrom_length_dict)
#         protein_file_coord = calculate_protein_coordinates(proteins, protein_db_info)
#         logger.info("5. 解析gff文件")
#         gff_info = parse_gff_file(gff_file, chrom_length_dict)
#         logger.info('6. 过滤N区域')
#         n_regions = detect_n_regions(genome_file)
#         filtered, removed = filter_n_overlaps(peptide_file_coord, n_regions)
#         logger.info("7. 注释肽段基因组位置")
#         annotated_peptides = annotate_peptide_loc(gff_info, filtered)
#         annotated_proteins = annotate_protein_loc(protein_file_coord, protein_db)
#         logger.info("8. 分类CDS肽段")
#         annotated_peptides_finally = cds_peptide_analysis(annotated_peptides, protein_db_info, chrom_aa_dict)
#         annotated_proteins_finally = cds_peptide_analysis(annotated_proteins, protein_db_info, chrom_aa_dict)
#         logger.info("9. 分析肽段性质")
#         peptide_finally = batch_analyze_peptides(annotated_peptides_finally)
#         removed_n = batch_analyze_peptides(removed)
#         protein_finally = batch_analyze_peptides(annotated_proteins_finally)
#         logger.info("10. 保存结果")
#         output_data = {
#             "Annotated_Peptides": peptide_finally,
#             "Annotated_Proteins": protein_finally,
#             "Genomic_Features": gff_info,
#             "removed_n": removed_n
#         }
#         output_file = os.path.join(output_dir, "peptide_analysis_results.xlsx")
#         save_to_excel(output_data, output_file)        
#         logger.info(f"===分析流程完成 ===")
#     except Exception as e:
#         logger.error(f"流程执行失败：{str(e)}")
#         raise
# if __name__ == "__main__":
#     input_files = {
#         "peptide_db": "/Volumes/caca/test_fractionation/00raw/Eu_peptide_database_customized_5.fa",
#         "protein_db": "/Volumes/caca/test_fractionation/00raw/Eu_peptide_database.fa",
#         "peptide_file": "/Volumes/caca/test_fractionation/00raw/Eu_peptide_raw.xlsx",
#         "protein_file": "/Volumes/caca/test_fractionation/00raw/Eu_protein_raw.xlsx",
#         "gff_file": "/Volumes/caca/test_fractionation/00raw/GWHBISF00000000.gff",
#         "genome_file": "/Volumes/caca/test_fractionation/00raw/Eu_genome.fasta",
#         "output_dir": "/Volumes/caca/test_fractionation/00raw/output_test"
#     }
#     try:
#         run_workflow(**input_files)
#     except Exception as e:
#         logger.critical(f"程序终止：{str(e)}")


# # 蛋白质性质预测错误的肽段重新预测
# import pandas as pd
# from Bio.SeqUtils.ProtParam import ProteinAnalysis
# from tqdm import tqdm
# import logging
# logging.basicConfig(level=logging.WARNING)
# logger = logging.getLogger()
# def analyze_peptide_properties(sequence):
#     if not sequence or not isinstance(sequence, str):
#         return {
#             "length": 0,
#             "molecular_weight": None,
#             "isoelectric_point": None,
#             "gravy": None,
#             "aromaticity": None,
#             "instability_index": None
#         }
#     try:
#         analysis = ProteinAnalysis(sequence)
#         properties = {
#             "length": len(sequence),
#             "molecular_weight": analysis.molecular_weight(),
#             "isoelectric_point": analysis.isoelectric_point(),
#             "gravy": analysis.gravy(),
#         }
#         try:
#             properties.update({
#                 "instability_index": analysis.instability_index()
#             })
#         except Exception as e:
#             logger.warning(f"部分性质计算失败: {str(e)}")
#         return properties  
#     except Exception as e:
#         logger.error(f"肽段分析失败: {str(e)}")
#         return {
#             "length": len(sequence),
#             "molecular_weight": None,
#             "isoelectric_point": None,
#             "gravy": None,
#             "aromaticity": None,
#             "instability_index": None
#         }
# def batch_analyze_peptides_from_excel(file_path):
#     df = pd.read_excel(file_path, sheet_name="NCP")
#     if 'sequence' not in df.columns:
#         logger.error("Excel 文件中未找到 'sequence' 列")
#         return None
#     results = []
#     for _, row in tqdm(df.iterrows(), total=len(df), desc="分析肽段性质"):
#         seq = row.get('sequence', '')
#         results.append(analyze_peptide_properties(seq))
#     result_df = pd.DataFrame(results)
#     return result_df
# file_path = r'G:\test_fractionation\00raw\output\Eu_sp_finally.xlsx'
# output_df = batch_analyze_peptides_from_excel(file_path)
# if output_df is not None:
#     output_df.to_excel(r'G:\test_fractionation\00raw\output\output_with_peptide_properties.xlsx', index=False)
#     print("分析结果已保存到 'output_with_peptide_properties.xlsx'")

# # 生成区间完全不重叠的肽段的GFF文件
# import pandas as pd
# from collections import defaultdict
# from tqdm import tqdm
# from datetime import datetime
# from intervaltree import Interval, IntervalTree
# def filter_overlapping_regions(df):
#     filtered_dfs = []
#     for chrom, group in df.groupby('chrom'):
#         group = group.sort_values('start').reset_index(drop=True)
#         filtered_indices = set(range(len(group)))
#         tree = IntervalTree()
#         for idx, row in group.iterrows():
#             tree.addi(row['start'], row['end'], idx)        
#         to_remove = set()        
#         for idx in range(len(group)):
#             if idx in to_remove:
#                 continue                
#             current = group.loc[idx]
#             start, end = current['start'], current['end']            
#             overlaps = tree.overlap(start, end)
#             overlapping_indices = {iv.data for iv in overlaps}            
#             if len(overlapping_indices) == 1:
#                 continue                
#             overlapping_indices = sorted(overlapping_indices)
#             if len(overlapping_indices) == 2:
#                 other_idx = overlapping_indices[1] if overlapping_indices[0] == idx else overlapping_indices[0]
#                 other_row = group.loc[other_idx]                
#                 if (other_row['end'] - other_row['start']) > (end - start):
#                     to_remove.add(idx)
#                 else:
#                     to_remove.add(other_idx)
#             else:
#                 best_remove = None
#                 max_non_overlap = 0                
#                 for candidate in overlapping_indices:
#                     temp_remove = to_remove | {candidate}
#                     remaining = [i for i in overlapping_indices if i not in temp_remove]
#                     has_overlap = False
#                     for i in range(len(remaining)):
#                         for j in range(i+1, len(remaining)):
#                             row1 = group.loc[remaining[i]]
#                             row2 = group.loc[remaining[j]]
#                             if row1['end'] > row2['start']:
#                                 has_overlap = True
#                                 break
#                         if has_overlap:
#                             break                    
#                     if not has_overlap:
#                         non_overlap_count = len(remaining)
#                         if non_overlap_count > max_non_overlap:
#                             max_non_overlap = non_overlap_count
#                             best_remove = candidate                
#                 if best_remove is not None:
#                     to_remove.add(best_remove)
#                 else:
#                     longest_idx = overlapping_indices[0]
#                     max_length = group.loc[longest_idx]['end'] - group.loc[longest_idx]['start']
#                     for i in overlapping_indices[1:]:
#                         current_length = group.loc[i]['end'] - group.loc[i]['start']
#                         if current_length > max_length:
#                             longest_idx = i
#                             max_length = current_length
#                     for i in overlapping_indices:
#                         if i != longest_idx:
#                             to_remove.add(i)        
#         filtered_dfs.append(group[~group.index.isin(to_remove)])    
#     return pd.concat(filtered_dfs, ignore_index=True)
# def remove_fully_contained(df):
#     df = df.copy()
#     df['start'] = df['start'].astype(int)
#     df['end'] = df['end'].astype(int)
#     filtered_rows = []
#     grouped = df.groupby(['chrom', 'strand'], group_keys=False)    
#     for (chrom, strand), group in grouped:     
#         tree = IntervalTree()
#         for idx, row in group.iterrows():
#             start, end = row['start'], row['end']
#             if start > end:
#                 start, end = end, start
#             if start < end:
#                 tree[start:end] = idx        
#         to_remove = set()
#         for interval in tree:
#             containing = tree.envelop(interval.begin, interval.end)
#             for other in containing:
#                 if other.data != interval.data:
#                     to_remove.add(interval.data)
#         filtered_group = group.drop(index=to_remove)
#         filtered_rows.append(filtered_group)    
#     return pd.concat(filtered_rows, ignore_index=True)
# def process_junctions(df, gff_df):
#     gff_dict = defaultdict(list)
#     for _, row in gff_df.iterrows():
#         key = (row['chrom'], row['strand'])
#         if int(row['start']) < int(row['end']):
#             gff_dict[key].append((
#                 int(row['start']),
#                 int(row['end']),
#                 row['type']
#             ))
#     cds_trees = defaultdict(IntervalTree)
#     for key in gff_dict:
#         for start, end, type_ in gff_dict[key]:
#             if type_ == 'CDS':
#                 if start < end:
#                     cds_trees[key].addi(start, end)
#     to_drop = []
#     junction_mask = df['type'] == 'junction'    
#     for idx, row in tqdm(df[junction_mask].iterrows(), total=len(df[junction_mask]), desc="Processing junctions"):
#         chrom = row['chrom']
#         strand = row['strand']
#         j_start = int(row['start'])
#         j_end = int(row['end'])
#         if j_start >= j_end:
#             continue
#         key = (chrom, strand)
#         if key in cds_trees and cds_trees[key].overlaps(j_start, j_end):
#             to_drop.append(idx)
#     df = df.drop(index=to_drop)
#     return df
# def excel_to_gff3(df, output_gff):
#     required_columns = ['ID', 'accessions', 'start', 'end', 'strand', 'chrom']
#     missing_cols = [col for col in required_columns if col not in df.columns]
#     if missing_cols:
#         print(f"错误: 缺少必要的列: {', '.join(missing_cols)}")
#         return
#     gff_header = f"""##gff-version 3
# ##date {datetime.now().strftime('%Y-%m-%d')}
# ##source {df}
# ##genome-build v1.0
# """
#     gff_lines = []
#     for idx, row in df.iterrows():
#         seqid = row['chrom']
#         source = "EuNCP"
#         start = int(row['start'])
#         end = int(row['end'])
#         strand = row['strand']
#         accession = row['accessions']
#         gene_id = row['ID']
#         gene_line = (
#             f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={gene_id};Name={gene_id}"
#         )
#         gff_lines.append(gene_line)
#         mrna_id = f"{gene_id}.t1"
#         mrna_line = (
#             f"{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={mrna_id};Parent={gene_id};product=predicted protein"
#         )
#         gff_lines.append(mrna_line)
#         exon_id = f"{mrna_id}.exon1"
#         exon_line = (
#             f"{seqid}\t{source}\texon\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={exon_id};Parent={mrna_id}"
#         )
#         gff_lines.append(exon_line)
#         cds_id = f"{mrna_id}.cds"
#         cds_line = (
#             f"{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
#             f"ID={cds_id};Parent={mrna_id}"
#         )
#         gff_lines.append(cds_line)
#     try:
#         with open(output_gff, 'w') as f:
#             f.write(gff_header)
#             f.write("\n".join(gff_lines))
#         print(f"成功生成GFF3文件: {output_gff}")
#         print(f"转换了 {len(df)} 条基因记录，共生成 {len(gff_lines)} 行GFF记录")
#     except Exception as e:
#         print(f"写入GFF文件失败: {e}")
# def main():
#     NCP_file = "/Volumes/caca/test_fractionation/00raw/sp_loc/Eu_sp_finally.xlsx"
#     output_gff = "/Volumes/caca/test_fractionation/00raw/sp_loc/Eu_NCPs.gff"
#     gff_file = "/Volumes/caca/test_fractionation/00raw/sp_loc/peptide_analysis_results.xlsx"
#     try:
#         print("步骤1: 读取原始数据...")
#         df_NCP = pd.read_excel(NCP_file, sheet_name="GFF")
#         df_gff = pd.read_excel(gff_file, sheet_name="Genomic_Features")
#         print(f"成功读取 {len(df_NCP)} 行数据")        
#         print("\n步骤2: 移除完全包含的行...")
#         df_remove = remove_fully_contained(df_NCP)
#         print(f"过滤后剩余 {len(df_remove)} 行数据")        
#         print("\n步骤3: 处理junction区域（去除CDS重叠）...")
#         df_r_CDS = process_junctions(df_remove, df_gff)        
#         print("\n步骤4: 处理重叠区域...")
#         df_filtered = filter_overlapping_regions(df_r_CDS)
#         print(f"重叠处理后剩余 {len(df_filtered)} 行数据")
#         print("\n步骤5: 生成gff文件")
#         excel_to_gff3(df_filtered, output_gff)
#         df_filtered.to_excel("/Volumes/caca/test_fractionation/00raw/sp_loc/gff_results_test.xlsx", index=False)
#     except Exception as e:
#         print(f"\n处理过程中出错: {e}")
# if __name__ == "__main__":
#     main()

# 筛选表达的小肽
import pandas as pd
import os
from tqdm import tqdm
def merge_count_files(input_dir, RNA_info_file, output_file, gene_id_col="Geneid"):
    count_files = pd.read_excel(RNA_info_file, sheet_name="group")
    merged_df = None
    for _, row in tqdm(count_files.iterrows(), desc="合并进度"):
        file = row['Sample']
        sample_name = f"{file}_counts.txt"
        file_path = os.path.join(input_dir, sample_name)
        try:
            df = pd.read_csv(file_path, sep='\t', comment='#')
            counts = df[[gene_id_col, df.columns[-1]]]
            counts.columns = ['GeneID', file]
            if merged_df is None:
                merged_df = counts
            else:
                merged_df = pd.merge(merged_df, counts, on='GeneID', how='outer')
        except Exception as e:
            print(f"\n处理失败 {file}: {str(e)}")
            continue
    if merged_df is not None:
        print(f"\n合并后数据维度: {merged_df.shape}")
        if merged_df.duplicated('GeneID').any():
            print(f"警告：存在重复基因ID，将取第一个出现的值")
            merged_df = merged_df.drop_duplicates('GeneID')
        merged_df.to_csv(output_file, index=False)
        return merged_df
    else:
        raise ValueError(f"错误：未成功合并任何数据")
def combine_gene_sp_data(gene_df, sp_df, output_file):
    combined_df = pd.concat([gene_df, sp_df], axis=0)
    combined_df.to_csv(output_file, index=False)
    print(f"\n合并后总数据维度: {combined_df.shape}")
    return combined_df
def compute_cpm(count_df, output_file):
    count_df = count_df.set_index(count_df.columns[0]) 
    count_df = count_df.apply(pd.to_numeric, errors='coerce')
    library_sizes = count_df.sum(axis=0)
    cpm = count_df.div(library_sizes, axis=1) * 1e6
    cpm.to_csv(output_file, index=True)
    return cpm
def filter_expressed_genes(count_df_cpm, RNA_info_file, output_prefix):
    group_df = pd.read_excel(RNA_info_file, sheet_name="group")
    srr_col = group_df.columns[0]
    group_col = group_df.columns[2]
    group_df = group_df[group_df[srr_col].isin(count_df_cpm.columns)]
    group_to_srrs = (group_df.groupby(group_col)[srr_col].apply(list).to_dict())
    groupwise_cv = {}
    for group, srr_list in group_to_srrs.items():
        group_cpm = count_df_cpm[srr_list]
        group_cv = group_cpm.std(axis=1) / group_cpm.mean(axis=1)
        group_cv[group_cpm.mean(axis=1) == 0] = 0
        groupwise_cv[group] = group_cv
    cv_df = pd.DataFrame(groupwise_cv)
    cv_df.index = count_df_cpm.index
    cv_df.to_csv(f"{output_prefix}_groupwise_cv.csv")
    # 筛选开始
    valid_genes_by_cv = []
    for gene in cv_df.index:
        gene_cv = cv_df.loc[gene]
        valid_groups = (gene_cv <= 0.5).sum()
        total_groups = len(gene_cv)
        if valid_groups / total_groups >= 0.95:
            valid_genes_by_cv.append(gene)
    filtered_cpm = count_df_cpm.loc[valid_genes_by_cv]
    valid_genes = []
    for group, srr_list in group_to_srrs.items():
        group_cpm = filtered_cpm[srr_list]
        valid_genes_in_group = group_cpm[(group_cpm >= 1).all(axis=1)].index
        valid_genes.extend(valid_genes_in_group)
    valid_genes = list(set(valid_genes))
    return valid_genes, filtered_cpm
def extracted_sp_info(valid_genes, filter_cpm, sp_info_file, output_file, output_prefix):
    sp_info_df = pd.read_excel(sp_info_file, sheet_name="NCP")
    filtered_sp_info = sp_info_df[sp_info_df['ID'].isin(valid_genes)].copy()
    filtered_sp_info['proteins'] = pd.to_numeric(filtered_sp_info['proteins'], errors='coerce')
    filtered_sp_info = filtered_sp_info[filtered_sp_info['proteins']==1]
    filtered_sp_info.to_excel(output_file, index=False)
    final_sp_ids = filtered_sp_info['ID'].tolist()
    final_filtered_genes = filter_cpm.loc[final_sp_ids] 
    final_filtered_genes.to_csv(f"{output_prefix}_filtered_genes.csv")
    print(f"\n提取表达小肽信息完成，结果保存至: {output_file}")
def main():
    base_dir = r"D:\Desktop\peptidemicro\00file\00raw\rnaseq"
    gene_input_dir = os.path.join(base_dir, "02count_gene")
    sp_input_dir = os.path.join(base_dir, "02count_sp")
    RNA_info_file = os.path.join(base_dir, "Total_rna_seq.xlsx")
    output_dir = os.path.join(base_dir, "03ouput")
    sp_info_file = r"D:\Desktop\peptidemicro\00file\00raw\sp_loc\Eu_sp_finally.xlsx"
    os.makedirs(output_dir, exist_ok=True)
    print("=== 步骤1/4: 合并计数文件 ===")
    gene_matrix = merge_count_files(
        gene_input_dir, RNA_info_file,
        output_file=os.path.join(output_dir, "total_gene_matrix.csv")
    )
    sp_matrix = merge_count_files(
        sp_input_dir, RNA_info_file,
        output_file=os.path.join(output_dir, "total_sp_matrix_all.csv")
    )
    print("\n=== 步骤2/4: 合并基因/小肽 ===")
    combined_matrix = combine_gene_sp_data(
        gene_matrix, sp_matrix,
        os.path.join(output_dir, "total_combined_matrix.csv")
    )
    combined_matrix_cpm = compute_cpm(combined_matrix, os.path.join(output_dir, "total_combined_matrix_cpm.csv"))
    print("\n=== 步骤3/4: 筛选表达基因/小肽 ===")
    valid_genes, filter_cpm = filter_expressed_genes(
        combined_matrix_cpm, RNA_info_file, os.path.join(output_dir, "finally")
    )
    print("\n=== 步骤4/4: 提取表达小肽信息 ===")
    extracted_sp_info(valid_genes, filter_cpm, sp_info_file, os.path.join(output_dir, "finally_expressed_sp_info.xlsx"), os.path.join(output_dir, "finally"))
if __name__ == "__main__":
    main()

