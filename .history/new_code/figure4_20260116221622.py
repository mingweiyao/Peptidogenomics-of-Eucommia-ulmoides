# # 将Excel中的基因坐标转换为GFF3格式：后续GFF文件生成以此代码为参考
# from datetime import datetime
# import pandas as pd
# def excel_to_gff3(df, output_gff):
#     required_columns = ['ID', 'start', 'end', 'strand', 'chrom']
#     missing_cols = [col for col in required_columns if col not in df.columns]
#     if missing_cols:
#         print(f"错误: 缺少必要的列: {', '.join(missing_cols)}")
#         return
#     gff_header = f"""##gff-version 3
# ##date {datetime.now().strftime('%Y-%m-%d')}
# ##source EuNCP
# ##genome-build v1.0
# """
#     gff_lines = []
#     used_ids = set()
#     for idx, row in df.iterrows():
#         seqid = str(row['chrom']).strip()
#         source = "EuNCP"
#         start = int(row['start'])
#         end = int(row['end'])
#         if start > end:
#             start, end = end, start
#         strand = str(row['strand']).strip()
#         if strand not in ['+', '-']:
#             strand = '.'
#         gene_id = str(row['ID']).strip()
#         if gene_id == "" or gene_id.lower() == "nan":
#             continue
#         base_gene_id = gene_id
#         n = 1
#         while gene_id in used_ids:
#             n += 1
#             gene_id = f"{base_gene_id}_{n}"
#         used_ids.add(gene_id)
#         gene_line = (
#             f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={gene_id};Name={gene_id}"
#         )
#         gff_lines.append(gene_line)
#         mrna_id = f"{gene_id}.t1"
#         mrna_line = (
#             f"{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
#             f"ID={mrna_id};Parent={gene_id};product=predicted_protein"
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
#         with open(output_gff, 'w', encoding='utf-8') as f:
#             f.write(gff_header)
#             f.write("\n".join(gff_lines))
#             f.write("\n")  # ✅ 修正：文件末尾补换行更稳
#         print(f"成功生成GFF3文件: {output_gff}")
#         print(f"转换了 {len(df)} 条记录（空ID已跳过），共生成 {len(gff_lines)} 行GFF记录")
#     except Exception as e:
#         print(f"写入GFF文件失败: {e}")
# excel_file = r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_info.xlsx"
# output_gff = r"F:\work_mechanism\new_file\02figure\figure4\sp_codon.gff"
# df = pd.read_excel(excel_file, sheet_name="unique")
# excel_to_gff3(df, output_gff)

# # 转录组数据基础的NCP分布veen图
# import pandas as pd
# gene_expression_df = pd.read_csv(r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_cpm.csv", index_col=0)
# tissue_mapping_df = pd.read_excel(r"D:\Desktop\Total_rna_seq.xlsx", sheet_name="extract")
# tissue_group_to_samples = (tissue_mapping_df.groupby(['Tissues', 'Group'])['Sample'].apply(list).to_dict())
# peptide_ids_by_tissue = {}
# for (tissue, group), samples in tissue_group_to_samples.items():
#     samples = [s for s in samples if s in gene_expression_df.columns]
#     group_expr = gene_expression_df.loc[:, samples]
#     group_mask = (group_expr >= 1).all(axis=1)
#     if tissue not in peptide_ids_by_tissue:
#         peptide_ids_by_tissue[tissue] = set()
#     peptide_ids_by_tissue[tissue].update(group_expr.index[group_mask].tolist())
# peptide_ids_by_tissue = {tissue: list(ids) for tissue, ids in peptide_ids_by_tissue.items()}
# max_length = max((len(ids) for ids in peptide_ids_by_tissue.values()), default=0)
# peptide_ids_df = pd.DataFrame({
#     tissue: ids + [None] * (max_length - len(ids))
#     for tissue, ids in peptide_ids_by_tissue.items()
# })
# output_excel_path = r"D:\Desktop\rnaseq_veen_leaf_location.xlsx"
# peptide_ids_df.to_excel(output_excel_path, index=False)

# # 数量-肽统计
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# from brokenaxes import brokenaxes
# import warnings
# warnings.filterwarnings('ignore')
# def plot_peptide_distribution_bar_broken_y(expression_matrix_path, output_pdf_path, output_csv_path, threshold=1, ylims=((0, 400), (600, 800), (2800, 3000)), bar_color="#4C72B0", edge_color="white"):
#     if expression_matrix_path.endswith('.csv'):
#         df = pd.read_csv(expression_matrix_path, index_col=0)
#     elif expression_matrix_path.endswith(('.xls', '.xlsx')):
#         df = pd.read_excel(expression_matrix_path, index_col=0)
#     else:
#         raise ValueError("仅支持 csv / xlsx")
#     print(f"读取矩阵：{df.shape[0]} 个ID × {df.shape[1]} 个列")
#     # ---------- 2. 二值化 & 统计 ----------
#     binary_df = (df >= threshold).astype(int)
#     counts_per_id = binary_df.sum(axis=1)
#     count_dist = counts_per_id.value_counts().sort_index()
#     total_ids = len(counts_per_id)
#     plot_data = pd.DataFrame({
#         "Number_of_Columns": count_dist.index.astype(int),
#         "Number_of_IDs": count_dist.values.astype(int),
#         "Percentage": (count_dist.values / total_ids * 100).round(2)
#     })
#     # ---------- 3) 输出 CSV ----------
#     plot_data.to_csv(output_csv_path, index=False)
#     print(f"✅ 统计结果已保存为：{output_csv_path}")
#     # ---------- 4) 画柱状图（broken y-axis） ----------
#     fig = plt.figure(figsize=(10, 8))
#     bax = brokenaxes(ylims=ylims, hspace=0.04, fig=fig)
#     x = plot_data["Number_of_Columns"].values
#     y = plot_data["Number_of_IDs"].values
#     bax.bar(x, y, width=0.8, color=bar_color, edgecolor=edge_color, linewidth=1.0)
#     bax.set_xlabel(f"Number of columns with CPM ≥ {threshold}")
#     bax.set_ylabel("Count of IDs")
#     bax.set_title("Peptide/ID Expression Distribution (Y-axis truncated)")
#     fig.savefig(output_pdf_path, bbox_inches="tight")
#     plt.close(fig)
#     print(f"✅ 截断柱状图已保存为：{output_pdf_path}")
#     return plot_data
# if __name__ == "__main__":
#     expression_file = r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_cpm.csv"
#     output_pdf = r"F:\work_mechanism\new_file\02figure\figure4\test_peptide_distribution_bar_truncated.pdf"
#     output_csv = r"F:\work_mechanism\new_file\02figure\figure4\test_peptide_distribution_stat.csv"
#     plot_peptide_distribution_bar_broken_y(
#         expression_file,
#         output_pdf,
#         output_csv,
#         threshold=1,
#         ylims=((0, 400), (600, 800), (2800, 3000))
#     )

# # 筛选有肽段支持的新转录本
# import os
# import pandas as pd
# import gffutils
# def build_or_load_db(anno_file, db_path="transcriptome.db"):
#     if os.path.exists(db_path):
#         return gffutils.FeatureDB(db_path, keep_order=True)
#     db = gffutils.create_db(
#         anno_file,
#         dbfn=db_path,
#         force=True,
#         keep_order=True,
#         merge_strategy="merge",
#         sort_attribute_values=True,
#         disable_infer_transcripts=True,
#         disable_infer_genes=True
#     )
#     return db
# def gtf_attr_str(attrs: dict):
#     parts = []
#     for k, vs in attrs.items():
#         if vs is None: continue
#         if not isinstance(vs, (list, tuple)): vs = [vs]
#         for v in vs: parts.append(f'{k} "{v}";')
#     return " ".join(parts)
# def feature_to_gtf_line(f):
#     seqid = f.seqid
#     source = f.source if f.source is not None else "."
#     featuretype = f.featuretype
#     start = int(f.start)
#     end = int(f.end)
#     score = f.score if f.score is not None else "."
#     strand = f.strand if f.strand is not None else "."
#     frame = f.frame if f.frame is not None else "."
#     attr = gtf_attr_str(f.attributes)
#     return "\t".join([seqid, source, featuretype, str(start), str(end), str(score), strand, str(frame), attr])
# def overlap(a_start, a_end, b_start, b_end):
#     return a_start <= b_start <= b_end <= a_end
# def annotate_peptides_and_export(peptide_df, db, out_peptide_xlsx, out_hit_gtf):
#     hit_tx_ids = set()
#     add_cols = {
#         "hit_transcript_ids": [],
#         "hit_tx_start": [],
#         "hit_tx_end": [],
#         "hit_tx_strand": [],
#     }
#     for _, row in peptide_df.iterrows():
#         chrom = str(row["chrom"])
#         pep_start = int(row["start"])
#         pep_end = int(row["end"])
#         pep_strand = str(row['strand'])
#         not_filter_transcripts = db.region(seqid=chrom, start=pep_start, end=pep_end, featuretype="transcript")
#         transcripts = [tx for tx in not_filter_transcripts if tx.strand == pep_strand or tx.strand == '.']
#         cur_ids, cur_starts, cur_ends, cur_strands = [], [], [], []
#         for tx in transcripts:
#             if not overlap(tx.start, tx.end, pep_start, pep_end): continue
#             exon_hit = False
#             for exon in db.children(tx, featuretype="exon", order_by="start"):
#                 if overlap(exon.start, exon.end, pep_start, pep_end):
#                     exon_hit = True
#                     break
#             if not exon_hit: continue
#             # 命中这个 transcript
#             hit_tx_ids.add(tx.id)
#             cur_ids.append(tx.id)
#             cur_starts.append(str(int(tx.start)))
#             cur_ends.append(str(int(tx.end)))
#             cur_strands.append(tx.strand if tx.strand is not None else ".")
#         add_cols["hit_transcript_ids"].append(";".join(cur_ids) if cur_ids else "")
#         add_cols["hit_tx_start"].append(";".join(cur_starts) if cur_starts else "")
#         add_cols["hit_tx_end"].append(";".join(cur_ends) if cur_ends else "")
#         add_cols["hit_tx_strand"].append(";".join(cur_strands) if cur_strands else "")
#     # 1) 输出：肽段表追加命中信息
#     out_df = peptide_df.copy()
#     for k, v in add_cols.items(): out_df[k] = v
#     out_df.to_excel(out_peptide_xlsx, index=False)
#     # 2) 输出：命中的转录本 + 全部 exon 到新的 GTF
#     with open(out_hit_gtf, "w") as fw:
#         fw.write("# subset GTF: transcripts validated by peptides (exon-overlap)\n")
#         for tx_id in sorted(hit_tx_ids):
#             tx = db[tx_id]
#             fw.write(feature_to_gtf_line(tx) + "\n")
#             for exon in db.children(tx, featuretype="exon", order_by="start"):
#                 fw.write(feature_to_gtf_line(exon) + "\n")
# def main():
#     excel_file = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\finally_expressed_sp_info.xlsx"
#     sheet_name = "unique"
#     anno_file = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\merged.gtf"
#     out_peptide_xlsx = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\peptide_annotated.xlsx"
#     out_hit_gtf = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\hit_transcripts.gtf"
#     db_path = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\transcriptome.db"
#     peptide_df = pd.read_excel(excel_file, sheet_name=sheet_name)
#     for col in ["chrom", "start", "end", "strand"]:
#         if col not in peptide_df.columns: raise ValueError(f"Excel 缺少必要列: {col}")
#     db = build_or_load_db(anno_file, db_path=db_path)
#     annotate_peptides_and_export(peptide_df, db, out_peptide_xlsx, out_hit_gtf)
#     print(f"完成：\n- {out_peptide_xlsx}\n- {out_hit_gtf}\n- (db) {db_path}")
# if __name__ == "__main__":
#     main()



# # 提取基因的Kozak序列
# import pandas as pd
# from Bio import SeqIO
# import re
# from tqdm import tqdm
# position_regex = re.compile(r'Position=([GWHBISF\d]+): (\d+).*?(\d+):\s([+-])')
# def extract_translation_start_positions(protein_file):
#     protein_start_positions = {}
#     for rec in SeqIO.parse(protein_file, "fasta"):
#         description = rec.description
#         match = position_regex.search(description)
#         if match:
#             chrom = match.group(1) 
#             start_position = match.group(2)
#             end_position = match.group(3)
#             strand = match.group(4)
#             protein_start_positions[rec.id] = (chrom, int(start_position), int(end_position), strand)
#     return protein_start_positions
# def extract_kozak_sequence_from_genome(genome_dict, chrom, start, end, strand, flank=6):
#     if chrom in genome_dict:
#         seq = genome_dict[chrom]
#         if strand == '-':
#             kozak_seq = seq[end-9:end+6].reverse_complement()
#         else:
#             kozak_seq = seq[start-7:start+8]
#     return kozak_seq
# def extract_kozak_sequences(protein_file, genome_file):
#     protein_positions = extract_translation_start_positions(protein_file)
#     kozak_sequences = {}
#     genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(genome_file, "fasta")}
#     for protein_id, (gene_id, start, end, strand) in tqdm(protein_positions.items(), desc="Extracting Kozak sequences", unit="protein"):
#         chrom = gene_id
#         kozak_seq = extract_kozak_sequence_from_genome(genome_dict, chrom, start, end, strand)
#         if kozak_seq:
#             kozak_sequences[protein_id] = kozak_seq
#         else:
#             print(f"警告: 未找到基因组序列 {protein_id} 对应的 Kozak 序列")
#     return kozak_sequences
# def main():
#     protein_file = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_peptide_database.fa"
#     genome_file = r"D:\Desktop\peptidemicro\00file\00raw\Raw_database\Eu_genome.fasta"
#     output_file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\gene_kozak_sequences.xlsx"
#     kozak_sequences = extract_kozak_sequences(protein_file, genome_file)
#     df_kozak = pd.DataFrame(list(kozak_sequences.items()), columns=['Protein ID', 'Kozak Sequence'])
#     df_kozak.to_excel(output_file, index=False)
# if __name__ == "__main__":
#     main()

# # 起始密码子预测
# import pandas as pd
# from Bio import SeqIO
# from multiprocessing import Pool
# from Bio.Seq import Seq
# import os, gffutils
# from Bio.SeqRecord import SeqRecord
# from collections import defaultdict
# # -----------------------------
# # 1) 常量与全局变量
# # -----------------------------
# START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
# STOP_CODONS  = {"TAA", "TAG", "TGA"}
# MINUS_START_CODONS = {"CAT", "CAG", "CAC", "CAA", "CGT"}
# MINUS_STOP_CODONS  = {"TTA", "CTA", "TCA"}
# MAX_SCAN_NT = 600
# THREADS = 10
# _max_scan_nt = None
# _hit_transcript_fasta = None
# # ========== 1. 汇总 accession 覆盖信息 ===========
# def extract_fasta_from_gtf(hit_transcript, db_path, genome_fasta):
#     if not os.path.exists(db_path):
#         gffutils.create_db(hit_transcript, dbfn=db_path, force=True,
#             keep_order=True, merge_strategy="merge", sort_attribute_values=True,
#             disable_infer_transcripts=True, disable_infer_genes=True)
#     db = gffutils.FeatureDB(db_path, keep_order=True)
#     genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(genome_fasta, "fasta")}
#     tx_iter = db.features_of_type("transcript")
#     fasta_dict = {}
#     exons_coords_dict = {}
#     hit_transcript_fasta_output = []
#     hit_transcript_fasta_output_extra = []
    # for tx in tx_iter:
    #     tid = tx.id
    #     tid_start = int(tx.start)
    #     tid_end = int(tx.end)
    #     chrom = tx.chrom
    #     exons = list(db.children(tx, featuretype="exon", order_by="start"))
    #     exons_coords = [(e.start, e.end) for e in exons]
    #     exons_coords.sort(key=lambda x: x[0])
    #     exons_coords_dict[tid] = exons_coords
    #     parts = [genome_dict[chrom][s - 1:e] for (s, e) in exons_coords]
    #     spliced = Seq("").join(parts)
    #     fasta_dict[tid] = str(spliced)
    #     L = len(genome_dict[chrom])
    #     left0  = max(0, tid_start-1-200)
    #     left1  = max(0, tid_start-1)
    #     right0 = min(L, tid_end)
    #     right1 = min(L, tid_end+200)
    #     left_flank  = genome_dict[chrom][left0:left1]
    #     right_flank = genome_dict[chrom][right0:right1]
    #     spliced_extra = Seq("").join([left_flank] + parts + [right_flank])
    #     rec = SeqRecord(Seq(spliced), id = tid, description="")
    #     rec_extra = SeqRecord(Seq(spliced_extra), id = tid, description="")
    #     hit_transcript_fasta_output.append(rec)
    #     hit_transcript_fasta_output_extra.append(rec_extra)
#     SeqIO.write(hit_transcript_fasta_output, "/Volumes/caca/work_mechanism/new_file/02figure/figure4/new_transcript/hit_transcripts.fa", "fasta")
#     SeqIO.write(hit_transcript_fasta_output_extra, "/Volumes/caca/work_mechanism/new_file/02figure/figure4/new_transcript/hit_transcripts_extra.fa", "fasta")
#     return fasta_dict, exons_coords_dict
# def build_transcript_mapper(exons_coords):
#     cum = 0
#     exon_blocks = []
#     for (s,e) in exons_coords:
#         exon_blocks.append((s,e,cum))
#         cum += (e-s+1)
#     def map_genomic_pos_to_tpos(gpos):
#         for (s, e, cum_before) in exon_blocks:
#             if s <= gpos <= e:
#                 offset = gpos - s
#                 return cum_before + offset + 1
#         return None
#     return map_genomic_pos_to_tpos        
# def filter_peptide_seq_cal_cov(peptide_df, database_dict, exons_coords_dict):
#     accession_segments = {}
#     for _, row in peptide_df.iterrows():
#         hit_trans_id = row["hit_transcript_ids"]
#         hit_trans_ids = hit_trans_id.split(";")
    #     start = int(row['start'])
    #     end = int(row['end'])
    #     strand = row['strand']
    #     chrom = row['chrom']
    #     for id in hit_trans_ids:
    #         mapper = build_transcript_mapper(exons_coords_dict[id])
    #         trans_start = mapper(start)
    #         trans_end = mapper(end)
    #         accession_segments.setdefault(id, []).append((trans_start, trans_end, strand, chrom))
    # stats = []
    # for accession, segments in accession_segments.items():
    #     if accession not in database_dict: continue
    #     frame_groups = defaultdict(list)
    #     for (ts, te, strand, chrom) in segments:
    #         frame = (ts - 1) % 3
    #         frame_groups[frame].append((ts, te, strand, chrom))
    #     for frame, segs in frame_groups.items():
    #         min_start = min(s[0] for s in segs)
    #         max_end   = max(s[1] for s in segs)
    #         strand = segs[0][2]
    #         chrom  = segs[0][3]
    #         peptide_count = len(segs)
    #         stats.append({
    #             'trans_id': accession,
#                 'trans_id_frame': f"{accession}|frame{frame}",  # 推荐：用于输出区分
#                 'frame': frame,
#                 'strand': strand,
#                 'chrom': chrom,
#                 'min_start': min_start,
#                 'max_end': max_end,
#                 'peptide_count': peptide_count
#             })
#     return stats
# # ========== 2. 多进程扫描并枚举候选 ==========
# def init_worker(hit_transcript_fasta, max_scan_nt):
#     global _max_scan_nt, _hit_transcript_fasta
#     _max_scan_nt = max_scan_nt
#     _hit_transcript_fasta = hit_transcript_fasta
# def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank=6):
#     genome_seq = Seq(str(genome_seq).upper())
#     if strand == '+': ctx = genome_seq[phy_start - 1 - flank: phy_start + flank + 2]
#     else: ctx = genome_seq[phy_end - 3 - flank: phy_end + flank].reverse_complement()
#     return {"context": str(ctx)}
# def _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         s = min_start - 3 * i
#         if s - 1 < 0 or s + 2 > L: continue
#         triplet = seq_str[s - 1:s + 2]
#         if triplet in START_CODONS: yield s, triplet
# def _find_stop_plus(seq_str, max_end, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for j in range(max_steps):
#         e = max_end + 3 * j
#         if e + 3 > L: break
#         triplet = seq_str[e:e + 3]
#         if triplet in STOP_CODONS: return e
#     return None
# def enumerate_orf_candidates_plus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
#     seq_str = str(Seq(str(genome_seq)).upper())
#     candidates = []
#     for phy_start, triplet in _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
#         phy_end = _find_stop_plus(seq_str, max_end, max_scan_nt)
#         if phy_end is None: continue
#         if not (phy_start <= min_start and phy_end >= max_end): continue
#         kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='+', flank=flank)
#         candidates.append({
#             "phy_start": phy_start,
#             "phy_end": phy_end,
#             "codon": triplet,
#             "kozak_seq": kozak["context"],
#             "start_to_peptide_nt": int(min_start - phy_start)
#         })
#     return candidates
# def _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         left = max_end + 3 * (i - 1)
#         right = max_end + 3 * i
#         if right <= 0 or left >= L: continue
#         triplet = seq_str[left:right]
#         if len(triplet) != 3: continue
#         if triplet in MINUS_START_CODONS:
#             phy_end = max_end + 3 * i
#             yield phy_end, triplet
# def _find_stop_minus(seq_str, min_start, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         left = (min_start - 1) - 3 * (i + 1)
#         right = (min_start - 1) - 3 * i
#         if right <= 0 or left >= L: continue
#         triplet = seq_str[left:right]
#         if len(triplet) != 3: continue
#         phy_start = min_start - 3 * i
#         if triplet in MINUS_STOP_CODONS: return phy_start
#     return None
# def enumerate_orf_candidates_minus(genome_seq, min_start, max_end, max_scan_nt=MAX_SCAN_NT, flank=6):
#     seq_str = str(Seq(str(genome_seq)).upper())
#     candidates = []
#     for phy_end, triplet_raw in _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
#         phy_start = _find_stop_minus(seq_str, min_start, max_scan_nt)
#         if phy_start is None: continue
#         if not (phy_start <= min_start and phy_end >= max_end): continue
#         kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='-', flank=flank)
#         triplet_rc = str(Seq(triplet_raw).reverse_complement())
#         candidates.append({
#             "phy_start": phy_start,
#             "phy_end": phy_end,
#             "codon": triplet_rc,
#             "kozak_seq": kozak["context"],
#             "start_to_peptide_nt": int(phy_end - max_end)
#         })
#     return candidates
# def run_scan_and_output_for_item(item):
#     global _max_scan_nt
#     strand = item['strand']
#     min_start = int(item['min_start'])
#     max_end = int(item['max_end'])
#     hit_id = item['trans_id']
#     gseq = _hit_transcript_fasta[hit_id]
#     if strand == '+':
#         candidates = enumerate_orf_candidates_plus(gseq, min_start, max_end, _max_scan_nt, flank=6)
#     else:
#         candidates = enumerate_orf_candidates_minus(gseq, min_start, max_end, _max_scan_nt, flank=6)
#     item['candidates'] = candidates
#     return item
# def run_scan_and_output(stats, hit_transcript_fasta, max_scan_nt=MAX_SCAN_NT, nproc=THREADS):
#     with Pool(processes=nproc, initializer=init_worker, initargs=(hit_transcript_fasta, max_scan_nt)) as pool:
#         stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
#     return stats_update
# def main():
#     peptide_info = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/new_transcript/finally_expressed_sp_info_annotated.xlsx"
#     hit_transcript_gtf = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/new_transcript/hit_transcripts.gtf"
#     genome_fasta = "/Volumes/caca/work_mechanism/new_file/02figure/Eu_genome_modified/Eu_genome.fasta"
#     out_cand = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/new_transcript/hit_transcript_predict_orf.csv"
#     db_path = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/new_transcript/hit_transcripts.db"
#     df = pd.read_excel(peptide_info, sheet_name="annotated")
#     # 1) 汇总 accession 覆盖信息
#     hit_transcript_fasta_dict, exons_coords_dict = extract_fasta_from_gtf(hit_transcript_gtf, db_path, genome_fasta)
#     stats = filter_peptide_seq_cal_cov(df, hit_transcript_fasta_dict, exons_coords_dict)
#     # 2) 多进程扫描并枚举候选
#     stats_update = run_scan_and_output(stats, hit_transcript_fasta_dict, max_scan_nt=MAX_SCAN_NT, nproc=THREADS)
#     # 3) 输出 candidates 表
#     cand_rows = []
#     for it in stats_update:
#         base = {k: it[k] for k in it.keys() if k not in ("candidates")}
#         for rank, cand in enumerate(it.get("candidates", []), start=1):
#             cand_rows.append({**base, "rank": rank, **cand})
#     df_cand = pd.DataFrame(cand_rows)
#     df_cand.to_csv(out_cand, index=False)  
# if __name__ == "__main__":
#     main()

# # 根据筛选的excel文件提取fasta和dna序列
# from Bio.SeqRecord import SeqRecord
# from Bio import SeqIO
# import pandas as pd
# from Bio.Seq import Seq
# excel_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/codon/codon_prediction/codon_prediction_v7/candidates_scored.xlsx"
# transcript_fa = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/new_transcript/hit_transcripts.fa"
# candidate_trans_dna = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/codon/codon_prediction/codon_prediction_v7/candidate_trans_dna.fa"
# candidate_trans_pro = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/codon/codon_prediction/codon_prediction_v7/candidate_trans_pro.fa"
# transcript_dict = {rec.id: rec.seq for rec in SeqIO.parse(transcript_fa, "fasta")}
# df = pd.read_excel(excel_file, sheet_name="finally")
# trans_dna = []
# trans_pro = []
# for _, row in df.iterrows():
#     candidate_key = row['candidate_key']
#     trans_id = row['trans_id']
#     start = row['phy_start']
#     end = row['phy_end']
#     strand = row['strand']
#     gseq = transcript_dict[trans_id]
#     target_seq_dna = gseq[start-1:end] if strand == "+" else Seq(gseq[start-1:end]).reverse_complement()
#     target_seq_pro = target_seq_dna.translate()
#     rec_dna = SeqRecord(Seq(target_seq_dna), id=candidate_key, description="")
#     rec_pro = SeqRecord(Seq(target_seq_pro), id=candidate_key, description="")
#     trans_dna.append(rec_dna)
#     trans_pro.append(rec_pro)
# SeqIO.write(trans_dna, candidate_trans_dna, "fasta")
# SeqIO.write(trans_pro, candidate_trans_pro, "fasta")

# ==================== 画图 =====================
# # figure1: pwm矩阵文件
# import json
# import numpy as np
# import pandas as pd
# def load_pwm(path):
#     with open(path, "r") as f: pwm = json.load(f)
#     return pd.DataFrame(pwm, index=[f"pos{i+1}" for i in range(13)]).T
# def main(tp_json, tn_json, out_xlsx, eps=1e-6):
#     tp = load_pwm(tp_json)
#     tn = load_pwm(tn_json)
#     log_mat = np.log2((tp + eps) / (tn + eps))
#     with pd.ExcelWriter(out_xlsx) as w:
#         tp.to_excel(w, sheet_name="TP")
#         tn.to_excel(w, sheet_name="TN")
#         log_mat.to_excel(w, sheet_name="LOG")
# if __name__ == "__main__":
#     OUT_DIR = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/codon/codon_prediction/codon_prediction_v4"
#     main(
#         f"{OUT_DIR}/pwm_tp.json",
#         f"{OUT_DIR}/pwm_tn.json",
#         f"/Volumes/caca/work_mechanism/new_file/02figure/figure4/figure4c_pwm_log_matrix.xlsx"
    # )

# import pandas as pd
# IN_XLSX  = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/codon/codon_prediction/codon_prediction_v4/candidates_scored.xlsx"
# OUT_XLSX = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/candidates_scored_with_delta.xlsx"
# df = pd.read_excel(IN_XLSX, sheet_name="Sheet1")
# delta = (
#     df[df["tis_rank"].isin([1, 2])]
#     .pivot(index="accession", columns="tis_rank", values="tis_scores")
# )
# delta["rank1_rank2_delta"] = delta[1] - delta[2]
# df = df.merge(
#     delta["rank1_rank2_delta"],
#     left_on="accession",
#     right_index=True,
#     how="left"
# )
# df.loc[df["tis_rank"] != 1, "rank1_rank2_delta"] = pd.NA
# df.to_excel(OUT_XLSX, index=False)

# # Kozak序列位置分析
# import pandas as pd
# file = r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\output_best.xlsx"
# df = pd.read_excel(file, sheet_name="GTG")
# seqs = df['kozak_seq'].dropna().astype(str)
# seqs = seqs[seqs.str.len() == 15]
# pos_counts = {i: {} for i in range(1, 16)}
# for s in seqs:
#     for i, ch in enumerate(s, start=1):
#         pos_counts[i][ch] = pos_counts[i].get(ch, 0) + 1
# all_chars = sorted({ch for counts in pos_counts.values() for ch in counts.keys()})
# result_df = pd.DataFrame(
#     index=all_chars,
#     columns=[f"X{i}" for i in range(1, 16)]
# )
# for pos in range(1, 16):
#     for ch, cnt in pos_counts[pos].items():
#         result_df.loc[ch, f"X{pos}"] = cnt
# result_df = result_df.fillna(0).astype(int)
# result_df.to_excel(r"D:\Desktop\peptidemicro\00file\01figure\figure4\codon\sp_GTG_kozak_seq_statistic.xlsx")

# import pandas as pd
# in_file = r"F:\work_mechanism\new_file\02figure\figure4\codon\codon_prediction\codon_prediction_v4\candidates_scored.xlsx"
# out_file = r"F:\work_mechanism\new_file\02figure\figure4\kozak_count.xlsx"
# codons = ["CTG", "ACG", "GTG", "ATG", "TTG"]
# SEQ_COL = "kozak_seq"
# CODON_COL = "codon"
# L = 13
# df = pd.read_excel(in_file, sheet_name="rank1")
# def count_matrix(seqs, L=13):
#     bases = ['A', 'C', 'G', 'T']
#     cols = [f"X{i}" for i in range(1, L + 1)]
#     m = pd.DataFrame(0, index=bases, columns=cols, dtype=int)
#     seqs = (pd.Series(seqs).dropna().astype(str).str.upper().str.slice(0, L))
#     seqs = seqs[seqs.str.len() == L]
#     seqs = seqs[seqs.str.fullmatch(rf"[ACGT]{{{L}}}")]
#     for s in seqs:
#         for i, ch in enumerate(s, start=1):
#             m.loc[ch, f"X{i}"] += 1
#     return m
# with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
#     count_matrix(df[SEQ_COL], L=L).to_excel(writer, sheet_name="ALL")
#     codon_series = df[CODON_COL].fillna("").astype(str).str.upper()
#     for c in codons:
#         sub = df.loc[codon_series.str.contains(c, regex=False), SEQ_COL]
#         count_matrix(sub, L=L).to_excel(writer, sheet_name=c)    