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
# _exons_coords_dict = None
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
#     for tx in tx_iter:
#         tid = tx.id
#         tid_start = int(tx.start)
#         tid_end = int(tx.end)
#         chrom = tx.chrom
#         exons = list(db.children(tx, featuretype="exon", order_by="start"))
#         exons_coords = [(e.start, e.end) for e in exons]
#         exons_coords.sort(key=lambda x: x[0])
#         exons_coords_dict[tid] = exons_coords
#         parts = [genome_dict[chrom][s - 1:e] for (s, e) in exons_coords]
#         spliced = Seq("").join(parts)
#         fasta_dict[tid] = str(spliced)
#         L = len(genome_dict[chrom])
#         left0  = max(0, tid_start-1-200)
#         left1  = max(0, tid_start-1)
#         right0 = min(L, tid_end)
#         right1 = min(L, tid_end+200)
#         left_flank  = genome_dict[chrom][left0:left1]
#         right_flank = genome_dict[chrom][right0:right1]
#         spliced_extra = Seq("").join([left_flank] + parts + [right_flank])
#         rec = SeqRecord(Seq(spliced), id = tid, description="")
#         rec_extra = SeqRecord(Seq(spliced_extra), id = tid, description="")
#         hit_transcript_fasta_output.append(rec)
#         hit_transcript_fasta_output_extra.append(rec_extra)
#     SeqIO.write(hit_transcript_fasta_output, "/Volumes/caca/work_mechanism/new_file/02figure/figure5/new_transcript/hit_transcripts.fa", "fasta")
#     SeqIO.write(hit_transcript_fasta_output_extra, "/Volumes/caca/work_mechanism/new_file/02figure/figure5/new_transcript/hit_transcripts_extra.fa", "fasta")
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
#         start = int(row['start'])
#         end = int(row['end'])
#         strand = row['strand']
#         chrom = row['chrom']
#         pep_id = row['ID']
#         for id in hit_trans_ids:
#             mapper = build_transcript_mapper(exons_coords_dict[id])
#             trans_start = mapper(start)
#             trans_end = mapper(end)
#             accession_segments.setdefault(id, []).append((trans_start, trans_end, strand, chrom, pep_id))
#     stats = []
#     for accession, segments in accession_segments.items():
#         if accession not in database_dict: continue
#         frame_groups = defaultdict(list)
#         for (ts, te, strand, chrom, pep_id) in segments:
#             frame = (ts - 1) % 3
#             frame_groups[frame].append((ts, te, strand, chrom, pep_id))
#         for frame, segs in frame_groups.items():
#             peptide_count = len(segs)
#             if peptide_count > 1: continue
#             min_start = segs[0][0]
#             max_end   = segs[0][1]
#             strand = segs[0][2]
#             chrom  = segs[0][3]
#             pep_id = segs[0][4]
#             stats.append({
#                 'trans_id': accession,
#                 'trans_id_frame': f"{accession}|frame{frame}",  # 推荐：用于输出区分
#                 'frame': frame,
#                 'strand': strand,
#                 'chrom': chrom,
#                 'min_start': min_start,
#                 'max_end': max_end,
#                 'pep_id': pep_id
#             })
#     return stats
# # ========== 2. 多进程扫描并枚举候选 ==========
# def init_worker(hit_transcript_fasta, exons_coords_dict, max_scan_nt):
#     global _max_scan_nt, _hit_transcript_fasta, _exons_coords_dict
#     _max_scan_nt = max_scan_nt
#     _hit_transcript_fasta = hit_transcript_fasta
#     _exons_coords_dict = exons_coords_dict
# def map_tpos_interval_to_genome_segments(exons_coords):
#     exons_coords = sorted(exons_coords, key=lambda x: x[0])
#     cum = 0
#     exon_blocks = []
#     for (gs, ge) in exons_coords:
#         L = ge - gs + 1
#         t_start = cum + 1
#         t_end = cum + L
#         exon_blocks.append((t_start, t_end, gs, ge))
#         cum += L
#     def map_tpos_to_gpos(tpos):
#         tpos = int(tpos)
#         for i in range(len(exon_blocks)):
#             if exon_blocks[i][0] <= tpos <= exon_blocks[i][1]:
#                 offset = tpos - exon_blocks[i][0]
#                 return offset + exon_blocks[i][2], i
#         return None
#     return map_tpos_to_gpos
# def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank=6):
#     genome_seq = Seq(str(genome_seq).upper())
#     L = len(genome_seq)
#     if strand == '+':
#         a = max(0, (phy_start - 1 - flank))
#         b = min(L, (phy_start + flank + 2))
#         ctx = genome_seq[a:b]
#     else:
#         a = max(0, (phy_end - 3 - flank))
#         b = min(L, (phy_end + flank))
#         ctx = genome_seq[a:b].reverse_complement()
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
#             "kozak_seq": kozak["context"]
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
#             "kozak_seq": kozak["context"]
#         })
#     return candidates
# def run_scan_and_output_for_item(item):
#     global _max_scan_nt, _exons_coords_dict
#     strand = item['strand']
#     min_start = int(item['min_start'])
#     max_end = int(item['max_end'])
#     hit_id = item['trans_id']
#     gseq = _hit_transcript_fasta[hit_id]
#     if strand == '+':
#         candidates = enumerate_orf_candidates_plus(gseq, min_start, max_end, _max_scan_nt, flank=6)
#     else:
#         candidates = enumerate_orf_candidates_minus(gseq, min_start, max_end, _max_scan_nt, flank=6)
#     exons_coords = _exons_coords_dict.get(hit_id)
#     if exons_coords:
#         mapper_gtf_genome = map_tpos_interval_to_genome_segments(exons_coords)
#         for c in candidates:
#             ps = int(c["phy_start"])
#             pe = int(c["phy_end"])
#             genome_start, genome_start_i = mapper_gtf_genome(ps)
#             genome_end, genome_end_i = mapper_gtf_genome(pe)
#             if genome_start_i == genome_end_i:
#                 c['cross_exon'] = f"{genome_start}-{genome_end}"
#             else:
#                 c['cross_exon'] = ";".join(f"{exons_coords[j][0]}-{exons_coords[j][1]}" for j in range(genome_start_i, genome_end_i + 1))     
#     item['candidates'] = candidates
#     return item
# def run_scan_and_output(stats, hit_transcript_fasta, exons_coords_dict, max_scan_nt=MAX_SCAN_NT, nproc=THREADS):
#     with Pool(processes=nproc, initializer=init_worker, initargs=(hit_transcript_fasta, exons_coords_dict, max_scan_nt)) as pool:
#         stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
#     return stats_update
# def main():
#     peptide_info = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/new_transcript/finally_expressed_sp_info_annotated.xlsx"
#     hit_transcript_gtf = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/new_transcript/hit_transcripts.gtf"
#     genome_fasta = "/Volumes/caca/work_mechanism/new_file/02figure/Eu_genome_modified/Eu_genome.fasta"
#     out_cand = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/new_transcript/hit_transcript_predict_orf.csv"
#     db_path = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/new_transcript/hit_transcripts.db"
#     df = pd.read_excel(peptide_info, sheet_name="annotated")
#     # 1) 汇总 accession 覆盖信息
#     hit_transcript_fasta_dict, exons_coords_dict = extract_fasta_from_gtf(hit_transcript_gtf, db_path, genome_fasta)
#     stats = filter_peptide_seq_cal_cov(df, hit_transcript_fasta_dict, exons_coords_dict)
#     # 2) 多进程扫描并枚举候选
#     stats_update = run_scan_and_output(stats, hit_transcript_fasta_dict, exons_coords_dict, max_scan_nt=MAX_SCAN_NT, nproc=THREADS)
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

# # 根据密码子预测的转录本定量并生成 GFF 文件
# # 将Excel中的基因坐标转换为GFF3格式：后续GFF文件生成以此代码为参考
# from datetime import datetime
# import pandas as pd
# def parse_cross_exon(cross_exon_value):
#     s = str(cross_exon_value).strip()
#     exon_coords = []
#     for block in s.split(";"):
#         block = block.strip()
#         a, b = block.split("-", 1)
#         a, b = int(a), int(b)
#         if a > b: a, b = b, a
#         exon_coords.append((a, b))
#     return exon_coords
# def select_longest_by_trans_id_frame(df):
#     required = ["trans_id_frame", "phy_start", "phy_end"]
#     missing = [c for c in required if c not in df.columns]
#     if missing:
#         raise ValueError(f"缺少用于筛选的列: {', '.join(missing)}")
#     d = df.copy()
#     d["phy_start"] = pd.to_numeric(d["phy_start"], errors="coerce")
#     d["phy_end"] = pd.to_numeric(d["phy_end"], errors="coerce")
#     d = d.dropna(subset=["trans_id_frame", "phy_start", "phy_end"])
#     d["aa_len"] = (d["phy_end"] - d["phy_start"] + 1) / 3
#     idx = d.groupby("trans_id_frame", sort=False)["aa_len"].idxmax()
#     d_best = d.loc[idx].copy()
#     return d_best
# def excel_to_gff3(df, output_gff, source="EuNCP"):
#     required_columns = ['ID', 'strand', 'chrom', 'cross_exon']
#     missing_cols = [col for col in required_columns if col not in df.columns]
#     if missing_cols:
#         print(f"错误: 缺少必要的列: {', '.join(missing_cols)}")
#         return
#     gff_header = f"""##gff-version 3
# ##date {datetime.now().strftime('%Y-%m-%d')}
# ##source {source}
# ##genome-build v1.0
# """
#     gff_lines = []
#     used_ids = set() 
#     for idx, row in df.iterrows():
#         seqid = str(row['chrom']).strip()
#         strand = str(row['strand']).strip()
#         if strand not in ['+', '-']: strand = '.'
#         gene_id = str(row['ID']).strip()
#         if gene_id == "" or gene_id.lower() == "nan": continue
#         exon_coords = parse_cross_exon(row['cross_exon'])
#         start = min(s for s, e in exon_coords)
#         end = max(e for s, e in exon_coords)
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
#         for i, (exon_start, exon_end) in enumerate(exon_coords, start=1):
#             exon_id = f"{mrna_id}.exon{i}"
#             exon_line = (
#                 f"{seqid}\t{source}\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t"
#                 f"ID={exon_id};Parent={mrna_id}"
#             )
#             gff_lines.append(exon_line)
#             cds_id = f"{mrna_id}.cds{i}"
#             cds_line = (
#                 f"{seqid}\t{source}\tCDS\t{exon_start}\t{exon_end}\t.\t{strand}\t0\t"
#                 f"ID={cds_id};Parent={mrna_id}"
#             )
#             gff_lines.append(cds_line)
#     try:
#         with open(output_gff, 'w', encoding='utf-8') as f:
#             f.write(gff_header)
#             f.write("\n".join(gff_lines))
#             f.write("\n")  # ✅ 修正：文件末尾补换行更稳
#         print(f"成功生成GFF3文件: {output_gff}")
#         print(f"转换了 {len(df)} 条记录（空ID已跳过），共生成 {len(gff_lines)} 行GFF记录")
#     except Exception as e:
#         print(f"写入GFF文件失败: {e}")
# excel_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/codon/codon_prediction/codon_prediction_v7/candidates_scored.xlsx"
# output_gff = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/codon/codon_prediction/codon_prediction_v7/EuNCP_trans.gff"
# df = pd.read_excel(excel_file)
# df_best = select_longest_by_trans_id_frame(df)
# excel_to_gff3(df_best, output_gff)

# # 合并文件
# import pandas as pd
# import os
# from tqdm import tqdm
# def merge_count_files(input_dir, RNA_info_file, output_file, gene_id_col="Geneid"):
#     count_files = pd.read_excel(RNA_info_file, sheet_name="condition")
#     merged_df = None
#     for _, row in tqdm(count_files.iterrows(), desc="合并进度"):
#         file = row['Sample']
#         sample_name = f"{file}_counts.txt"
#         file_path = os.path.join(input_dir, sample_name)
#         try:
#             df = pd.read_csv(file_path, sep='\t', comment='#')
#             counts = df[[gene_id_col, df.columns[-1]]]
#             counts.columns = ['GeneID', file]
#             if merged_df is None:
#                 merged_df = counts
#             else:
#                 merged_df = pd.merge(merged_df, counts, on='GeneID', how='outer')
#         except Exception as e:
#             print(f"\n处理失败 {file}: {str(e)}")
#             continue
#     if merged_df is not None:
#         print(f"\n合并后数据维度: {merged_df.shape}")
#         if merged_df.duplicated('GeneID').any():
#             print(f"警告：存在重复基因ID，将取第一个出现的值")
#             merged_df = merged_df.drop_duplicates('GeneID')
#         merged_df.to_csv(output_file, index=False)
#         return merged_df
#     else:
#         raise ValueError(f"错误：未成功合并任何数据")
# def main():
#     count_dir = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rnaseq"
#     RNA_info_file = "/Volumes/caca/work_mechanism/new_file/02figure/Total_rna_seq.xlsx"
#     output_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/total_candidate_sp_count_matrix.csv"
#     merge_count_files(count_dir, RNA_info_file, output_file)
# if __name__ == "__main__":
#     main()

# # 按照FC和Pvalue筛选表达基因
# import pandas as pd
# def filter_excel_by_log2fc_pvalue(input_excel, output_excel, log2fc_col="log2FC", 
#                                   pvalue_col="PValue", log2fc_cutoff=1, pvalue_cutoff=0.05):
#     xls = pd.ExcelFile(input_excel)
#     result_dfs = []
#     for sheet in xls.sheet_names:
#         filtered = []
#         df = pd.read_excel(xls, sheet_name=sheet)
#         if log2fc_col not in df.columns or pvalue_col not in df.columns:
#             print(f"[跳过] sheet '{sheet}' 缺少 {log2fc_col} 或 {pvalue_col}")
#             continue
#         # 转为数值，避免字符串/空值问题
#         df[log2fc_col] = pd.to_numeric(df[log2fc_col], errors="coerce")
#         df[pvalue_col] = pd.to_numeric(df[pvalue_col], errors="coerce")
#         # 筛选条件
#         df_f = df[(df[log2fc_col].abs() > log2fc_cutoff) & (df[pvalue_col] < pvalue_cutoff)].copy()
#         if df_f.empty: continue
#         filtered = df_f[["name", log2fc_col, pvalue_col]].copy()
#         filtered["source_sheet"] = sheet
#         result_dfs.append(filtered)
#     if not result_dfs:
#         print("没有任何 sheet 满足筛选条件")
#         return
#     final_df = pd.concat(result_dfs, ignore_index=True)
#     # 写入新的 Excel
#     final_df.to_excel(output_excel, index=False)
#     print(f"筛选完成，共输出 {len(final_df)} 行")
#     print(f"结果已写入: {output_excel}")
# input_excel = r"F:\work_mechanism\new_file\02figure\figure5\rubber\deseq_analysis.xlsx"
# output_excel = r"F:\work_mechanism\new_file\02figure\figure5\rubber\deseq_analysis_diff.xlsx"
# filter_excel_by_log2fc_pvalue(input_excel, output_excel)

# # 提取相关肽的信息，用于肽分类的构建
# import pandas as pd
# id_file = r"F:\work_mechanism\new_file\02figure\figure5\rubber\deseq_analysis_diff.xlsx"
# temp_file = r"F:\work_mechanism\new_file\02figure\figure5\codon\codon_prediction\codon_prediction_v7\candidates_scored.xlsx"
# expression_info = r"F:\work_mechanism\new_file\01location\rnaseq\03ouput\finally_expressed_sp_info.xlsx"
# df_id = pd.read_excel(id_file, sheet_name="Sheet2")
# temp_file_df = pd.read_excel(temp_file)
# df_express = pd.read_excel(expression_info, sheet_name="unique")
# id_list = df_id['name']
# id_temp = temp_file_df[temp_file_df['ID'].isin(id_list)]
# id_express = df_express[df_express['ID'].isin(id_temp['pep_id'])]
# id_express.to_excel(r"F:\work_mechanism\new_file\02figure\figure5\rubber\deseq_analysis_diff_pep_info.xlsx", index=False)

# # TPM标准化
# import os
# import pandas as pd
# import gffutils
# GFF_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/codon/codon_prediction/codon_prediction_v7/sp_codon.gff"
# COUNT_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/rubber/rubber_count.xlsx"
# OUT_DIR = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/rubber/"
# def prepare_length_data(gff_file):
#     db_path = gff_file + ".db"
#     if not os.path.exists(db_path):
#         print("🔄 正在创建 GFF 数据库...")
#         gffutils.create_db(
#             gff_file,
#             dbfn=db_path,
#             force=True,
#             keep_order=True,
#             merge_strategy="merge",
#             id_spec={"gene": "ID", "mRNA": "ID", "CDS": "Parent"},
#             disable_infer_genes=True,
#             disable_infer_transcripts=True
#         )
#     db = gffutils.FeatureDB(db_path)
#     gene_lengths = {}
#     for gene in db.features_of_type("gene"):
#         total_length = 0
#         for mrna in db.children(gene, featuretype="mRNA", order_by="start"):
#             exons = list(db.children(mrna, featuretype="exon", order_by="start"))
#             if exons:
#                 mrna_len = sum(e.end - e.start + 1 for e in exons)
#                 total_length += mrna_len
#         if total_length > 0:
#             gene_id = gene.id.replace("evm.model.", "evm.TU.")
#             gene_lengths[gene_id] = total_length
#     length_df = pd.DataFrame(
#         gene_lengths.items(),
#         columns=["GeneID", "length"]
#     )
#     out_len = os.path.join(OUT_DIR, "gene_lengths.csv")
#     length_df.to_csv(out_len, index=False)
#     print(f"✅ 基因长度表已生成：{out_len}")
#     return length_df
# def read_counts(count_file):
#     if count_file.lower().endswith(".csv"):
#         df = pd.read_csv(count_file)
#     else:
#         df = pd.read_excel(count_file)
#     if df.columns[0] != "GeneID":
#         df = df.rename(columns={df.columns[0]: "GeneID"})
#     return df
# def normalize_tpm(count_df, length_df, output_file):
#     df = pd.merge(count_df, length_df, on='GeneID', how='inner')
#     df = df[df['length'] > 0]
#     sample_cols = [col for col in df.columns if col not in ['GeneID', 'length']]
#     tpm_data = {}
#     for sample in sample_cols:
#         rpk = (df[sample] * 10**3) / df['length']
#         per_million_scaling_factor = rpk.sum() / 10**6
#         tpm = rpk / per_million_scaling_factor
#         tpm_data[sample] = tpm
#     tpm_df = pd.concat([df[['GeneID']], pd.DataFrame(tpm_data)], axis=1)
#     tpm_df.to_excel(output_file, index=False)
# if __name__ == "__main__":
#     length_df = prepare_length_data(GFF_FILE)
#     count_df = read_counts(COUNT_FILE)
#     normalize_tpm(count_df, length_df, os.path.join(OUT_DIR, "rubber_count_tpm.xlsx"))

# # FPKM标准化
# import os
# import pandas as pd
# import gffutils
# GFF_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/codon/codon_prediction/codon_prediction_v7/sp_codon.gff"
# COUNT_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/rubber_count_5.xlsx"
# OUT_DIR = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/"
# def prepare_length_data(gff_file):
#     db_path = gff_file + ".db"
#     if not os.path.exists(db_path):
#         print("🔄 正在创建 GFF 数据库...")
#         gffutils.create_db(
#             gff_file,
#             dbfn=db_path,
#             force=True,
#             keep_order=True,
#             merge_strategy="merge",
#             id_spec={"gene": "ID", "mRNA": "ID", "CDS": "Parent"},
#             disable_infer_genes=True,
#             disable_infer_transcripts=True
#         )
#     db = gffutils.FeatureDB(db_path)
#     gene_lengths = {}
#     for gene in db.features_of_type("gene"):
#         total_length = 0
#         for mrna in db.children(gene, featuretype="mRNA", order_by="start"):
#             exons = list(db.children(mrna, featuretype="exon", order_by="start"))
#             if exons:
#                 mrna_len = sum(e.end - e.start + 1 for e in exons)
#                 total_length += mrna_len
#         if total_length > 0:
#             gene_id = gene.id.replace("evm.model.", "evm.TU.")
#             gene_lengths[gene_id] = total_length
#     length_df = pd.DataFrame(gene_lengths.items(), columns=["GeneID", "length"])
#     out_len = os.path.join(OUT_DIR, "gene_lengths.csv")
#     length_df.to_csv(out_len, index=False)
#     print(f"✅ 基因长度表已生成：{out_len}")
#     return length_df
# def read_counts(count_file):
#     if count_file.lower().endswith(".csv"):
#         df = pd.read_csv(count_file)
#     else:
#         df = pd.read_excel(count_file)
#     if df.columns[0] != "GeneID":
#         df = df.rename(columns={df.columns[0]: "GeneID"})
#     # 确保 GeneID 为字符串且去空格
#     df["GeneID"] = df["GeneID"].astype(str).str.strip()
#     return df
# def normalize_fpkm(count_df, length_df, output_file):
#     length_df = length_df.copy()
#     length_df["GeneID"] = length_df["GeneID"].astype(str).str.strip()
#     df = pd.merge(count_df, length_df, on="GeneID", how="inner")
#     df = df[df["length"] > 0].copy()
#     sample_cols = [c for c in df.columns if c not in ["GeneID", "length"]]
#     fpkm_data = {}
#     for sample in sample_cols:
#         counts = pd.to_numeric(df[sample], errors="coerce").fillna(0)
#         N = counts.sum()
#         if N == 0:
#             # 全是0时，直接输出0
#             fpkm = counts * 0.0
#         else:
#             L = df["length"].astype(float)  # bp
#             fpkm = (1e9 * counts) / (N * L)
#         fpkm_data[sample] = fpkm
#     fpkm_df = pd.concat([df[["GeneID"]], pd.DataFrame(fpkm_data)], axis=1)
#     fpkm_df.to_excel(output_file, index=False)
#     print(f"✅ FPKM 输出完成：{output_file}")
# if __name__ == "__main__":
#     os.makedirs(OUT_DIR, exist_ok=True)
#     length_df = prepare_length_data(GFF_FILE)
#     count_df = read_counts(COUNT_FILE)
#     normalize_fpkm(count_df, length_df, os.path.join(OUT_DIR, "rubber_count_5_fpkm.xlsx"))

# # 提取特定基因表达量
# import pandas as pd
# id_mapping_df = pd.read_excel("/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/rubber.xlsx", sheet_name="Sheet2")
# data_df = pd.read_excel("/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/rubber_count_5_fpkm.xlsx")
# mapped_ids = id_mapping_df['evm_old']
# mapped_df = data_df[data_df['GeneID'].isin(mapped_ids)]
# mapped_df.to_excel("/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/rubber_count_5_fpkm_filter.xlsx", index=False)

# 按group提取pearson子矩阵
import os
import re
import pandas as pd
PEARSON_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/pearson_r_p.xlsx"
R_SHEET = "R"
P_SHEET = "P"
OUT_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/position_evm_vs_NCPT_RP.xlsx"
EVM_PREFIX = "evm"
NCPT_PREFIX = "NCPT"
MAP_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/rubber.xlsx"
def sanitize_sheetname(name):
    name = str(name)
    name = re.sub(r"[:\\/?*\[\]]+", "_", name).strip()
    if not name: name = "sheet"
    return name[:31]
def is_prefix(x, prefix):
    return str(x).strip().lower().startswith(prefix.lower())
# =========================
# Main
# =========================
def main():
    map_df = pd.read_excel(MAP_FILE, sheet_name="Sheet2", dtype=str)
    map_df.columns = ["evm_old", "evm_new"]
    map_df["evm_old"] = map_df["evm_old"].str.strip()
    map_df["evm_new"] = map_df["evm_new"].str.strip()
    evm_map = dict(zip(map_df["evm_old"], map_df["evm_new"]))
    evm_map_filtered = {k: v for k, v in evm_map.items() if is_prefix(k, EVM_PREFIX)}
    rmat = pd.read_excel(PEARSON_FILE, sheet_name=R_SHEET, index_col=0)
    pmat = pd.read_excel(PEARSON_FILE, sheet_name=P_SHEET, index_col=0)
    rmat = rmat.rename(index=evm_map_filtered, columns=evm_map_filtered)
    pmat = pmat.rename(index=evm_map_filtered, columns=evm_map_filtered)
    evm_new_set = set(map_df["evm_new"])
    evm_ids = [rid for rid in rmat.index if rid in evm_new_set]
    ncpt_cols = [cid for cid in rmat.columns if is_prefix(cid, NCPT_PREFIX)]
    # Prepare output dir
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    used_sheetnames = set()
    with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
        for evm in evm_ids:
            r_row = rmat.loc[evm, ncpt_cols]
            p_row = pmat.loc[evm, ncpt_cols]
            out_df = pd.DataFrame(
                {
                    "var": ncpt_cols,
                    "R": r_row.values,
                    "P": p_row.values,
                }
            )
            out_df["R"] = pd.to_numeric(out_df["R"], errors="coerce")
            out_df["P"] = pd.to_numeric(out_df["P"], errors="coerce")
            sig_df = out_df[(out_df["R"] >= 0.7) & (out_df["P"] <= 0.05)].copy()
            sheet = sanitize_sheetname(evm)
            base = sheet
            k = 1
            while sheet in used_sheetnames:
                suffix = f"_{k}"
                sheet = sanitize_sheetname(base[: (31 - len(suffix))] + suffix)
                k += 1
            used_sheetnames.add(sheet)
            sig_df.to_excel(writer, sheet_name=sheet, index=False)
    print(f"Done. Output written to: {OUT_FILE}")
    print(f"evm sheets: {len(evm_ids)} | NCPT columns used: {len(ncpt_cols)}")
if __name__ == "__main__":
    main()

# # 上个文件的summary
# import pandas as pd
# import os
# IN_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/evm_vs_NCPT_RP.xlsx"
# OUT_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/evm_vs_NCPT_RP_summary.xlsx"
# R_COL = "R"
# def main():
#     xls = pd.ExcelFile(IN_FILE)
#     summary_rows = []
#     for sheet in xls.sheet_names:
#         df = pd.read_excel(xls, sheet_name=sheet)
#         if R_COL not in df.columns: continue
#         r = pd.to_numeric(df[R_COL], errors="coerce")
#         pos_count = (r > 0).sum()
#         neg_count = (r < 0).sum()
#         total_count = r.notna().sum()
#         summary_rows.append({
#             "sheet": sheet,
#             "R_positive(>0)": int(pos_count),
#             "R_negative(<0)": int(neg_count),
#             "R_total(non-NA)": int(total_count),
#         })
#     summary_df = pd.DataFrame(summary_rows)
#     os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
#     summary_df.to_excel(OUT_FILE, index=False)
#     print(f"Done. Summary written to:\n{OUT_FILE}")
# if __name__ == "__main__":
#     main()

# # 上个文件的var
# import pandas as pd
# import os
# IN_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/evm_vs_NCPT_RP.xlsx"
# OUT_FILE = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/evm_vs_NCPT_RP_var.xlsx"
# VAR_COL = "var"
# def main():
#     xls = pd.ExcelFile(IN_FILE)
#     cols = {}
#     for sheet in xls.sheet_names:
#         df = pd.read_excel(xls, sheet_name=sheet)
#         if VAR_COL not in df.columns: continue
#         s = df[VAR_COL].dropna().astype(str).reset_index(drop=True)
#         cols[sheet] = s
#     out_df = pd.DataFrame(cols)
#     out_df = out_df.reindex(sorted(out_df.columns), axis=1)
#     os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
#     out_df.to_excel(OUT_FILE, index=False)
#     print(f"Done. Output written to:\n{OUT_FILE}")
#     print(f"Sheets included: {len(out_df.columns)}")
# if __name__ == "__main__":
#     main()

# # 提取不同sheet的表达量
# import pandas as pd
# import os
# import numpy as np
# in_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/evm_vs_NCPT_RP.xlsx"
# expression_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/rubber_count_5_fpkm_filter.xlsx"
# out_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure5/rubber/correlation/evm_vs_NCPT_RP_fpkm.xlsx"
# # ========= 读取表达矩阵 =========
# df_expr = pd.read_excel(expression_file, sheet_name="Sheet1")
# df_expr["GeneID"] = df_expr["GeneID"].astype(str).str.strip()
# # ========= 读取相关性文件 =========
# xls = pd.ExcelFile(in_file)
# os.makedirs(os.path.dirname(out_file), exist_ok=True)
# with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
#     for sheet in xls.sheet_names:
#         if sheet.upper() == "SUMMARY": continue
#         df = pd.read_excel(xls, sheet_name=sheet)
#         if "var" not in df.columns: continue
#         # 目标基因列表（mubiaoID）
#         target_ids = (df["var"].astype(str).str.strip().dropna().unique())
#         # 从表达矩阵中提取
#         sub_expr = df_expr[df_expr["GeneID"].isin(target_ids)].copy()
#         # 如果你希望表达量顺序和 var 一致（很重要，默认推荐）
#         sub_expr["GeneID"] = pd.Categorical(sub_expr["GeneID"], categories=target_ids, ordered=True)
#         sub_expr = sub_expr.sort_values("GeneID")
#         expr_cols = sub_expr.columns.difference(["GeneID"])
#         sub_expr[expr_cols] = np.log2(sub_expr[expr_cols] + 0.1)        
#         sub_expr.to_excel(writer, sheet_name=sheet[:31], index=False)
# print(f"Done. Expression data written to:\n{out_file}")




























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

# # 根据筛选的excel文件提取fasta和dna序列
# from Bio.SeqRecord import SeqRecord
# from Bio import SeqIO
# import pandas as pd
# from Bio.Seq import Seq
# excel_file = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/codon/codon_prediction/codon_prediction_v7/candidates_scored.xlsx"
# transcript_fa = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/new_transcript/hit_transcripts.fa"
# candidate_trans_dna = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/codon/codon_prediction/codon_prediction_v7/test_candidate_trans_dna.fa"
# candidate_trans_pro = "/Volumes/caca/work_mechanism/new_file/02figure/figure4/codon/codon_prediction/codon_prediction_v7/test_candidate_trans_pro.fa"
# transcript_dict = {rec.id: rec.seq for rec in SeqIO.parse(transcript_fa, "fasta")}
# df = pd.read_excel(excel_file, sheet_name="Sheet1")
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
#     if "*" in target_seq_pro: continue
#     rec_dna = SeqRecord(Seq(target_seq_dna), id=candidate_key, description="")
#     rec_pro = SeqRecord(Seq(target_seq_pro), id=candidate_key, description="")
#     trans_dna.append(rec_dna)
#     trans_pro.append(rec_pro)
# SeqIO.write(trans_dna, candidate_trans_dna, "fasta")
# SeqIO.write(trans_pro, candidate_trans_pro, "fasta")

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

# # TPM标准化
# import os
# import pandas as pd
# import gffutils
# GFF_FILE = r"D:\Desktop\peptidemicro\00file\01figure\figure5\NCP_codon.gff"
# COUNT_FILE = r"D:\Desktop\peptidemicro\00file\01figure\total_all_matrix.xlsx"
# OUT_DIR = r"D:\Desktop\peptidemicro\00file\01figure"
# def prepare_length_data(gff_file):
#     db_path = gff_file + ".db"
#     if not os.path.exists(db_path):
#         print("🔄 正在创建 GFF 数据库...")
#         gffutils.create_db(
#             gff_file,
#             dbfn=db_path,
#             force=True,
#             keep_order=True,
#             merge_strategy="merge",
#             id_spec={"gene": "ID", "mRNA": "ID", "CDS": "Parent"},
#             disable_infer_genes=True,
#             disable_infer_transcripts=True
#         )
#     db = gffutils.FeatureDB(db_path)
#     gene_lengths = {}
#     for gene in db.features_of_type("gene"):
#         total_length = 0
#         for mrna in db.children(gene, featuretype="mRNA", order_by="start"):
#             exons = list(db.children(mrna, featuretype="exon", order_by="start"))
#             if exons:
#                 mrna_len = sum(e.end - e.start + 1 for e in exons)
#                 total_length += mrna_len
#         if total_length > 0:
#             gene_id = gene.id.replace("evm.model.", "evm.TU.")
#             gene_lengths[gene_id] = total_length
#     length_df = pd.DataFrame(
#         gene_lengths.items(),
#         columns=["GeneID", "length"]
#     )
#     out_len = os.path.join(OUT_DIR, "gene_lengths.csv")
#     length_df.to_csv(out_len, index=False)
#     print(f"✅ 基因长度表已生成：{out_len}")
#     return length_df
# def read_counts(count_file):
#     if count_file.lower().endswith(".csv"):
#         df = pd.read_csv(count_file)
#     else:
#         df = pd.read_excel(count_file)
#     if df.columns[0] != "GeneID":
#         df = df.rename(columns={df.columns[0]: "GeneID"})
#     return df
# def normalize_tpm(count_df, length_df, output_file):
#     df = pd.merge(count_df, length_df, on='GeneID', how='inner')
#     df = df[df['length'] > 0]
#     sample_cols = [col for col in df.columns if col not in ['GeneID', 'length']]
#     tpm_data = {}
#     for sample in sample_cols:
#         rpk = (df[sample] * 10**3) / df['length']
#         per_million_scaling_factor = rpk.sum() / 10**6
#         tpm = rpk / per_million_scaling_factor
#         tpm_data[sample] = tpm
#     tpm_df = pd.concat([df[['GeneID']], pd.DataFrame(tpm_data)], axis=1)
#     tpm_df.to_excel(output_file, index=False)
# if __name__ == "__main__":
#     length_df = prepare_length_data(GFF_FILE)
#     count_df = read_counts(COUNT_FILE)
#     normalize_tpm(count_df, length_df, os.path.join(OUT_DIR, "total_all_matrix_tpm.xlsx"))






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