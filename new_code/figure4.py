# # 转录组数据基础的NCP分布veen图
# import pandas as pd
# gene_expression_df = pd.read_csv(r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\03ouput\finally_filtered_genes.csv", index_col=0)
# tissue_mapping_df = pd.read_excel(r"D:\Desktop\peptidemicro\00file\00raw\rnaseq\Total_rna_seq.xlsx", sheet_name="group")
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
# output_excel_path = r"D:\Desktop\peptidemicro\00file\01figure\figure4\rnaseq_veen.xlsx"
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
#     expression_file = r"D:\Desktop\peptidemicro\00file\01figure\finally_filtered_genes.csv"
#     output_pdf = r"D:\Desktop\peptidemicro\00file\01figure\figure4\peptide_distribution_bar_truncated.pdf"
#     output_csv = r"D:\Desktop\peptidemicro\00file\01figure\figure4\peptide_distribution_stat.csv"
#     plot_peptide_distribution_bar_broken_y(
#         expression_file,
#         output_pdf,
#         output_csv,
#         threshold=1,
#         ylims=((0, 400), (600, 800), (2800, 3000))
#     )

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
# from Bio.Seq import Seq
# from tqdm import tqdm
# from multiprocessing import Pool
# import warnings
# warnings.filterwarnings("ignore")
# # -----------------------------
# # 1) 常量与全局变量
# # -----------------------------
# START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
# STOP_CODONS  = {"TAA", "TAG", "TGA"}
# # 负链：沿用你原脚本的“在正向基因组序列上匹配的三联体集合”
# MINUS_START_CODONS = {"CAT", "CAG", "CAC", "CAA", "CGT"}
# MINUS_STOP_CODONS  = {"TTA", "CTA", "TCA"}
# CODON_PRIOR = {
#     "ATG": 0.0,      # 95%以上，基准
#     "CTG": -0.5,     # ~2%，相对常见
#     "GTG": -1.0,     # ~1%
#     "TTG": -1.5,     # <1%
#     "ACG": -3.0     # 极罕见
# }
# MINUS_CODON_PRIOR = {
#     "CAT": 0.0,    # 对应ATG
#     "CAG": -0.5,   # 对应CTG
#     "CAC": -1.0,   # 对应GTG
#     "CAA": -1.5,   # 对应TTG
#     "CGT": -3.0   # 对应ACG
# }
# _genome_dict = None
# _max_scan_nt = None
# def init_worker(genome_file, max_scan_nt):
#     global _genome_dict, _max_scan_nt
#     _max_scan_nt = max_scan_nt
#     _genome_dict = {}
#     for rec in SeqIO.parse(genome_file, "fasta"):
#         _genome_dict[rec.id] = rec.seq
# # -----------------------------
# # 2) 读取 peptide 表、按 accession 合并区间、算 coverage
# # -----------------------------
# def merge_segments(segments):
#     min_start = min(segments, key=lambda x: x[0])[0]
#     max_end = max(segments, key=lambda x: x[1])[1]
#     return min_start, max_end
# def filter_peptide_seq_cal_cov(peptide_file, database_file):
#     database_dict = {}
#     for rec in SeqIO.parse(database_file, "fasta"):
#         database_dict[rec.id] = rec.seq
#     df_NCP = pd.read_excel(peptide_file, sheet_name="unique")
#     accession_segments = {}
#     for _, row in df_NCP.iterrows():
#         accession = row['accessions']
#         start = int(row['start'])
#         end = int(row['end'])
#         chrom = row['chrom']
#         strand = row['strand']
#         accession_segments.setdefault(accession, []).append((start, end, chrom, strand))
#     stats = []
#     for accession, segments in accession_segments.items():
#         if accession not in database_dict:
#             continue
#         seq = database_dict[accession]
#         length_aa = len(seq)
#         min_start, max_end = merge_segments(segments)
#         cov_length_aa = (max_end - min_start + 1) / 3.0
#         coverage = cov_length_aa / length_aa if length_aa > 0 else None
#         chrom = segments[0][2]
#         strand = segments[0][3]
#         peptide_count = len(segments)
#         stats.append({
#             'accession': accession,
#             'chrom': chrom,
#             'strand': strand,
#             'min_start': min_start,
#             'max_end': max_end,
#             'coverage': coverage,
#             'peptide_count': peptide_count,
#             'sequence_length_aa': length_aa
#         })
#     return stats
# # -----------------------------
# # 3) Kozak 评分
# # -----------------------------
# def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank=6):
#     genome_seq = Seq(str(genome_seq).upper())
#     if strand == '+':
#         ctx = genome_seq[phy_start - 1 - flank: phy_start + flank + 2]
#     else:
#         ctx = genome_seq[phy_end - 3 - flank: phy_end + flank].reverse_complement()
#     codon_pos = flank
#     if len(ctx) < codon_pos + 4:
#         return {"core": None, "extended": None, "total": None, "context": None}
#     core = 0.0
#     if codon_pos - 3 >= 0 and ctx[codon_pos - 3] in ('A', 'G'):
#         if ctx[codon_pos - 3] == 'A':
#             core += 2.5
#         else:
#             core += 2.0
#     if codon_pos + 3 < len(ctx) and ctx[codon_pos + 3] == 'G':
#         core += 2.0
#     return {"core": core, "extended": None, "total": core, "context": str(ctx)}
# # -----------------------------
# # 4) 枚举候选 ORF（+）
# # -----------------------------
# def _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         s = min_start - 3 * i
#         if s - 1 < 0 or s + 2 > L:
#             continue
#         triplet = seq_str[s - 1:s + 2]
#         if triplet in START_CODONS:
#             yield s, triplet, CODON_PRIOR.get(triplet, -10.0)
# def _find_stop_plus(seq_str, max_end, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for j in range(max_steps):
#         e = max_end + 3 * j
#         if e + 3 > L:
#             break
#         triplet = seq_str[e:e + 3]
#         if triplet in STOP_CODONS:
#             return e
#     return None
# def enumerate_orf_candidates_plus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
#     seq_str = str(Seq(str(genome_seq)).upper())
#     candidates = []
#     for phy_start, triplet, prior_value in _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
#         phy_end = _find_stop_plus(seq_str, max_end, max_scan_nt)
#         if phy_end is None:
#             continue
#         if not (phy_start <= min_start and phy_end >= max_end):
#             continue
#         kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='+', flank=flank)
#         if kozak["total"] is None:
#             continue
#         total_score = prior_value * 0.6 + kozak["total"] * 0.4
#         candidates.append({
#             "phy_start": phy_start,
#             "phy_end": phy_end,
#             "prior": triplet,
#             "prior_value": prior_value,
#             "kozak_total": kozak["total"],
#             "kozak_seq": kozak["context"],
#             "total_score": total_score,
#             "start_to_peptide_nt": int(min_start - phy_start)
#         })
#     candidates.sort(key=lambda d: d["total_score"], reverse=True)
#     return candidates
# # -----------------------------
# # 5) 枚举候选 ORF（-）——沿用你原坐标/切片定义
# # -----------------------------
# def _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         left = max_end + 3 * (i - 1)
#         right = max_end + 3 * i
#         if right <= 0 or left >= L:
#             continue
#         triplet = seq_str[left:right]
#         if len(triplet) != 3:
#             continue
#         if triplet in MINUS_START_CODONS:
#             phy_end = max_end + 3 * i
#             prior_value = MINUS_CODON_PRIOR.get(triplet, -10.0)
#             yield phy_end, triplet, prior_value
# def _find_stop_minus(seq_str, min_start, max_scan_nt):
#     max_steps = int(max_scan_nt / 3)
#     L = len(seq_str)
#     for i in range(max_steps):
#         left = (min_start - 1) - 3 * (i + 1)
#         right = (min_start - 1) - 3 * i
#         if right <= 0 or left >= L:
#             continue
#         triplet = seq_str[left:right]
#         if len(triplet) != 3:
#             continue
#         phy_start = min_start - 3 * i
#         if triplet in MINUS_STOP_CODONS:
#             return phy_start
#     return None
# def enumerate_orf_candidates_minus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
#     seq_str = str(Seq(str(genome_seq)).upper())
#     candidates = []
#     for phy_end, triplet_raw, prior_value in _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
#         phy_start = _find_stop_minus(seq_str, min_start, max_scan_nt)
#         if phy_start is None:
#             continue
#         if not (phy_start <= min_start and phy_end >= max_end):
#             continue
#         kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='-', flank=flank)
#         if kozak["total"] is None:
#             continue
#         total_score = prior_value * 0.6 + kozak["total"] * 0.4
#         triplet_rc = str(Seq(triplet_raw).reverse_complement())
#         candidates.append({
#             "phy_start": phy_start,
#             "phy_end": phy_end,
#             "prior": triplet_rc,
#             "prior_value": prior_value,
#             "kozak_total": kozak["total"],
#             "kozak_seq": kozak["context"],
#             "total_score": total_score,
#             "start_to_peptide_nt": int(phy_end - max_end)
#         })
#     candidates.sort(key=lambda d: d["total_score"], reverse=True)
#     return candidates
# # -----------------------------
# # 6) worker：返回 best + candidates
# # -----------------------------
# def run_scan_and_output_for_item(item):
#     global _genome_dict, _max_scan_nt
#     chrom = item['chrom']
#     strand = item['strand']
#     min_start = int(item['min_start'])
#     max_end = int(item['max_end'])
#     if chrom not in _genome_dict:
#         item['note'] = 'chrom_not_found'
#         item['best'] = None
#         item['candidates'] = []
#         return item
#     gseq = _genome_dict[chrom]
#     if strand == '+':
#         candidates = enumerate_orf_candidates_plus(gseq, min_start, max_end, _max_scan_nt, flank=6)
#     else:
#         candidates = enumerate_orf_candidates_minus(gseq, min_start, max_end, _max_scan_nt, flank=6)
#     item['best'] = candidates[0] if candidates else None
#     item['candidates'] = candidates
#     return item
# def run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=20):
#     with Pool(processes=nproc, initializer=init_worker, initargs=(genome_file, max_scan_nt)) as pool:
#         stats_update = list(pool.imap_unordered(run_scan_and_output_for_item, stats, chunksize=50))
#     return stats_update
# # -----------------------------
# # 7) main：输出两张表
# # -----------------------------
# def main():
#     peptide_file  = "finally_expressed_sp_info.xlsx"
#     database_file = "Eu_peptide_database_customized_5.fa"
#     genome_file   = "Eu_genome.fasta"
#     out_best = "output_best.csv"
#     out_cand = "output_candidates.csv"
#     # 1) 汇总 accession 覆盖信息
#     stats = filter_peptide_seq_cal_cov(peptide_file, database_file)
#     # 2) 多进程扫描并枚举候选
#     stats_update = run_scan_and_output(stats, genome_file, max_scan_nt=300, nproc=100)
#     # 3) 拆成 best 表与 candidates 表
#     best_rows = []
#     cand_rows = []
#     for it in stats_update:
#         base = {k: it[k] for k in it.keys() if k not in ("best", "candidates")}
#         # best：每个 accession 一行
#         if it.get("best") is not None:
#             best_rows.append({**base, **it["best"]})
#         else:
#             best_rows.append({
#                 **base,
#                 "phy_start": None, "phy_end": None, "prior": None,
#                 "prior_value": None, "kozak_total": None, "kozak_seq": None,
#                 "total_score": None, "start_to_peptide_nt": None
#             })
#         # candidates：每个 accession 多行（rank 从 1 开始）
#         for rank, cand in enumerate(it.get("candidates", []), start=1):
#             cand_rows.append({**base, "rank": rank, **cand})
#     df_best = pd.DataFrame(best_rows)
#     df_cand = pd.DataFrame(cand_rows)
#     df_best.to_csv(out_best, index=False)
#     df_cand.to_csv(out_cand, index=False)
#     print(f"✅ 已输出 best：{out_best}（每个 accession 最优 ORF）")
#     print(f"✅ 已输出 candidates：{out_cand}（每个 accession 所有候选 ORF，含 rank）")
# if __name__ == "__main__":
#     main()

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

import pandas as pd
in_file = r"F:\work_mechanism\new_file\02figure\figure4\codon\codon_prediction\codon_prediction_v4\candidates_scored.xlsx"
out_file = r"F:\work_mechanism\new_file\02figure\figure4\kozak_count.xlsx"
codons = ["CTG", "ACG", "GTG", "ATG", "TTG"]
SEQ_COL = "kozak_seq"
CODON_COL = "codon"
L = 13
df = pd.read_excel(in_file, sheet_name="rank1")
def count_matrix(seqs, L=13):
    bases = ['A', 'C', 'G', 'T']
    cols = [f"X{i}" for i in range(1, L + 1)]
    m = pd.DataFrame(0, index=bases, columns=cols, dtype=int)
    seqs = (pd.Series(seqs).dropna().astype(str).str.upper().str.slice(0, L))
    seqs = seqs[seqs.str.len() == L]
    seqs = seqs[seqs.str.fullmatch(rf"[ACGT]{{{L}}}")]
    for s in seqs:
        for i, ch in enumerate(s, start=1):
            m.loc[ch, f"X{i}"] += 1
    return m
with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
    count_matrix(df[SEQ_COL], L=L).to_excel(writer, sheet_name="ALL")
    codon_series = df[CODON_COL].fillna("").astype(str).str.upper()
    for c in codons:
        sub = df.loc[codon_series.str.contains(c, regex=False), SEQ_COL]
        count_matrix(sub, L=L).to_excel(writer, sheet_name=c)    