import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio.Seq import Seq
import os, gffutils

# -----------------------------
# 1) 常量与全局变量
# -----------------------------
START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
STOP_CODONS  = {"TAA", "TAG", "TGA"}
MINUS_START_CODONS = {"CAT", "CAG", "CAC", "CAA", "CGT"}
MINUS_STOP_CODONS  = {"TTA", "CTA", "TCA"}
MAX_SCAN_NT = 600
THREADS = 3

def prepare_gffutils_db(trans_file, db_path):
    if os.path.exists(db_path): return
    gffutils.create_db(trans_file, dbfn=db_path, force=True,
        keep_order=True, merge_strategy="merge", sort_attribute_values=True,
        disable_infer_transcripts=True, disable_infer_genes=True)
def run_scan_and_output(records, hit_transcript_gtf, genome_fasta, db_path, max_scan_nt=MAX_SCAN_NT, nproc=THREADS):
    prepare_gffutils_db(hit_transcript_gtf, db_path)
    with Pool(processes=nproc, initializer=init_worker, initargs=(genome_fasta, db_path, max_scan_nt)) as pool:
        stats = list(pool.imap_unordered(run_scan_and_output_for_item, records, chunksize=50))
    return stats
def main():
    peptide_info = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\finally_expressed_sp_info_annotated.xlsx"
    hit_transcript_gtf = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\hit_transcripts.gtf"
    genome_fasta = r"F:\work_mechanism\new_file\02figure\Eu_genome_modified\Eu_genome.fasta"
    out_cand = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\output.csv"
    db_path = r"F:\work_mechanism\new_file\02figure\figure4\new_transcript\hit_transcripts.db"
    df = pd.read_excel(peptide_info, sheet_name="annotated")
    records = df.to_dict("records")
    stats = run_scan_and_output(records, hit_transcript_gtf, genome_fasta, db_path, max_scan_nt=MAX_SCAN_NT, nproc=THREADS)

if __name__ == "__main__":
    main()