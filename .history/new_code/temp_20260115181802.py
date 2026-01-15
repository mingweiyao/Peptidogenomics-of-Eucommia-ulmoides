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
def init_worker(genome_fasta, db_path, max_scan_nt):
    global _genome_dict, _max_scan_nt, _trans_db
    _max_scan_nt = max_scan_nt
    _genome_dict = {}
    for rec in SeqIO.parse(genome_fasta, "fasta"): _genome_dict[rec.id] = rec.seq
    _trans_db = gffutils.FeatureDB(db_path, keep_order=True)

def _get_exons_for_transcript(db, transcript_id):
    tx = db[transcript_id]
    strand = tx.strand
    chrom = tx.chrom
    exons = list(db.children(tx, featuretype="exon", order_by="start"))
    exons_coords = [(e.start, e.end) for e in exons]
    if strand is None or chrom is None or not exons_coords:
        raise KeyError(f"Cannot find exons for transcript_id={transcript_id}")
    return strand, chrom, exons_coords
def _build_transcript_seq_and_mapper(genome_seq, exons_coords, strand):
    cum = 0
    exon_blocks = []
    parts = []
    for (s, e) in exons_coords:
        exon_blocks.append((s, e, cum))
        parts.append(genome_seq[s - 1 : e])
        cum += (e - s + 1)
    spliced = Seq("").join(parts)
    if strand == "-": spliced = spliced.reverse_complement()
    def map_genomic_pos_to_tpos(gpos):
        for (s, e, cum_before) in exon_blocks:
            if s <= gpos <= e:
                offset = gpos - s
                return cum_before + offset + 1
        return None
    return spliced, map_genomic_pos_to_tpos
def get_hit_transcript_dna_and_coords(hit_id, gstart, gend):
    global _genome_dict, _trans_db
    strand, chrom, exons_coords = _get_exons_for_transcript(_trans_db, hit_id)
    if chrom not in _genome_dict: raise KeyError(f"Chrom {chrom} not found in genome fasta ids")
    genome_seq = _genome_dict[chrom]
    tx_seq, mapper = _build_transcript_seq_and_mapper(genome_seq, exons_coords, strand)
    t1 = mapper(gstart)
    t2 = mapper(gend)
    tstart, tend = (t1, t2) if t1 <= t2 else (t2, t1)
    return str(tx_seq), (tstart, tend), strand, chrom
def run_scan_and_output_for_item(item):
    global _max_scan_nt
    gstart = int(item["start"])
    gend = int(item["end"])
    gstrand = str(item['strand'])
    hit_trans_id = item["hit_transcript_ids"]
    hit_trans_ids = hit_trans_id.split(";")
    per_hit = []
    for hit_id in hit_trans_ids:
        tx_seq, (tstart, tend), tx_strand, chrom = get_hit_transcript_dna_and_coords(hit_id, gstart, gend)           

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