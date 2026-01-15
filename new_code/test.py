import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from Bio.Seq import Seq
import os, gffutils

MAX_SCAN_NT = 600
THREADS = 3

# -----------------------------
# 1) 常量与全局变量
# -----------------------------
START_CODONS = {"ATG", "CTG", "GTG", "TTG", "ACG"}
STOP_CODONS  = {"TAA", "TAG", "TGA"}
MINUS_START_CODONS = {"CAT", "CAG", "CAC", "CAA", "CGT"}
MINUS_STOP_CODONS  = {"TTA", "CTA", "TCA"}

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

def kozak_score_cal(genome_seq, phy_start, phy_end, strand, flank=6):
    genome_seq = Seq(str(genome_seq).upper())
    if strand == '+': ctx = genome_seq[phy_start - 1 - flank: phy_start + flank + 2]
    else: ctx = genome_seq[phy_end - 3 - flank: phy_end + flank].reverse_complement()
    return {"context": str(ctx)}
def _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        s = min_start - 3 * i
        if s - 1 < 0 or s + 2 > L: continue
        triplet = seq_str[s - 1:s + 2]
        if triplet in START_CODONS:
            yield s, triplet
def _find_stop_plus(seq_str, max_end, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for j in range(max_steps):
        e = max_end + 3 * j
        if e + 3 > L: break
        triplet = seq_str[e:e + 3]
        if triplet in STOP_CODONS: return e
    return None
def enumerate_orf_candidates_plus(genome_seq, min_start, max_end, max_scan_nt=300, flank=6):
    seq_str = str(Seq(str(genome_seq)).upper())
    candidates = []
    for phy_start, triplet in _iter_start_candidates_plus(seq_str, min_start, max_scan_nt):
        phy_end = _find_stop_plus(seq_str, max_end, max_scan_nt)
        if phy_end is None: continue
        if not (phy_start <= min_start and phy_end >= max_end): continue
        kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='+', flank=flank)
        candidates.append({
            "phy_start": phy_start,
            "phy_end": phy_end,
            "prior": triplet,
            "kozak_seq": kozak["context"],
            "start_to_peptide_nt": int(min_start - phy_start)
        })
    return candidates
def _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        left = max_end + 3 * (i - 1)
        right = max_end + 3 * i
        if right <= 0 or left >= L: continue
        triplet = seq_str[left:right]
        if len(triplet) != 3: continue
        if triplet in MINUS_START_CODONS:
            phy_end = max_end + 3 * i
            yield phy_end, triplet
def _find_stop_minus(seq_str, min_start, max_scan_nt):
    max_steps = int(max_scan_nt / 3)
    L = len(seq_str)
    for i in range(max_steps):
        left = (min_start - 1) - 3 * (i + 1)
        right = (min_start - 1) - 3 * i
        if right <= 0 or left >= L: continue
        triplet = seq_str[left:right]
        if len(triplet) != 3: continue
        phy_start = min_start - 3 * i
        if triplet in MINUS_STOP_CODONS: return phy_start
    return None
def enumerate_orf_candidates_minus(genome_seq, min_start, max_end, max_scan_nt=MAX_SCAN_NT, flank=6):
    seq_str = str(Seq(str(genome_seq)).upper())
    candidates = []
    for phy_end, triplet_raw in _iter_start_candidates_minus(seq_str, max_end, max_scan_nt):
        phy_start = _find_stop_minus(seq_str, min_start, max_scan_nt)
        if phy_start is None: continue
        if not (phy_start <= min_start and phy_end >= max_end): continue
        kozak = kozak_score_cal(seq_str, phy_start, phy_end, strand='-', flank=flank)
        triplet_rc = str(Seq(triplet_raw).reverse_complement())
        candidates.append({
            "phy_start": phy_start,
            "phy_end": phy_end,
            "prior": triplet_rc,
            "kozak_seq": kozak["context"],
            "start_to_peptide_nt": int(phy_end - max_end)
        })
    return candidates

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
        if gstrand == '+':
            candidates = enumerate_orf_candidates_plus(tx_seq, tstart, tend, _max_scan_nt, flank=6)
        else:
            candidates = enumerate_orf_candidates_minus(tx_seq, tstart, tend, _max_scan_nt, flank=6)   
            per_hit.append({
                "hit_id": hit_id,
                "chrom": chrom,
                "tx_strand": tx_strand,
                "tstart": tstart,
                "tend": tend,
                "candidates": candidates,
                "error": None
            })
    item["per_hit"] = per_hit
    return item     

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
    cand_rows = []
    for it in stats:
        base_item = {k: it[k] for k in it.keys() if k != "per_hit"}
        for hit in it.get("per_hit", []):
            base_hit = {
                **base_item,
                "hit_id": hit.get("hit_id"),
                "chrom": hit.get("chrom"),
                "tx_strand": hit.get("tx_strand"),
                "tstart": hit.get("tstart"),
                "tend": hit.get("tend"),
                "hit_error": hit.get("error"),
            }
            for rank, cand in enumerate(hit.get("candidates", []), start=1):
                cand_rows.append({**base_hit, "rank": rank, **cand})
    df_cand = pd.DataFrame(cand_rows)
    df_cand.to_csv(out_cand, index=False)  
if __name__ == "__main__":
    main()