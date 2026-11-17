from pathlib import Path
import os
import re
import csv
import time
import bisect
from collections import defaultdict

import pandas as pd
from Bio.Seq import Seq


BASE = str(Path(__file__).resolve().parents[1])
OUT = f"{BASE}/output/02filter_trans"          
RESULT_XLSX = f"{BASE}/output/01location/NCP_peptides.xlsx"   
NOVEL_GFF = f"{OUT}/novel.gff3"
GENOME_FA = f"{BASE}/00raw_data/Eu_genome.fasta"
CLASS_TSV = f"{OUT}/classification.tsv"
TPM_MATRIX = f"{OUT}/TPM_matrix.tsv"

MIN_LEN = 5                       
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
EXPRESS_MIN_TPM = 1.0             
EXPRESS_MIN_SAMPLES = 3

_T0 = time.time()

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def rc(s):
    return str(Seq(s).reverse_complement())

# ================================================================== Stage 0
def load_peptide_intervals():
    """Return unique peptide intervals with recalculated peptide and nucleotide lengths."""
    log("Stage0a Reading peptide table...")
    df = pd.read_excel(
        RESULT_XLSX, sheet_name="NCP",
        usecols=["sequence", "chrom", "start", "end", "strand", "frame",
                 "psms", "accessions"])
    n_raw = len(df)

    df["sequence"] = df["sequence"].astype(str).str.upper().str.strip()
    df = df[df["strand"].isin(["+", "-"])]
    df = df.dropna(subset=["chrom", "start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    n_coord = len(df)

    df["pep_len"] = df["sequence"].str.len()
    n_too_short = int((df["pep_len"] < MIN_LEN).sum())
    df = df[df["pep_len"] >= MIN_LEN]
    df = df[df["sequence"].apply(lambda s: set(s) <= STANDARD_AA)]
    n_usable = len(df)


    df["nt_len"] = df["end"] - df["start"] + 1
    bad = df[df["nt_len"] != 3 * df["pep_len"]]
    if len(bad):
        bad.to_csv(f"{OUT}/anomaly_interval_length.tsv", sep="\t", index=False)
    df = df[df["nt_len"] == 3 * df["pep_len"]]


    df["frame"] = df["frame"].astype(str)
    dedup = (df.groupby(["sequence", "chrom", "start", "end", "strand", "frame"],
                        as_index=False)
               .agg(psms=("psms", "max"), accessions=("accessions", "first")))

    dedup["pep_len"] = dedup["sequence"].str.len()
    dedup["nt_len"] = dedup["end"] - dedup["start"] + 1
    log(f"  Raw rows: {n_raw}; with coordinates: {n_coord}; usable: {n_usable} "
        f"(too short: {n_too_short}); invalid interval lengths: {len(bad)}; "
        f"unique intervals: {len(dedup)}")
    return dedup

def load_novel_gff():
    """Load novel gff."""
    # The GFF geneID attribute is authoritative for locus assignment; deriving
    # locus IDs by trimming transcript IDs can misidentify assembled isoforms.
    log("Stage0b Reading novel.gff3...")
    tx, exons, gene = {}, defaultdict(list), {}
    with open(NOVEL_GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] in ("transcript", "mRNA"):
                tid = re.search(r"ID=([^;]+)", f[8]).group(1)
                tx[tid] = (f[0], f[6])



                m_g = re.search(r'geneID="?([^";]+)"?', f[8])
                gene[tid] = m_g.group(1) if m_g else tid.rsplit(".", 1)[0]
            elif f[2] == "exon":
                tid = re.search(r"Parent=([^;]+)", f[8]).group(1)
                exons[tid].append((int(f[3]), int(f[4])))

    index = defaultdict(lambda: ([], [], [], []))
    by_key = defaultdict(list)
    for tid, lst in exons.items():
        if tid not in tx:
            continue                       
        chrom, strand = tx[tid]
        for s, e in lst:
            by_key[(chrom, strand)].append((s, e, tid))
    for key, lst in by_key.items():
        lst.sort()
        starts = [x[0] for x in lst]
        ends = [x[1] for x in lst]
        tids = [x[2] for x in lst]
        premax = [0] * len(lst)
        m = 0
        for i, e in enumerate(ends):       
            if e > m:
                m = e
            premax[i] = m
        index[key] = (starts, ends, tids, premax)

    exons_g = {t: sorted(v) for t, v in exons.items()}
    log(f"  Transcripts: {len(tx)}; exons: {sum(len(v) for v in exons_g.values())}, "
        f"index keys (chrom,strand): {len(index)}")
    return tx, gene, exons_g, dict(index)

def load_genome():
    log("Stage0c Loading genome...")
    g, name, buf = {}, None, []
    with open(GENOME_FA) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    g[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        g[name] = "".join(buf)
    g = {k: v.upper() for k, v in g.items()}
    log(f"  {len(g)} chromosomes; total {sum(map(len, g.values()))/1e9:.2f} Gb")
    return g

def load_class_codes():
    cls = {}
    with open(CLASS_TSV) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r)
        for row in r:
            cls[row[0]] = row[2]
    log(f"Stage0d class codes: {len(cls)} entries")
    return cls

# ================================================================== Stage 1
def covered_by_exons(exon_list, ps, pe):
    """Covered by exons."""
    cur = ps
    for s, e in exon_list:
        if e < cur:            
            continue
        if s > cur:            
            return False
        cur = e + 1
        if cur > pe:
            return True
    return cur > pe

def find_candidate_transcripts(intervals, index, exons_g):
    """Find candidate transcripts."""
    log("Stage1 Genomic matching (same strand and full exon coverage)...")
    result, reason = {}, {}
    n_with_cand = 0
    n_no_ss, n_no_ovl, n_not_contained = 0, 0, 0
    for i, row in intervals.iterrows():
        entry = index.get((row["chrom"], row["strand"]))
        if entry is None:           
            result[i] = []
            reason[i] = "no_same_strand_transcript"
            n_no_ss += 1
            continue
        starts, ends, tids, premax = entry
        ps, pe = int(row["start"]), int(row["end"])


        j = bisect.bisect_right(starts, pe) - 1
        ovl = set()
        while j >= 0 and premax[j] >= ps:
            if ends[j] >= ps:
                ovl.add(tids[j])
            j -= 1
        if not ovl:
            result[i] = []
            reason[i] = "no_exon_overlap"
            n_no_ovl += 1
            continue

        cand = [t for t in ovl if covered_by_exons(exons_g[t], ps, pe)]
        result[i] = cand
        if cand:
            n_with_cand += 1
        else:
            reason[i] = "exon_not_fully_contained"
            n_not_contained += 1
    log(f"  {n_with_cand}/{len(intervals)} intervals have transcript candidates; unmatched: "
        f"no same-strand transcript: {n_no_ss}; no exon overlap: {n_no_ovl}, "
        f"insufficient exon coverage: {n_not_contained}")
    return result, reason

# ================================================================== Stage 2
def extract_cdna_segment(genome, chrom, exons_tx_order, strand, ps, pe):
    """Extract cdna segment."""
    pieces = []
    for s, e in exons_tx_order:
        lo, hi = max(s, ps), min(e, pe)
        if lo > hi:
            continue
        seg = genome[chrom][lo - 1:hi]
        pieces.append(seg if strand == "+" else rc(seg))
    return "".join(pieces)

def verify_frames(intervals, candidates, tx, exons_g, genome, reason):
    """Verify frames."""
    log("Stage2 Reading-frame validation (translated cDNA comparison)...")
    pairs = []
    n_fail = 0
    n_int_fail = 0
    passed_interval = set()
    cache_order = {}
    for i, row in intervals.iterrows():
        cands = candidates.get(i, [])
        if not cands:
            continue
        pseq = row["sequence"]
        ok_any = False
        for tid in cands:
            chrom, strand = tx[tid]
            if tid not in cache_order:
                lst = exons_g[tid]
                cache_order[tid] = lst if strand == "+" else lst[::-1]
            nt = extract_cdna_segment(genome, chrom, cache_order[tid],
                                      strand, int(row["start"]), int(row["end"]))
            aa = str(Seq(nt).translate()) if len(nt) == len(pseq) * 3 else ""
            if aa == pseq:
                pairs.append((i, tid))
                ok_any = True
            else:
                n_fail += 1
        if ok_any:
            passed_interval.add(i)
        else:
            reason[i] = "translation_mismatch"
            n_int_fail += 1
    log(f"  Passed: {len(pairs)} pairs; removed due to translation mismatch: {n_fail} pairs "
        f"(involving {n_int_fail} intervals)")
    return pairs

# ================================================================== Stage 3
def build_pair_table(intervals, pairs, locus_of_tid, cls):
    log("Stage3 Building match pairs and ambiguity labels...")
    rows = []
    for i, tid in pairs:
        row = intervals.loc[i]
        rows.append({
            "peptide": row["sequence"],
            "peptide_len": int(row["pep_len"]),
            "chrom": row["chrom"],
            "start": int(row["start"]),
            "end": int(row["end"]),
            "nt_len": int(row["nt_len"]),
            "strand": row["strand"],
            "frame_xlsx": row["frame"],
            "psms": row["psms"],
            "transcript_id": tid,
            "locus": locus_of_tid[tid],
            "class_code": cls.get(tid, "NA"),
        })
    pdf = pd.DataFrame(rows)



    amb = (pdf.groupby(["peptide", "chrom", "start", "end", "strand"])["locus"]
              .agg(n_tid="count", n_locus="nunique"))
    def _amb(r):
        if r["n_locus"] > 1:
            return "multi_locus"
        return "same_locus" if r["n_tid"] > 1 else "unique"
    amb["ambiguity_type"] = amb.apply(_amb, axis=1)
    pdf = pdf.merge(amb["ambiguity_type"],
                    left_on=["peptide", "chrom", "start", "end", "strand"],
                    right_index=True, how="left")
    return pdf, locus_of_tid

def write_peptide_summary(pair_df):



    npos = (pair_df.drop_duplicates(["peptide", "chrom", "start", "end", "strand"])
                  .groupby("peptide").size().rename("n_genomic_positions"))
    nchrom = (pair_df.drop_duplicates(["peptide", "chrom"])
                    .groupby("peptide")["chrom"].nunique().rename("n_chroms"))
    agg = pair_df.groupby("peptide").agg(
        peptide_len=("peptide_len", "first"),
        psms=("psms", "max"),
        n_transcripts=("transcript_id", "nunique"),
        n_loci=("locus", "nunique"),
        class_codes=("class_code", lambda s: ",".join(sorted(set(s)))),
        ambiguity=("ambiguity_type",
                   lambda s: ("multi_locus" if (s == "multi_locus").any()
                              else "same_locus" if (s == "same_locus").any()
                              else "unique")),
        example_transcripts=("transcript_id",
                             lambda s: ",".join(sorted(s)[:5])),
    ).join(npos).join(nchrom)
    agg = agg.reset_index()
    agg.to_csv(f"{OUT}/peptide_summary.tsv", sep="\t", index=False)
    return agg

def write_unmatched(intervals, pairs):
    """Write unmatched."""
    matched_iv = {i for i, _ in pairs}
    matched_pep = {intervals.loc[i]["sequence"] for i in matched_iv}
    rows = []
    for i, row in intervals.iterrows():
        if i in matched_iv:
            continue
        rows.append((row["sequence"], int(row["pep_len"]), row["psms"],
                     row["chrom"], int(row["start"]), int(row["end"]),
                     row["strand"], row["frame"], row["reason"]))
    with open(f"{OUT}/unmatched_peptides.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["peptide", "peptide_len", "psms", "chrom", "start", "end",
                    "strand", "frame", "reason"])
        w.writerows(rows)
    dist = defaultdict(int)
    for r in rows:
        dist[r[8]] += 1
    n_unmatched_pep = len({r[0] for r in rows} - matched_pep)
    return len(rows), dict(dist), n_unmatched_pep

# ================================================================== Stage 4
def link_expression(pair_df, locus_of_tid):
    log("Stage4 Joining expression matrix...")
    verified_tids = set(pair_df["transcript_id"])
    verified_loci = {locus_of_tid[t] for t in verified_tids}
    want = {}
    gene_checked = gene_mismatch = 0
    with open(TPM_MATRIX) as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        i_tid, i_gene = header.index("transcript_id"), header.index("gene_id")
        samples = header[3:]
        for row in r:
            tid = row[i_tid]
            gene = row[i_gene]
            if tid in verified_tids:
                gene_checked += 1
                if locus_of_tid.get(tid) != gene:
                    gene_mismatch += 1
            if gene in verified_loci:
                want[tid] = (gene, [float(x) for x in row[3:]])
    if gene_mismatch:
        log(f"  [FATAL] TPM matrix gene_id disagrees with GFF3 locus: {gene_mismatch} entries, "
            f"aborting (verify TPM_matrix.tsv and novel.gff3 originate from the same assembly)")
        raise SystemExit(2)
    if gene_checked != len(verified_tids):
        log(f"  [FATAL] {len(verified_tids) - gene_checked} verified transcripts are missing "
            f"from TPM matrix; aborting")
        raise SystemExit(2)
    missing = verified_loci - {g for g, _ in want.values()}
    if missing:
        log(f"  [FATAL] {len(missing)} verified loci are missing from TPM matrix: "
            f"{sorted(missing)[:5]}; aborting")
        raise SystemExit(2)
    log(f"  gene_id consistency check passed: {gene_checked} verified transcripts agree with GFF3 loci"
        f"; retrieved {len(want)} transcripts covering {len(verified_loci)} verified loci")


    meta = pair_df.groupby("transcript_id").agg(
        n_peptides=("peptide", "nunique"),
        class_codes=("class_code", lambda s: ",".join(sorted(set(s)))),
    )
    with open(f"{OUT}/verified_transcript_TPM.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["transcript_id", "locus", "class_code", "n_peptides",
                    "mean_TPM", "max_TPM", "expressed(TPM>=1 in >=3 samples)"]
                   + samples)
        for tid in sorted(verified_tids):
            gene, vals = want[tid]
            ok = sum(1 for v in vals if v >= EXPRESS_MIN_TPM) >= EXPRESS_MIN_SAMPLES
            w.writerow([tid, gene,
                        meta.loc[tid, "class_codes"],
                        meta.loc[tid, "n_peptides"],
                        f"{sum(vals)/len(vals):.3f}", f"{max(vals):.3f}",
                        "Y" if ok else "N"] + [f"{v:.3f}" for v in vals])


    loc_agg = defaultdict(lambda: [0.0] * len(samples))
    loc_n_iso = defaultdict(int)
    for tid, (gene, vals) in want.items():
        loc_n_iso[gene] += 1
        acc = loc_agg[gene]
        for k, v in enumerate(vals):
            acc[k] += v
    pep_by_locus = pair_df.groupby("locus")["peptide"].nunique().to_dict()
    codes_by_locus = (pair_df.groupby("locus")["class_code"]
                        .agg(lambda s: ",".join(sorted(set(s)))).to_dict())
    with open(f"{OUT}/locus_expression.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["locus", "n_isoforms_total", "n_peptides_supported",
                    "class_codes", "mean_TPM(sum_isoforms)", "max_TPM(sum_isoforms)",
                    "expressed(TPM>=1 in >=3 samples)"] + samples)
        for gene in sorted(loc_agg):
            vals = loc_agg[gene]
            ok = sum(1 for v in vals if v >= EXPRESS_MIN_TPM) >= EXPRESS_MIN_SAMPLES
            w.writerow([gene, loc_n_iso[gene], pep_by_locus.get(gene, 0),
                        codes_by_locus.get(gene, "NA"),
                        f"{sum(vals)/len(vals):.3f}", f"{max(vals):.3f}",
                        "Y" if ok else "N"] + [f"{v:.3f}" for v in vals])
    return samples

# ================================================================== main
def main():
    os.makedirs(OUT, exist_ok=True)
    intervals = load_peptide_intervals()
    tx, locus_of_tid, exons_g, index = load_novel_gff()
    genome = load_genome()
    cls = load_class_codes()

    candidates, reason = find_candidate_transcripts(intervals, index, exons_g)
    pairs = verify_frames(intervals, candidates, tx, exons_g, genome, reason)

    n_unmatched_iv, unmatched_dist, n_unmatched_pep = write_unmatched(
        intervals.assign(reason=pd.Series(reason)), pairs)

    if pairs:
        pair_df, locus_of_tid = build_pair_table(intervals, pairs,
                                                  locus_of_tid, cls)
        pair_df.to_csv(f"{OUT}/peptide_transcript_pairs.tsv", sep="\t",
                       index=False)
        pep_agg = write_peptide_summary(pair_df)
        link_expression(pair_df, locus_of_tid)

        n_pep = len(pep_agg)
        n_loci = pair_df["locus"].nunique()
        n_tx = pair_df["transcript_id"].nunique()
        loc_any_u = (pair_df.assign(is_u=pair_df["class_code"] == "u")
                           .groupby("locus")["is_u"].any().sum())
        amb_counts = pep_agg["ambiguity"].value_counts().to_dict()
        n_multi_pos = int((pep_agg["n_genomic_positions"] > 1).sum())
        # Uniqueness is conditional on the supplied genomic coordinate; the
        # upstream peptide database records only the first sequence match.
        print(f"""
# Coordinate-based peptide/transcript matching summary
- Unique peptide intervals: {len(intervals)}
- Pairs passing same-strand exon coverage and cDNA translation checks: {len(pair_df)}
- Verified peptides (unique sequences): {n_pep}; unmatched unique peptides: {n_unmatched_pep}
- Unmatched intervals: {n_unmatched_iv}; reason counts: {unmatched_dist}
- Verified transcripts: {n_tx}; verified loci: {n_loci} (including {loc_any_u} novel class-u loci)
- Peptide ambiguity counts: {amb_counts}
- Peptides with multiple recorded genomic positions: {n_multi_pos}
- Interpretation: Excel positions come from the first match returned by db_seq.find().
  Therefore, unique means unique among novel transcripts at the supplied genomic
  coordinate, not a unique protein origin across the genome. Translation validation
  supports the reported spliced transcript segment but does not establish the true ORF/TIS.
- Output directory: {OUT}
""")
    else:
        log("No validated peptide-transcript matches")
    log(f"Total runtime: {time.time()-_T0:.0f}s")

if __name__ == "__main__":
    main()
