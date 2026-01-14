import os
import pandas as pd
import gffutils
def build_or_load_db(anno_file, db_path="transcriptome.db"):
    if os.path.exists(db_path):
        return gffutils.FeatureDB(db_path, keep_order=True)
    db = gffutils.create_db(
        anno_file,
        dbfn=db_path,
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True
    )
    return db
def gtf_attr_str(attrs: dict):
    parts = []
    for k, vs in attrs.items():
        if vs is None: continue
        if not isinstance(vs, (list, tuple)): vs = [vs]
        for v in vs: parts.append(f'{k} "{v}";')
    return " ".join(parts)
def feature_to_gtf_line(f):
    seqid = f.seqid
    source = f.source if f.source is not None else "."
    featuretype = f.featuretype
    start = int(f.start)
    end = int(f.end)
    score = f.score if f.score is not None else "."
    strand = f.strand if f.strand is not None else "."
    frame = f.frame if f.frame is not None else "."
    attr = gtf_attr_str(f.attributes)
    return "\t".join([seqid, source, featuretype, str(start), str(end), str(score), strand, str(frame), attr])
def overlap(a_start, a_end, b_start, b_end):
    return a_start <= b_start <= b_end <= a_end
def annotate_peptides_and_export(peptide_df, db, out_peptide_xlsx, out_hit_gtf):
    hit_tx_ids = set()
    add_cols = {
        "hit_transcript_ids": [],
        "hit_tx_start": [],
        "hit_tx_end": [],
        "hit_tx_strand": [],
    }
    for _, row in peptide_df.iterrows():
        chrom = str(row["chrom"])
        pep_start = int(row["start"])
        pep_end = int(row["end"])
        transcripts = db.features_of_type("transcript", seqid=chrom)
        cur_ids, cur_starts, cur_ends, cur_strands = [], [], [], []
        for tx in transcripts:
            if not overlap(tx.start, tx.end, pep_start, pep_end): continue
            exon_hit = False
            for exon in db.children(tx, featuretype="exon", order_by="start"):
                if overlap(exon.start, exon.end, pep_start, pep_end):
                    exon_hit = True
                    break
            if not exon_hit: continue
            # 命中这个 transcript
            hit_tx_ids.add(tx.id)
            cur_ids.append(tx.id)
            cur_starts.append(str(int(tx.start)))
            cur_ends.append(str(int(tx.end)))
            cur_strands.append(tx.strand if tx.strand is not None else ".")
        add_cols["hit_transcript_ids"].append(";".join(cur_ids) if cur_ids else "")
        add_cols["hit_tx_start"].append(";".join(cur_starts) if cur_starts else "")
        add_cols["hit_tx_end"].append(";".join(cur_ends) if cur_ends else "")
        add_cols["hit_tx_strand"].append(";".join(cur_strands) if cur_strands else "")
    # 1) 输出：肽段表追加命中信息
    out_df = peptide_df.copy()
    for k, v in add_cols.items(): out_df[k] = v
    out_df.to_excel(out_peptide_xlsx, index=False)
    # 2) 输出：命中的转录本 + 全部 exon 到新的 GTF
    with open(out_hit_gtf, "w") as fw:
        fw.write("# subset GTF: transcripts validated by peptides (exon-overlap)\n")
        for tx_id in sorted(hit_tx_ids):
            tx = db[tx_id]
            fw.write(feature_to_gtf_line(tx) + "\n")
            for exon in db.children(tx, featuretype="exon", order_by="start"):
                fw.write(feature_to_gtf_line(exon) + "\n")
def main():
    excel_file = "peptide_data.xlsx"
    sheet_name = "unique"
    anno_file = "merged.gtf"
    out_peptide_xlsx = "peptide_annotated.xlsx"
    out_hit_gtf = "hit_transcripts.gtf"
    db_path = "transcriptome.db"
    peptide_df = pd.read_excel(excel_file, sheet_name=sheet_name)
    for col in ["chrom", "start", "end"]:
        if col not in peptide_df.columns: raise ValueError(f"Excel 缺少必要列: {col}")
    db = build_or_load_db(anno_file, db_path=db_path)
    annotate_peptides_and_export(peptide_df, db, out_peptide_xlsx, out_hit_gtf)
    print(f"完成：\n- {out_peptide_xlsx}\n- {out_hit_gtf}\n- (db) {db_path}")
if __name__ == "__main__":
    main()