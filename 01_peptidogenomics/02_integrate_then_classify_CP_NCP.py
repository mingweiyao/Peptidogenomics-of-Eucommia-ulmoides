import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

NON_CDS = {"intergenic", "intron", "junction", "three_prime_UTR", "five_prime_UTR"}


def parse_args():
    base = Path(__file__).resolve().parent.parent
    p = argparse.ArgumentParser(description="Integrate exon-spanning coordinates and export CP/NCP separately")
    p.add_argument("--input", default=str(base / "output" / "01location" / "peptide_analysis_results.xlsx"),
                   help="Localization Excel containing Annotated_Proteins and Annotated_Peptides")
    p.add_argument("--protein-fa", default=str(base / "00raw_data" / "Eu_peptide_database.fa"),
                   help="Reference protein FASTA with genomic CDS blocks in the Position header field")
    p.add_argument("--output-dir", default=None, help="Output directory (defaults to input Excel directory)")
    return p.parse_args()


def load_cds_blocks(fasta_path):
    """Parse reference-protein CDS blocks in ascending genomic order."""
    cds = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):

        match = re.search(r"Position=[^:]*:\s*([0-9,\s-]+?)\s*:\s*[+-](?=\s|$)", record.description)
        if match is None:
            raise ValueError(f"Reference protein FASTA lacks a valid Position field: {record.id}, header={record.description[:200]}")
        blocks = []
        for item in match.group(1).split(","):
            item = item.strip()
            if not item:
                continue
            m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", item)
            if m is None:
                raise ValueError(f"Cannot parse CDS block: {record.id}, {item!r}")
            s, e = int(m.group(1)), int(m.group(2))
            if s < 1 or s > e:
                raise ValueError(f"Invalid CDS block coordinates: {record.id}, {s}-{e}")
            blocks.append((s, e))
        if not blocks:
            raise ValueError(f"No CDS blocks found: {record.id}")
        blocks.sort()
        if any(blocks[j][0] <= blocks[j - 1][1] for j in range(1, len(blocks))):
            raise ValueError(f"Overlapping CDS blocks: {record.id}, {blocks}")
        cds[record.id] = blocks
    if not cds:
        raise ValueError("No CDS blocks parsed from reference protein FASTA")
    return cds


def exact_int(value, label):
    if pd.isna(value):
        raise ValueError(f"{label} is missing")
    try:
        num = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{label} is not numeric: {value!r}") from None
    if not num.is_integer():
        raise ValueError(f"{label} is not an integer: {value!r}")
    return int(num)


def merge_coding_coordinates(pro, cds):
    """Reconstruct exon-spanning CDS coordinates without modifying the input data on errors."""
    need = {"accessions", "sequence", "start", "end", "first_boundary", "second_boundary"}
    missing = sorted(need - set(pro.columns))
    if missing:
        raise ValueError(f"Missing columns in Annotated_Proteins: {missing}")
    starts, ends, errors = [], [], []
    changed, multi, no_change = 0, 0, 0
    for rownum, row in enumerate(pro.itertuples(index=False, name=None), start=2):

        if rownum == 2:
            idx = {name: k for k, name in enumerate(pro.columns)}
        s_val, e_val = row[idx["start"]], row[idx["end"]]
        b1, b2 = row[idx["first_boundary"]], row[idx["second_boundary"]]
        accession = row[idx["accessions"]]
        seq = row[idx["sequence"]]
        has_one_boundary = pd.notna(b1) != pd.notna(b2)
        if has_one_boundary:
            errors.append((rownum, accession, "only_one_boundary", s_val, e_val))
            starts.append(s_val)
            ends.append(e_val)
            continue

        if pd.isna(s_val) and pd.isna(e_val) and pd.isna(b1) and pd.isna(b2):
            starts.append(s_val)
            ends.append(e_val)
            no_change += 1
            continue
        try:
            s0 = exact_int(s_val, "start")
            e0 = exact_int(e_val, "end")
            if not 1 <= s0 <= e0:
                raise ValueError(f"Reversed or invalid coordinates {s0}-{e0}")
            if not isinstance(seq, str) or not seq:
                raise ValueError(f"Invalid sequence: {seq!r}")
            length = len(seq.strip())

            if e0 - s0 + 1 == length * 3 and pd.isna(b1) and pd.isna(b2):
                starts.append(s_val)
                ends.append(e_val)
                no_change += 1
                continue
            if not isinstance(accession, str) or accession not in cds:
                raise ValueError(f"Accession not found in FASTA: {accession!r}")
            blocks = []
            for s, e in cds[accession]:
                lo, hi = max(s, s0), min(e, e0)
                if lo <= hi:
                    blocks.append((lo, hi))
            if not blocks:
                raise ValueError("Coordinates do not overlap the corresponding protein CDS")
            block_nt = sum(e - s + 1 for s, e in blocks)
            if block_nt != 3 * length:
                raise ValueError(f"CDS block length {block_nt} != peptide length * 3 = {3 * length}; blocks={blocks}")
            if pd.notna(b1):
                if len(blocks) < 2:
                    raise ValueError(f"Original boundary indicates exon spanning, but only {len(blocks)} block(s) found")
                ob1 = exact_int(b1, "first_boundary")
                ob2 = exact_int(b2, "second_boundary")
                if (blocks[0][1], blocks[1][0]) != (ob1, ob2):
                    raise ValueError(
                        f"Original boundary ({ob1},{ob2}) differs from first CDS junction "
                        f"({blocks[0][1]},{blocks[1][0]}); blocks={blocks}"
                    )
            if len(blocks) == 1:
                starts.append(blocks[0][0])
                ends.append(blocks[0][1])
            else:
                starts.append("|".join(str(s) for s, _ in blocks))
                ends.append("|".join(str(e) for _, e in blocks))
                changed += 1
                if len(blocks) >= 3:
                    multi += 1
        except ValueError as exc:
            starts.append(s_val)
            ends.append(e_val)
            errors.append((rownum, accession, str(exc), s_val, e_val))

    result = pro.copy()
    if errors:
        return result, pd.DataFrame(errors, columns=["excel_row", "accessions", "reason", "old_start", "old_end"]), (changed, multi, no_change)
    result["start"] = starts
    result["end"] = ends
    result.drop(columns=["first_boundary", "second_boundary"], inplace=True)
    return result, pd.DataFrame(columns=["excel_row", "accessions", "reason", "old_start", "old_end"]), (changed, multi, no_change)


def classify(pro, pep):
    for label, df, needed in (
        ("Annotated_Proteins", pro, {"type", "proteins"}),
        ("Annotated_Peptides", pep, {"type", "cds_type", "proteins"}),
    ):
        missing = sorted(needed - set(df.columns))
        if missing:
            raise ValueError(f"{label} is missing classification fields: {missing}")
    pro = pro.copy()
    pep = pep.copy()
    is_cp_pro = pro["type"].eq("CPs") & pd.to_numeric(pro["proteins"], errors="coerce").eq(1)
    pro["Class"] = "excluded"
    pro.loc[is_cp_pro, "Class"] = "CP"

    protein_one = pd.to_numeric(pep["proteins"], errors="coerce").eq(1)
    cp_mask = pep["type"].eq("CDS") & pep["cds_type"].eq("CDS(verified by protein_db)") & protein_one
    ncp_mask = (
        (pep["type"].isin(NON_CDS) |
         (pep["type"].eq("CDS") & pep["cds_type"].eq("CDS(out of frame)")))
        & protein_one
    )
    if (cp_mask & ncp_mask).any():
        raise ValueError("CP/NCP classification rules overlap; check the input data")
    pep["Class"] = "excluded"
    pep.loc[cp_mask, "Class"] = "CP"
    pep.loc[ncp_mask, "Class"] = "NCP"
    return pro, pep


DEDUP_COLS = ["material", "sequence", "chrom", "start", "end", "strand"]


def normalized_dedup_keys(df):
    """Normalize coordinate strings for deduplication without altering output values."""
    missing = [c for c in DEDUP_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing mandatory deduplication columns: {missing}. Material and genomic coordinates must not be omitted from deduplication.")
    keys = pd.DataFrame(index=df.index)
    for col in DEDUP_COLS:
        v = df[col].astype("string").str.strip()
        v = v.replace("", pd.NA)
        if col in {"start", "end"}:

            def coordinate_key(x):
                if pd.isna(x):
                    return pd.NA
                parts = str(x).split("|")
                try:
                    nums = [int(float(z.strip())) for z in parts]
                    if any(float(z.strip()) != n for z, n in zip(parts, nums)):
                        return str(x)
                    return "|".join(map(str, nums))
                except (ValueError, OverflowError):
                    return str(x)
            v = v.map(coordinate_key)
        keys[col] = v
    return keys


def distinct_join(values):
    """Combine unique semicolon-separated identifiers in first-seen order."""
    seen = {}
    for value in values:
        if pd.isna(value):
            continue
        for item in str(value).split(";"):
            item = item.strip()
            if item:
                seen.setdefault(item, None)
    return ";".join(seen)


def deduplicate_cp(cp_all):
    """Deduplicate CP observations by tissue and position, keeping the strongest PSM record."""
    keys = normalized_dedup_keys(cp_all)
    complete = keys.notna().all(axis=1)
    invalid_count = int((~complete).sum())
    valid = cp_all.loc[complete].copy()
    valid["__row_order"] = valid.index
    valid["__numeric_psms"] = pd.to_numeric(valid["psms"], errors="coerce").fillna(-1) if "psms" in valid else -1
    internal_keys = [f"__dedup_{col}" for col in DEDUP_COLS]
    for col, internal in zip(DEDUP_COLS, internal_keys):
        valid[internal] = keys.loc[complete, col].to_numpy()

    duplicate_mask = valid.duplicated(internal_keys, keep=False)
    duplicates = valid.loc[duplicate_mask].copy()

    if not duplicates.empty:
        agg_specs = {"Source": ("Source", distinct_join),
                     "merged_records": ("Source", "size")}
        if "accessions" in duplicates.columns:
            agg_specs["accessions"] = ("accessions", distinct_join)
        aggregates = duplicates.groupby(internal_keys, sort=False).agg(**agg_specs).reset_index()
    else:
        aggregates = None

    valid = (valid.sort_values(["__numeric_psms", "__row_order"],
                               ascending=[False, True], kind="stable")
                  .drop_duplicates(internal_keys, keep="first"))
    valid["merged_records"] = 1
    if aggregates is not None:
        valid = valid.merge(aggregates, on=internal_keys, how="left", suffixes=("", "__merged"), sort=False)
        has_duplicate = valid["merged_records__merged"].notna()
        valid.loc[has_duplicate, "Source"] = valid.loc[has_duplicate, "Source__merged"]
        valid.loc[has_duplicate, "merged_records"] = valid.loc[has_duplicate, "merged_records__merged"].astype(int)
        if "accessions__merged" in valid.columns:
            valid.loc[has_duplicate, "accessions"] = valid.loc[has_duplicate, "accessions__merged"]
        valid = valid.drop(columns=[c for c in valid.columns if c.endswith("__merged")])

    missing_key_rows = cp_all.loc[~complete].copy()
    missing_key_rows["merged_records"] = 1  
    valid = valid.drop(columns=internal_keys + ["__row_order", "__numeric_psms"])
    result = pd.concat([valid, missing_key_rows], ignore_index=False, sort=False)
    result = result.sort_index(kind="stable").reset_index(drop=True)
    result["merged_records"] = result["merged_records"].astype(int)
    stats = {
        "before": len(cp_all), "after": len(result),
        "removed": len(cp_all) - len(result),
        "duplicate_groups": len(aggregates) if aggregates is not None else 0,
        "missing_key_rows": invalid_count,
    }
    return result, stats


def write_workbooks(outdir, pro, pep):
    """Write deduplicated CP and unmodified NCP workbooks."""
    if "Source" in pro.columns or "Source" in pep.columns:
        raise ValueError("The input already contains a Source column; rename it to avoid overwriting provenance")

    cp_pro = pro.loc[pro["Class"].eq("CP")].copy()
    cp_pep = pep.loc[pep["Class"].eq("CP")].copy()
    ncp_pep = pep.loc[pep["Class"].eq("NCP")].copy()
    cp_pro["Source"] = "Annotated_Proteins"
    cp_pep["Source"] = "Annotated_Peptides"

    cp_all = pd.concat([cp_pro, cp_pep], axis=0, ignore_index=True, sort=False)
    cp_unique, cp_stats = deduplicate_cp(cp_all)
    ncp_keys = normalized_dedup_keys(ncp_pep)
    ncp_complete = ncp_keys.notna().all(axis=1)
    ncp_dup = int(ncp_keys.loc[ncp_complete].duplicated(DEDUP_COLS).sum())

    cp_path = outdir / "CP_peptides.xlsx"
    ncp_path = outdir / "NCP_peptides.xlsx"
    with pd.ExcelWriter(cp_path, engine="openpyxl") as writer:
        cp_unique.to_excel(writer, sheet_name="CP", index=False)
    with pd.ExcelWriter(ncp_path, engine="openpyxl") as writer:
        ncp_pep.to_excel(writer, sheet_name="NCP", index=False)
    return cp_path, ncp_path, len(cp_pro), len(cp_pep), len(ncp_pep), cp_stats, ncp_dup


def main():
    args = parse_args()
    source = Path(args.input).expanduser().resolve()
    fasta = Path(args.protein_fa).expanduser().resolve()
    outdir = Path(args.output_dir).expanduser().resolve() if args.output_dir else source.parent
    if not source.is_file() or not fasta.is_file():
        raise FileNotFoundError(f"Missing input file: Excel={source} (exists={source.is_file()}); "
                                f"FASTA={fasta} (exists={fasta.is_file()})")
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Input Excel: {source}", flush=True)
    print(f"Reference protein: {fasta}", flush=True)
    print("[1/3] Reading Excel and reference protein CDS blocks...", flush=True)
    pro = pd.read_excel(source, sheet_name="Annotated_Proteins")
    pep = pd.read_excel(source, sheet_name="Annotated_Peptides")
    blocks = load_cds_blocks(fasta)
    print(f"  Annotated_Proteins {len(pro)} rows; Annotated_Peptides {len(pep)} rows; protein CDS entries: {len(blocks)}", flush=True)

    print("[2/3] Integrating exon-spanning coordinates in Annotated_Proteins...", flush=True)
    pro, issues, (n_spliced, n_multi, n_plain) = merge_coding_coordinates(pro, blocks)
    if not issues.empty:
        report = outdir / "boundary_integration_errors.tsv"
        issues.to_csv(report, sep="\t", index=False, encoding="utf-8-sig")
        print(issues.head(10).to_string(index=False), file=sys.stderr)
        raise RuntimeError(f"{len(issues)} rows failed CDS coordinate validation; error report: {report}."
                           "CP/NCP files were not generated; review the report and rerun.")
    print(f"  Split {n_spliced} rows (including {n_multi} spanning at least two introns); {n_plain} unchanged or unlocalized rows", flush=True)

    print("[3/3] Classifying and exporting using the original rules...", flush=True)
    pro, pep = classify(pro, pep)
    cp, ncp, n_cp_pro, n_cp_pep, n_ncp, cp_stats, ncp_dup = write_workbooks(outdir, pro, pep)
    n_ex_pro = int(pro["Class"].eq("excluded").sum())
    n_ex_pep = int(pep["Class"].eq("excluded").sum())
    print(f"  CP: before deduplication: {cp_stats['before']} rows; deduplication keys: material+sequence+chrom+start+end+strand; "
          f"after deduplication: {cp_stats['after']} rows; duplicates removed: {cp_stats['removed']} rows; "
          f"duplicate groups: {cp_stats['duplicate_groups']}; {cp}")
    if cp_stats["missing_key_rows"]:
        print(f"  [Notice] CP has {cp_stats['missing_key_rows']} rows with missing deduplication fields; retained without deduplication")
    print(f"  NCP: not deduplicated; total {n_ncp} rows; the six-key definition identifies {ncp_dup} duplicate records; {ncp}")
    print(f"  excluded rows not exported: Annotated_Proteins {n_ex_pro} rows; Annotated_Peptides {n_ex_pep} rows")
    print("  The source Excel was not modified; exon-spanning CP coordinates use '|' to separate CDS blocks.")


if __name__ == "__main__":
    main()
