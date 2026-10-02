"""Merge per-sample gene count files and build the novel-gene expression matrix.

Steps: (1) outer-join the per-sample featureCounts outputs listed in the RNA-seq info
workbook; (2) keep genes whose sequence contains no stop codon (CDS FASTA pre-filtered
by 02_map_ms_peptides_to_transcripts.py --filter-stop-codon); (3) keep expressed genes
(nonzero-mean >= 5 and expressed in at least one sample); (4) TPM-normalise; (5) classify
each gene by the genomic elements it spans.

Inputs:  --count-dir (per-sample *_counts.txt), --rna-info (xlsx, Sheet2 with Sample column),
         --cds-fasta (CDS FASTA with chrom/strand/start/end in FASTA descriptions),
         --element-file (aa-pipeline results xlsx, sheet "Genomic_Features"), --out-dir
Outputs: <out-dir>/total_expressed.csv, total_expressed_tpm.csv, total_expression_matrix.csv
Usage:
  py 03_merge_count_matrices.py --count-dir ./data/processed/counts \
      --rna-info ./data/raw/Total_rna_seq.xlsx \
      --cds-fasta ./output/11_transcript_model_construction/codon_filter_gtf_cds.fasta \
      --element-file ./data/raw/peptide_analysis_results.xlsx \
      --out-dir ./output/11_transcript_model_construction
"""
import argparse
import pandas as pd
import os
from tqdm import tqdm
import numpy as np
from Bio import SeqIO
import re

def merge_count_files(input_dir, rna_info_file, gene_id_col="Geneid"):
    count_files = pd.read_excel(rna_info_file, sheet_name="Sheet2")
    merged_df = None
    for _, row in tqdm(count_files.iterrows(), desc="Merging"):
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
            print(f"\nFailed on {file}: {str(e)}")
            continue
    if merged_df is not None:
        print(f"\nMerged matrix shape: {merged_df.shape}")
        return merged_df
    else:
        raise ValueError("No count files were merged")

def filter_codon_gene(expression_matrix, cds_file):
    kept_id = []
    for rec in SeqIO.parse(cds_file, "fasta"):
        kept_id.append(re.sub(r'\.[^.]*$', '', rec.id))
    return expression_matrix[expression_matrix['GeneID'].isin(kept_id)]

def filter_expressed_genes(count_df, output_prefix, mean_threshold=5):
    expr_matrix = count_df.drop(columns=['GeneID'])
    nonzero_mean = expr_matrix.replace(0, np.nan).mean(axis=1, skipna=True)
    has_expression = (expr_matrix > 0).any(axis=1)
    condition = (nonzero_mean >= mean_threshold) & (has_expression)
    expressed_genes = count_df[condition]
    expressed_file = f"{output_prefix}_expressed.csv"
    expressed_genes.to_csv(expressed_file, index=False)
    print("\nGene expression filter:")
    print(f"Total genes: {len(count_df)}")
    print(f"Expressed genes (nonzero mean >= {mean_threshold}, expressed in >=1 sample): {len(expressed_genes)}")
    print(f"Saved to: {expressed_file}")
    return expressed_genes

def prepare_length_data(cds_file):
    gene_lengths = {}
    for rec in SeqIO.parse(cds_file, "fasta"):
        _, _, start, end = rec.description.split('\t')
        rec_id = re.sub(r'\.[^.]*$', '', rec.id)
        gene_lengths[rec_id] = int(end) - int(start) + 1
    gene_length_df = pd.DataFrame(list(gene_lengths.items()), columns=['GeneID', 'length'])
    return gene_length_df

def normalize_tpm(count_df, length_df, output_file):
    df = pd.merge(count_df, length_df, on='GeneID', how='inner')
    df = df[df['length'] > 0]
    sample_cols = [col for col in df.columns if col not in ['GeneID', 'length']]
    tpm_data = {}
    for sample in sample_cols:
        rpk = (df[sample] * 10**3) / df['length']
        per_million_scaling_factor = rpk.sum() / 10**6
        tpm = rpk / per_million_scaling_factor
        tpm_data[sample] = tpm
    tpm_df = pd.concat([df[['GeneID']], pd.DataFrame(tpm_data)], axis=1)
    tpm_df.to_csv(output_file, index=False)
    print(f"TPM normalisation done: {output_file} ({len(tpm_df)} rows)")

def extract_gene_locations(expression_df, cds_file):
    gene_info = {}
    for rec in SeqIO.parse(cds_file, "fasta"):
        chrom_temp, strand, start, end = rec.description.split('\t')
        chrom = chrom_temp.split(" ")[-1]
        rec_id = re.sub(r'\.[^.]*$', '', rec.id)
        if rec_id in expression_df['GeneID'].values:
            gene_info[rec_id] = (chrom, strand, start, end)
    return gene_info

def gene_location(expression_df, element_file, cds_file, output_file):
    gene_locations = extract_gene_locations(expression_df, cds_file)
    element_df = pd.read_excel(element_file, sheet_name="Genomic_Features")
    results = []
    for gene_id, gene_info in gene_locations.items():
        chrom, strand, start, end = gene_info
        chrom_elements = element_df[(element_df['chrom'] == chrom) & (element_df['strand'] == strand)].sort_values('start').reset_index(drop=True)
        overlapping_elements = []
        for _, elem in chrom_elements.iterrows():
            if not (elem['end'] < int(start) or elem['start'] > int(end)):
                overlapping_elements.append(elem['type'])
        unique_elements = []
        for elem_type in overlapping_elements:
            if elem_type not in unique_elements:
                unique_elements.append(elem_type)
        elements_string = "-".join(unique_elements) if unique_elements else None
        results.append({
            'gene_id': gene_id,
            'chrom': chrom,
            'strand': strand,
            'start': start,
            'end': end,
            'spanning_elements': elements_string,
            'element_count': len(unique_elements)
        })
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_file, index=False)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--count-dir", required=True, help="directory of per-sample *_counts.txt files")
    ap.add_argument("--rna-info", required=True, help="RNA-seq info xlsx (Sheet2 lists sample names)")
    ap.add_argument("--cds-fasta", required=True, help="stop-codon-free CDS FASTA (from 02 --filter-stop-codon)")
    ap.add_argument("--element-file", required=True, help="aa-pipeline results xlsx (sheet 'Genomic_Features')")
    ap.add_argument("--out-dir", required=True, help="output directory")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print("=== Step 1/4: merge count files ===")
    expression_matrix = merge_count_files(args.count_dir, args.rna_info)
    # keep genes whose sequence contains no in-frame stop codon
    expression_matrix_filter = filter_codon_gene(expression_matrix, args.cds_fasta)
    print("\n=== Step 2/4: expression filter ===")
    expression_genes = filter_expressed_genes(expression_matrix_filter, os.path.join(args.out_dir, "total"))
    print("\n=== Step 3/4: TPM normalisation ===")
    length_df = prepare_length_data(args.cds_fasta)
    normalize_tpm(
        expression_genes,
        length_df,
        os.path.join(args.out_dir, "total_expressed_tpm.csv")
    )
    print("\n=== Step 4/4: gene location classification ===")
    gene_location(expression_genes, args.element_file, args.cds_fasta,
                  os.path.join(args.out_dir, "total_expression_matrix.csv"))


if __name__ == "__main__":
    main()
