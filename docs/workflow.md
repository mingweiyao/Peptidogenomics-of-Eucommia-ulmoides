# Workflow and module dependencies

This document describes the recommended execution order, the inputs and outputs of each
module, and how the Python, R and Shell parts connect. Module and script numbering
(`NN_`) expresses the recommended order; several modules can be run independently once
their inputs exist.

```
00_database ──► 01_gff_cleanup ──► 02_feature_filtering ──► 03_rnaseq_quantification
                                                                     │
                                                                     ▼
               09_conservation   10_n_to_c_substitution   04_expression_analysis
                                             │
                     05_terminal_enrichment, 06_reading_frame (feature-based analyses)
                                             │
                     07_stress_response (R WGCNA ─► candidate screen)
                     08_rubber_biosynthesis (R WGCNA ─► candidate screen)
```

## 00_database — six-frame peptide database

| | |
|---|---|
| Scripts | `01_build_six_frame_peptide_database.py`, `02_verify_input_checksums.py` |
| Inputs | genome FASTA, GFF, six-frame peptide FASTA databases, LC–MS/MS search result tables |
| Outputs | `peptide_analysis_results.xlsx` (per-peptide classification, coordinates, properties) |
| Notes | script 01 uses module-level CONFIG paths (edit before running, or set `AA_RAW_DIR` / `AA_OUT_DIR`); script 02 recomputes SHA-256 checksums of key inputs and compares them against the archived reference values |

## 01_gff_cleanup — coordinate curation

| | |
|---|---|
| Script | `01_clean_gff_coordinates.py` |
| Inputs | source caches produced by `02_feature_filtering/01_prepare_source_tables.py` |
| Outputs | `_base_data_chain/gff_223870.pkl` / `.csv` (curated non-overlapping features) |
| Notes | deduplication → remove fully-contained intervals → junction/CDS overlap removal → overlap filtering; requires `intervaltree` |

## 02_feature_filtering — tables and high-confidence NCPs

| | |
|---|---|
| Scripts | `01_prepare_source_tables.py` → `02_build_cp_ncp_tables.py` → `03_filter_high_confidence_ncps.py` → `04_compute_peptide_properties.py` → `05_flag_low_complexity_sequences.py` |
| Inputs | `peptide_analysis_results.xlsx` (01); source caches (02, 03); CV/CPM expression matrices and sample groups (03) |
| Outputs | supplementary-table-format CP and NCP tables; `_base_data_chain/supplementary_table2_NCP.csv`; transcript-supported high-confidence NCPs (`supplementary_table4_high_confidence.csv`, then `..._with_recomputed_props.csv`); low-complexity flags |
| Notes | 01 converts all sheets of the aa-pipeline output to pkl caches and verifies the file SHA-256; 03 applies the groupwise CV filter followed by the CPM filter and intersects with unique NCPs; 05 flags candidates with Shannon < 2.5 bits or top-residue fraction ≥ 50% |

## 03_rnaseq_quantification — RNA-seq

| | |
|---|---|
| Scripts | `01_build_sample_metadata.py` (Python), `02_align_reads_and_count_features.sh` (Shell) |
| Inputs | sample/group spreadsheet; raw reads (FASTA/FASTQ), genome + GFF for the Shell pipeline |
| Outputs | sample metadata table; per-sample gene count files |
| Notes | the Shell pipeline runs QC (FastQC), trimming (fastp), alignment (STAR) and gene-level counting (featureCounts), with read-length-aware index building |

## 04_expression_analysis — expression support and SMD

| | |
|---|---|
| Script | `01_analyze_expression_support_and_smd.py` |
| Inputs | high-confidence NCP table (02 module); universe CPM matrix; sample groups; archived Venn workbook and SMD workbooks for cross-checks |
| Outputs | `output/04_expression_analysis/`: MS support counts, transcript-support Venn reconciliation, detection breadth, recomputed SMD statistics |

## 05_terminal_enrichment — terminal enrichment under five null models

| | |
|---|---|
| Scripts | `00_run_terminal_enrichment.py` runs `01`–`05` (null models) then `06` (figure) |
| Inputs | high-confidence NCP table; genome FASTA; GFF; chromosome lengths and scaffold↔chromosome mapping |
| Outputs | `output/05_terminal_enrichment/`: one CSV per null model (baseline permutation, non-coding-restricted, GC-exact, GC-matched, expression-restricted) + `SuppFig5_terminal_enrichment.png/pdf` |
| Notes | 2,000 within-chromosome, length-preserving permutations per model; empirical p = (≥ obs + 1)/2001 |

## 06_reading_frame — six-frame ORF distribution

| | |
|---|---|
| Script | `01_test_reading_frame_distribution.py` |
| Inputs | CP/NCP table (02 module); high-confidence NCP table; chromosome lengths; chromosome mapping |
| Outputs | `output/06_reading_frame/`: frame × type counts; frame-matched null test |

## 07_stress_response — drought/salt module

| | |
|---|---|
| Scripts | `01_run_wgcna_network.R` → `02_scan_wgcna_stability.R` (R); `03_screen_stress_candidates.py`, `04_validate_wgcna_results.py` (Python) |
| Inputs | drought–salt expression matrix and trait table (R); DESeq2 result tables, module lists (Python) |
| Outputs | WGCNA network, module assignments, hub genes, stability grid; stress candidate tables; reconstruction report validating archived module/candidate lists |
| Notes | signed network; soft power 16, minModuleSize 200, mergeCutHeight 0.4, deepSplit 0; needs ≥ 16 GB RAM |

## 08_rubber_biosynthesis — rubber module

| | |
|---|---|
| Scripts | `01_run_wgcna_network.R` → `02_scan_wgcna_stability.R` → `03_scan_tom_threshold.R` (R); `04_screen_rubber_candidates.py` (Python) |
| Inputs | rubber expression matrix and trait table (R); rubber expression matrix, gene list, high-TOM pair references (Python) |
| Outputs | WGCNA network, stability and TOM-threshold scans; NCP–gene association tables with Pearson verification and tissue-adjusted partial correlations |
| Notes | signed network; soft power 22, minModuleSize 150, mergeCutHeight 0.3, deepSplit 2; needs ≥ 32 GB RAM |

## 09_conservation — cross-species screen

| | |
|---|---|
| Script | `01_screen_kmer_conservation.py` |
| Inputs | NCP protein sequences; annotated *E. ulmoides* proteome; UniProt reference proteomes (see `input_data.md`) |
| Outputs | `conservation_kmer_screen.csv` (exact 8-mer match rates per reference) |

## 10_n_to_c_substitution — N→C quantification

| | |
|---|---|
| Script | `01_quantify_n_to_c_substitutions.py` |
| Inputs | genome FASTA |
| Outputs | N-region table; per-region counts; peptides whose translation requires N→C substitution |
