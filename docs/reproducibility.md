# Reproducibility

## Software

| Component | Version | Used by |
|---|---|---|
| Python | ≥ 3.9 | all `scripts/` modules except the WGCNA network steps |
| R | ≥ 4.2 | modules 07, 08 (WGCNA) |
| WGCNA | CRAN release compatible with R ≥ 4.2 | network construction, TOM, stability scans |
| reshape2, dplyr, stringr | CRAN | data handling in R scripts |
| DESeq2 | Bioconductor | stress differential expression (upstream of module 07) |
| STAR | recent release (a `STAR-avx2` build works identically) | alignment (module 03) |
| fastp, FastQC | recent releases | QC and trimming (module 03) |
| samtools | ≥ 1.9 | BAM filtering/indexing (module 03) |
| featureCounts (Subread) | ≥ 2.0 | gene-level counting (module 03) |
| bioawk | any | read-length detection (module 03) |
| EMBOSS (`sixpack`) | ≥ 6.6 | six-frame translation (upstream of module 00) |

`requirements.txt` lists the Python dependencies and their minimum tested versions (pip). For a single-step setup use
`environment.yml` (conda: Python + R + Bioconda tools). WGCNA network steps are memory
intensive: ≥ 16 GB RAM for the stress module (14 samples), ≥ 32 GB RAM for the rubber
module (15 samples). The feature-filtering module reads CSV matrices of several hundred MB
and needs a few GB of RAM; reading the six-frame peptide database (≈ 7 GB FASTA) requires
corresponding disk space.

## Randomness

All permutation-based null models use a fixed seed, so runs are bit-for-bit reproducible:

| Analysis | Seed | Permutations |
|---|---|---|
| Terminal enrichment, five null models (`05_terminal_enrichment/01`–`05`) | 20260906 | 2,000 per model |
| Reading-frame matched null (`06_reading_frame/01`) | 20260909 | 2,000 |

Empirical p-values are computed as (number of permutations ≥ observed + 1)/2001.

## Main commands

```bash
# environment
conda env create -f environment.yml && conda activate eucommia-ncp

# database + integrity check
py scripts/00_database/02_verify_input_checksums.py --out ./output --rej ./data/raw

# curation and filtering chain
py scripts/02_feature_filtering/01_prepare_source_tables.py --peptide_xlsx ./data/raw/peptide_analysis_results.xlsx --out ./output
py scripts/02_feature_filtering/02_build_cp_ncp_tables.py --out ./output
py scripts/02_feature_filtering/03_filter_high_confidence_ncps.py --out ./output --rej ./data/processed

# feature-based analyses
py scripts/05_terminal_enrichment/00_run_terminal_enrichment.py --out ./output --rej ./data/processed --gen13 ./data/reference
py scripts/06_reading_frame/01_test_reading_frame_distribution.py --out ./output --rej ./data/reference --gen13 ./data/reference

# WGCNA networks and stability scans (R)
Rscript scripts/07_stress_response/01_run_wgcna_network.R     # env: REPRO_REJ, REPRO_OUT
Rscript scripts/07_stress_response/02_scan_wgcna_stability.R
Rscript scripts/08_rubber_biosynthesis/01_run_wgcna_network.R
Rscript scripts/08_rubber_biosynthesis/02_scan_wgcna_stability.R
Rscript scripts/08_rubber_biosynthesis/03_scan_tom_threshold.R

# downstream screens
py scripts/07_stress_response/03_screen_stress_candidates.py --out ./output --rej ./data/processed --pcr ./data/processed
py scripts/08_rubber_biosynthesis/04_screen_rubber_candidates.py --out ./output --rej ./data/processed --pcr ./data/processed
py scripts/09_conservation/01_screen_kmer_conservation.py --out ./output
py scripts/10_n_to_c_substitution/01_quantify_n_to_c_substitutions.py --out ./output --rej ./data/reference
```

Every script also prints a set of verification targets (expected row counts, anchor
statistics) when it runs, so partial results can be checked against the manuscript
without rerunning the whole chain. R scripts take their roots from the environment
variables `REPRO_REJ` (data root), `REPRO_OUT` (output root) and, where needed,
`REPRO_ROOT` (reference root).

## Line endings

The repository enforces LF for text files via `.gitattributes`. After cloning on Windows,
run `git add --renormalize .` once if git reports line-ending warnings.
