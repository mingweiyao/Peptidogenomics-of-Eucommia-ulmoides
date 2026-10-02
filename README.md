# Eucommia ulmoides NCP workflow

Companion code for a peptidogenomic study of non-conventional peptides (NCPs) in
*Eucommia ulmoides* (manuscript under review). The workflow covers six-frame peptide
database construction, LC–MS/MS peptide classification (conventional vs. non-conventional),
transcriptome-based filtering, and downstream expression, distribution, conservation and
substitution analyses.

## Workflow overview

```
Six-frame genome translation → LC–MS/MS peptide search → CP/NCP classification
→ genomic-coordinate curation → transcriptome-based filtering → transcript-supported NCPs
→ expression & distribution analyses (SMD · reading frame · terminal enrichment)
→ WGCNA + DESeq2 candidate prioritization → tissue-adjusted partial correlation
→ cross-species conservation screen → N→C substitution quantification
```

## Repository structure

```
Eucommia-NCP-workflow/
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt          # Python dependencies
├── environment.yml           # full conda environment (Python + R + command-line tools)
├── .gitignore
├── .gitattributes
├── docs/
│   ├── workflow.md           # execution order, per-module inputs/outputs
│   ├── input_data.md         # external data sources, versions, checksums
│   └── reproducibility.md    # software versions, commands, runtime resources
├── example_data/
│   ├── sample_count_matrix.txt
│   └── sample_genome_491_head.fa
└── scripts/
    ├── 00_database/              # six-frame peptide database + input checksums
    ├── 01_gff_cleanup/           # genomic-coordinate curation
    ├── 02_feature_filtering/     # CP/NCP tables, high-confidence NCPs, properties, low complexity
    ├── 03_rnaseq_quantification/ # sample metadata + STAR/featureCounts pipeline
    ├── 04_expression_analysis/   # expression support, detection breadth, SMD statistics
    ├── 05_terminal_enrichment/   # terminal-10% enrichment under five null models
    ├── 06_reading_frame/         # six-frame ORF distribution + frame-matched null
    ├── 07_stress_response/       # WGCNA network + drought/salt candidate screen
    ├── 08_rubber_biosynthesis/   # WGCNA network + rubber-gene association screen
    ├── 09_conservation/          # cross-species k-mer conservation screen
    └── 10_n_to_c_substitution/   # N→C substitution quantification
```

Scripts are named `NN_snake_case` within each module; `NN` is the recommended execution
order inside that module. Modules are logically ordered but many can be run independently —
see [docs/workflow.md](docs/workflow.md) for the dependency graph and per-module I/O.

## Requirements and installation

- **Python ≥ 3.9**: `pip install -r requirements.txt`
- **R ≥ 4.2**: `WGCNA`, `reshape2`, `dplyr`, `stringr` (CRAN); `DESeq2` (Bioconductor)
- **Command-line tools**: `STAR`, `fastp`, `FastQC`, `samtools`, `featureCounts` (Subread), `EMBOSS` (`sixpack`)

Or create the full environment with conda:

```bash
conda env create -f environment.yml
conda activate eucommia-ncp
```

WGCNA steps require ≥ 16 GB RAM (stress module) or ≥ 32 GB RAM (rubber module).

## Input data and reference resources

Most analysis scripts accept data locations through CLI arguments (`--out` plus
module-specific roots; run any script with `--help` to list them); the database-construction and
RNA-seq scripts require configuring local paths before use. A typical layout:

```
data/
├── raw/           # LC–MS/MS search results, RNA-seq inputs
├── reference/     # genome FASTA, GFF annotation, reference proteomes
└── processed/     # expression matrices, sample metadata, derived tables
output/            # all script outputs (created automatically)
```

External resources (sources, versions, accessions and checksums in
[docs/input_data.md](docs/input_data.md)):

- *E. ulmoides* genome and annotation: accession **GWHBISF00000000** (NGDC Genome Warehouse)
- RNA-seq: 122 samples from NCBI SRA (Run accessions listed in the sample metadata table)
- Reference proteomes: UniProt **UP000006548** (*A. thaliana*), **UP000006729** (*P. trichocarpa*),
  **UP001227230** (*V. vinifera*)

## Reproducing the analyses

Modules are numbered by analysis stage rather than strict execution dependency — follow
[docs/workflow.md](docs/workflow.md) for the actual order; within a module, run scripts in their `NN` order.
Examples:

```bash
# six-frame peptide database (edit CONFIG paths inside the script, then run)
py scripts/00_database/01_build_six_frame_peptide_database.py

# terminal enrichment under five null models (runs all models + figure)
py scripts/05_terminal_enrichment/00_run_terminal_enrichment.py --out ./output --rej ./data/processed

# WGCNA networks (R)
Rscript scripts/07_stress_response/01_run_wgcna_network.R
Rscript scripts/08_rubber_biosynthesis/01_run_wgcna_network.R
```

| Module | Key scripts | Purpose |
|---|---|---|
| 00_database | `01_build_six_frame_peptide_database.py`, `02_verify_input_checksums.py` | six-frame translated peptide database; input integrity check |
| 01_gff_cleanup | `01_clean_gff_coordinates.py` | curated non-overlapping genomic features |
| 02_feature_filtering | `01`–`05` | CP/NCP tables, transcript-supported high-confidence NCPs, peptide properties, low-complexity flags |
| 03_rnaseq_quantification | `01_build_sample_metadata.py`, `02_align_reads_and_count_features.sh` | sample metadata; alignment + gene-level counts |
| 04_expression_analysis | `01_analyze_expression_support_and_smd.py` | tissue/development expression support and SMD statistics |
| 05_terminal_enrichment | `00_run_terminal_enrichment.py` (+ `01`–`06`) | terminal enrichment: baseline + non-coding-, GC-exact-, GC-matched-, expression-restricted nulls |
| 06_reading_frame | `01_test_reading_frame_distribution.py` | six-frame ORF distribution with frame-matched null |
| 07_stress_response | `01`–`04` | stress WGCNA network, stability scan, drought/salt candidate screen, archived-result validation |
| 08_rubber_biosynthesis | `01`–`04` | rubber WGCNA network, stability and TOM-threshold scans, rubber-gene association screen |
| 09_conservation | `01_screen_kmer_conservation.py` | exact k-mer conservation against reference proteomes |
| 10_n_to_c_substitution | `01_quantify_n_to_c_substitutions.py` | peptides created by N→C substitution at genomic N-regions |

## Expected outputs

All scripts write under the `--out` root (default `./output`), one folder per module —
for example `output/05_terminal_enrichment/` (per-model CSV results + figure),
`output/07_stress_response/` (candidate tables), `output/09_conservation/`
(`conservation_kmer_screen.csv`). See [docs/workflow.md](docs/workflow.md) for the full list.

## Data and code availability

All analysis code is contained in this repository. External data sources (genome,
annotation, proteomes, RNA-seq accessions) are listed in
[docs/input_data.md](docs/input_data.md).

## Citation

See [CITATION.cff](CITATION.cff). The associated manuscript is under review; the formal
citation and DOI will be added after publication.

## License

Released under the MIT License — see [LICENSE](LICENSE).
