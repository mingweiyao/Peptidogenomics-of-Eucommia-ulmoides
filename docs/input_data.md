# Input data and reference resources

This file lists every external dataset used by the workflow: what it is, where it comes
from, and how to verify it. Place the files in any local layout and pass the locations to
the scripts through their CLI arguments (see the main README).

## Genome and annotation

| Item | Value |
|---|---|
| Species | *Eucommia ulmoides* (EU17 assembly) |
| Genome accession | **GWHBISF00000000** (NGDC Genome Warehouse, https://ngdc.cncb.ac.cn/gwh/) |
| Genome FASTA | `Eu_genome.fasta` — 957,359,030 bytes |
| SHA-256 | `ab26acede9214c0f1a06ed6778e1a2b2a7b65f85175a8ad2ae450a67e24ae0b6` |
| Annotation GFF | `GWHBISF00000000.gff` — 56,436,010 bytes |
| SHA-256 | `b81b8884e574991d7798c18e060acc18cc6f9db09870511eb6f5cc4dfa01c544` |

Chromosome metadata used by the distribution analyses (chromosome lengths and the
scaffold-ID ↔ chromosome-number mapping) are derived from this assembly.

## RNA-seq

- 122 transcriptome samples (41 tissue/group combinations), NCBI SRA.
- SRA Run accessions for every sample are listed in the sample metadata table built by
  `scripts/03_rnaseq_quantification/01_build_sample_metadata.py`.

## Reference proteomes (conservation screen)

Downloaded from UniProt (https://www.uniprot.org/) in September 2026:

| Species | UniProt proteome accession |
|---|---|
| *Arabidopsis thaliana* | UP000006548 |
| *Populus trichocarpa* | UP000006729 |
| *Vitis vinifera* | UP001227230 |

The annotated *E. ulmoides* proteome used as a same-species control is derived from the
GFF annotation above.

## LC–MS/MS search outputs

| Item | Value |
|---|---|
| Six-frame peptide database (searched) | `Eu_peptide_database_customized_5.fa` — 7,285,964,522 bytes |
| SHA-256 | `ec636aa0952748776a502723e74e4626d5f82b5d9ebc6205f5c54307936192d5` |
| Peptide/protein search result workbook | `peptide_analysis_results.xlsx` — 97,764,606 bytes |
| SHA-256 | `16a365f9bf02056e392bbc8d8089d10058c39cb8d29cbb10a1ecbce34ebc781d` |

`scripts/00_database/02_verify_input_checksums.py` recomputes these checksums and reports
whether the local copies match; identical checksums guarantee that downstream scripts
reproduce the published numbers.

## Rubber-biosynthesis gene set

A curated list of 52 rubber-biosynthesis genes (MVA/MEP pathway, initiation and extension
enzymes) is used with the rubber expression matrix in module 08. The list is distributed
with the expression matrix input; gene symbols follow the manuscript annotation.

## Expression matrices

Two sample sets are used downstream:

| Set | Samples | Matrix | Used by |
|---|---|---|---|
| Stress (drought/salt) | 14 | stress expression matrix + trait table | module 07 |
| Rubber | 15 | rubber FPKM matrix + gene annotations | module 08 |

Both matrices are gene × sample tables (FPKM); the stress trait table encodes
control/drought/salt conditions, the rubber trait table encodes latex-yield groups.
