# Domain-Specific Categories: Plant Small Peptides Identification (2016–2029)

This file provides a domain-specific taxonomy for writing a review on **plant small peptide identification**,
covering peptides from **canonical coding regions** and **non-canonical ORFs** (smORFs in non-coding contexts),
with a focus on **evidence chains** and **method comparison**.

---

## 1) Definition & Scope (Recommended for Section 2)

### 1.1 What counts as a "plant small peptide"?
Common practical scopes (choose one and state clearly in the manuscript):
- **<100 aa** small peptides / microproteins
- **<150 aa** short proteins (more inclusive)

### 1.2 Key origin categories
1. **Canonical coding-region peptides**
   - annotated short proteins
   - secreted peptide precursor proteins (processed into mature peptides)
2. **Non-canonical ORF peptides (smORF-derived)**
   - uORFs (5'UTR)
   - dORFs (3'UTR)
   - lncRNA-encoded smORFs
   - intergenic transcripts containing smORFs
   - alternative ORFs within annotated genes (overlapping ORFs)

---

## 2) Evidence-Level Framework (Core Comparison Axis)

Use this unified evidence framework throughout the review:

- **L1: Computational prediction**
  - ORF prediction, conservation, motif or signal peptide predictions
- **L2: Transcription evidence**
  - RNA-seq / RT-qPCR supports transcript existence
- **L3: Translation evidence**
  - ribosome profiling (Ribo-seq), initiation site profiling, ribosome occupancy
- **L4: Peptide detection evidence**
  - LC–MS/MS proteomics or peptidomics detects peptide/fragment
- **L5: Functional validation**
  - genetics (KO/CRISPR/OE), synthetic peptide assays, rescue experiments, receptor binding evidence

---

## 3) Identification & Discovery Pipelines (Core Methods Section)

### 3.1 Computational discovery (L1)
1. ORF prediction from genome/transcriptome
2. length threshold filtering (e.g., <100 aa)
3. conservation-based prioritization (multi-species alignment)
4. motif scanning (family motifs, cleavage motifs)
5. signal peptide prediction (secreted peptide candidate enrichment)
6. subcellular localization prediction (optional)
7. secondary structure/disorder prediction (optional)

**Typical outputs:**
- candidate lists (low-to-medium confidence)
- prioritization scores

**Common limitations:**
- high false positive rate (translation not guaranteed)
- lineage-specific smORFs may lack conservation signal
- annotation incompleteness in some crops

---

### 3.2 Coding-region peptide identification (Canonical ORFs)
1. mining annotated genomes for short ORFs
2. identifying peptide precursor genes (prepropeptides)
3. transcript evidence across tissues/conditions (L2)
4. proteomics confirmation (L4)
5. functional validation (L5)

**Key practical note:**
- Many secreted peptides undergo processing; detected mature peptides may not directly match precursor length.

---

### 3.3 Non-coding region / non-canonical ORF identification (smORFs)
This is the core discovery frontier for "hidden coding potential".

#### 3.3.1 Translation-based discovery (L3)
1. ribosome occupancy patterns supporting translation
2. initiation site profiling / start-site capture approaches
3. disambiguating scanning noise vs real translation
4. reproducibility checks (datasets, tissues, conditions)

**Typical outputs:**
- translated smORF catalogs
- inferred translation efficiency trends

**Limitations:**
- mapping ambiguity for short ORFs
- condition-specific translation may be missed
- hard to resolve exact start codon without specialized profiling

#### 3.3.2 Proteomics / peptidomics-based discovery (L4)
1. shotgun proteomics detection of short proteins
2. peptidomics enrichment for secreted/processed peptides
3. targeted MS validation of candidate peptides

**Limitations:**
- low abundance peptides often undetectable
- post-translational modifications (PTMs) complicate search
- sample preparation bias (tissue, extraction, enrichment)

#### 3.3.3 Integrated multi-omics discovery (L1–L4)
Recommended integrated pipeline (best practice):
- transcript evidence (L2) + translation evidence (L3) + peptide evidence (L4)
- then prioritize for L5 functional validation

---

## 4) Translation Initiation & Start Codon Diversity (Supporting Section)

### 4.1 Start codon categories to discuss
1. **AUG (canonical start)**
2. **Near-cognate start codons**
   - examples: CUG, GUG, UUG (and others depending on literature evidence)

### 4.2 Factors affecting initiation and translation capacity
- Kozak-like initiation context (plant-specific preferences)
- 5'UTR structure and upstream AUGs
- uORF-mediated translation regulation
- alternative transcription start sites (TSS)
- alternative splicing changes ORF boundaries
- tissue- and stress-dependent translation remodeling

### 4.3 Implications for identification bias
- annotation pipelines often assume AUG starts → may miss non-AUG initiated peptides
- proteomics pipelines may miss low-translation candidates
- Ribo-seq helps but has false-positive risks

---

## 5) Functional Characterization (Downstream Validation Section)

### 5.1 Functional categories (What small peptides do)
1. Development and morphogenesis
2. Immune regulation and pathogen response
3. Abiotic stress response (salt, drought, temperature)
4. Nutrient signaling and metabolic control
5. Cell-to-cell communication (secreted peptide signaling)
6. Reproductive development (optional, depending on coverage)

### 5.2 Mechanistic modes (How peptides act)
- secreted peptide → receptor-mediated signaling (RLK/RLP systems)
- intracellular peptides → protein interaction / pathway modulation
- peptide processing and maturation as functional switch

### 5.3 Functional evidence checklist (L5)
- mutant phenotype (KO/CRISPR)
- overexpression / knockdown
- synthetic peptide application (dose + tissue specificity)
- rescue experiments
- binding assays / receptor interaction evidence (if applicable)

**Common limitations:**
- functional redundancy and compensation
- non-physiological peptide concentrations in treatments
- context-specific phenotypes (tissue/condition dependence)

---

## 6) Databases & Resources (Make a Dedicated Table)

Resource categories to include in the review:
1. Plant genome annotation portals (Arabidopsis + crops)
2. Transcriptomic expression resources (tissue/stress atlases)
3. Ribo-seq datasets and browsers
4. Proteomics / peptidomics repositories
5. Small peptide databases (if available)
6. Tools for ORF prediction and peptide annotation

Recommended outputs for the manuscript:
- **Table: Resources for plant small peptide discovery**
  - Name | Data type | Species coverage | What it supports | Notes/limitations

---

## 7) Suggested Comparison Tables (Directly reusable)

### Table A: Identification method comparison (core)
| Method | Input | Output | Typical Evidence Level | Strengths | Limitations | Best Use Case |

### Table B: Coding vs non-coding origin peptide discovery
| Origin | Candidate Features | Discovery Strategy | Evidence Bottleneck | Typical Validation Route |

### Table C: Start codons and translation evidence
| Start codon | Common context | Translation evidence | Expected detectability | Notes |

### Table D: Case studies (representative peptides)
| Candidate/Peptide | Species | Origin | Evidence Level | Key assays | Functional category |
