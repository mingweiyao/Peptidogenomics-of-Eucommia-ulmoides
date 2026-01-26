# Identification of Plant Small Peptides (2016–2029): Strategies from Coding and Non-coding Regions, Evidence Frameworks, and Functional Validation

## Key Points
- Plant small peptides can arise from canonical coding genes and non-canonical smORFs embedded in transcripts previously annotated as non-coding.
- Integrating translation evidence (e.g., ribosome profiling) with proteomics-based peptide detection may reduce false positives in peptide discovery.
- Non-AUG (near-cognate) initiation may contribute to hidden coding potential and influence the detectability of smORF-derived peptides.
- A standardized evidence-level framework (L1–L5) can help compare discovery confidence and prioritize candidates for downstream functional studies.

---

## Abstract
[Write a 200–250 word abstract summarizing:
(1) identification challenge,
(2) major pipelines,
(3) evidence framework contribution,
(4) key gaps and future directions.]

---

## 1. Introduction
### 1.1 Biological importance of plant small peptides
[Discuss signaling, development, immunity, stress response; use hedging language.]

### 1.2 Why identification remains challenging
- Low abundance and spatiotemporal specificity
- Precursor processing and potential PTMs
- Annotation bias against short ORFs and non-AUG initiation

### 1.3 Scope and contributions of this review
- Time span: 2016–2029
- Species: Arabidopsis + crops
- Focus: identification pipelines + evidence framework

---

## 2. Definitions, Taxonomy, and Conceptual Boundaries
### 2.1 What counts as a "small peptide"?
[Define length threshold; explain rationale.]

### 2.2 Canonical vs non-canonical origins
- Canonical coding-region peptides
- smORFs in lncRNAs / UTRs / intergenic transcripts

### 2.3 Evidence levels for peptide identification (L1–L5)
**Table 3. Evidence-level framework for plant small peptide identification**

| Level | Evidence Type | Typical Data | Strength | Common Pitfalls |
|------|---------------|--------------|----------|-----------------|
| L1 | Computational prediction | ORF prediction, conservation | Low–Medium | high false positives |
| L2 | Transcription evidence | RNA-seq, RT-qPCR | Medium | transcript ≠ peptide |
| L3 | Translation evidence | Ribo-seq, initiation profiling | Medium–High | scanning artifacts |
| L4 | Peptide detection | LC–MS/MS, peptidomics | High | low abundance/PTM |
| L5 | Functional validation | genetics, peptide assays | Highest | redundancy, context |

---

## 3. Identification Pipelines: From Candidates to Validated Peptides (Core Section)
### 3.1 Computational discovery of candidate small peptides (L1)
#### 3.1.1 ORF prediction and size thresholds
[Summarize algorithms and heuristics; define sensitivity vs specificity trade-offs.]

#### 3.1.2 Conservation-based prioritization
[Cross-species conservation; Ka/Ks if used; limitations.]

#### 3.1.3 Motif and signal peptide prediction
[Secreted peptides vs intracellular peptides.]

**Limitations:**
[Discuss annotation bias, false positives, missing condition-specific peptides.]

---

### 3.2 Identification from canonical coding regions
#### 3.2.1 Annotated precursor proteins and processing-derived peptides
[Discuss precursor processing concept.]

#### 3.2.2 Transcript evidence and expression specificity (L2)
[RNA-seq, tissue specificity, stress induction patterns.]

**Limitations:**
[Transcript abundance may not correlate with mature peptide abundance.]

---

### 3.3 Identification from non-coding regions and non-canonical ORFs
#### 3.3.1 smORFs in lncRNAs, UTRs, and intergenic transcripts
[Define classes and typical discovery routes.]

#### 3.3.2 Translation evidence via ribosome profiling (L3)
- ribosome occupancy patterns
- initiation site profiling
- distinguishing true translation from scanning

**Limitations:**
[Ribo-seq artifacts, mapping ambiguity, reproducibility across labs.]

---

### 3.4 Proteomics-based peptide detection (L4)
#### 3.4.1 Shotgun proteomics vs targeted strategies
[Summarize coverage and biases.]

#### 3.4.2 Peptidomics for secreted and processed peptides
[Extraction bias, enrichment methods.]

**Limitations:**
[Low abundance peptides, PTMs, sequence redundancy.]

---

### 3.5 Integrative scoring: multi-omics evidence fusion
[Propose a rubric: L1+L2+L3+L4 convergence → high confidence.]

**Table 2. Method comparison table for plant small peptide identification**

| Pipeline/Method | Input | Output | Evidence Level | Strengths | Limitations | Best Use Case |
|----------------|-------|--------|----------------|-----------|-------------|---------------|
| ORF prediction | genome/transcript | candidate list | L1 | high recall | false positives | first-pass scan |
| conservation | orthologs | prioritized candidates | L1 | improves precision | lineage-specific loss | cross-species peptides |
| Ribo-seq | footprints | translated smORFs | L3 | translation signal | scanning artifacts | non-canonical ORFs |
| LC–MS/MS | proteins/peptides | detected peptides | L4 | direct evidence | low abundance | validation stage |
| genetics | mutants/OE | phenotype evidence | L5 | strongest | redundancy | functional proof |

---

## 4. Translation Initiation and Start Codon Diversity (Supporting Section)
### 4.1 AUG initiation and Kozak-like context in plants
[Discuss initiation context; link to detectability.]

### 4.2 Near-cognate start codons and hidden coding potential
[Discuss non-AUG start codons; potential translation efficiency differences.]

### 4.3 Implications for discovery biases
[How annotation and search strategies miss non-AUG peptides.]

---

## 5. Functional Validation and Annotation (Downstream Section)
### 5.1 Functional categories commonly associated with plant small peptides
- Developmental regulation
- Immune responses
- Abiotic stress responses
- Nutrient signaling

### 5.2 Experimental validation strategies (L5)
- genetics (CRISPR KO, OE)
- synthetic peptide application
- receptor binding assays (if relevant)
- rescue experiments

**Limitations:**
[Functional redundancy, context-specific phenotypes.]

**Table 4. Representative plant small peptides and evidence levels (case studies)**

| Peptide/Candidate | Species | Origin (canonical/non-canonical) | Proposed Function | Evidence Level | Key Assays |
|------------------|---------|-----------------------------------|------------------|----------------|-----------|
| [Example] | [Arabidopsis] | [smORF] | [stress response] | [L3/L4/L5] | [Ribo-seq/MS/KO] |

---

## 6. Databases, Tools, and Community Resources
**Table 5. Resources for plant small peptide identification and annotation**

| Resource Type | Name | What It Provides | Best For | Notes |
|--------------|------|------------------|----------|------|
| Genome portal | [ ] | ORF annotation | candidates | |
| Proteomics repo | [ ] | MS datasets | L4 evidence | |
| Ribo-seq repo | [ ] | translation datasets | L3 evidence | |
| Peptide DB | [ ] | curated peptides | reference | |

---

## 7. Challenges, Pitfalls, and Best Practices
### 7.1 Major challenges in discovery
- Low abundance / tissue specificity
- PTMs and processing events
- short sequence mapping ambiguity

### 7.2 Best practices for rigorous identification
- Multi-evidence confirmation (L3 + L4 preferred)
- Clear reporting standards for ORF boundaries and start codons
- Replication across conditions or datasets

---

## 8. Future Directions (2016–2029 and beyond)
- Standardized benchmarking datasets for peptide identification
- Higher-resolution initiation profiling
- Improved peptide enrichment and detection sensitivity
- Community annotation practices for non-canonical ORFs

---

## 9. Conclusion
[Summarize evidence-driven identification and the need for standardized pipelines.]

---

## References
[To be compiled, 80–120 references, grouped by topic:
- computational prediction
- Ribo-seq / translation evidence
- proteomics / peptidomics
- functional validation case studies
- resources & databases]
