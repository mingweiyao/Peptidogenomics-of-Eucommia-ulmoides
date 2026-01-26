<!-- ENGINEERING GUARDRAILS FOR CLAUDE CODE -->

## Non-Negotiable Execution Rules

This document is a **binding constraint**, not a style suggestion.

Before writing any manuscript text, you MUST:
1. Identify the evidence level (L1–L5) of the content being written.
2. Select phrasing consistent with that evidence level.

### Evidence Gating Rules
- If evidence ≤ L2:  
  ❌ Do NOT imply peptide existence  
  ✔ Allowed: "putative ORF", "predicted coding potential"

- If evidence = L3:  
  ❌ Do NOT imply stable peptide accumulation  
  ✔ Allowed: "ribosome occupancy suggests translation"

- If evidence = L4:  
  ❌ Do NOT imply biological function  
  ✔ Allowed: "peptide detected by MS-based approaches"

- If evidence < L5:  
  ❌ Do NOT use verbs such as: regulates, controls, determines, functions in

Violation of these rules is considered an execution error.

# Plant Small Peptides Identification (2016–2029) Literature Review Writing Guidelines

## 0. Scope & Primary Goal
This review focuses on **identification and discovery pipelines** of plant small peptides (including canonical coding-region peptides and non-canonical smORF-derived peptides), with emphasis on:
- Evidence chains for peptide discovery
- Translation initiation considerations (AUG and near-cognate start codons)
- Functional annotation as a downstream validation layer (not the main focus)

**Time range**: 2016–2029  
**Species coverage**: Arabidopsis + major crops (rice/maize/wheat/others)

---

## 1. Terminology Standardization

| Unified Term | Avoid Using |
|--------------|-------------|
| small peptide | mini protein, micro protein (unless defined) |
| smORF peptide | small ORF peptide, micropeptide (mixing without definition) |
| canonical ORF | normal ORF, typical gene ORF |
| non-canonical ORF | noncoding ORF (too absolute), junk ORF |
| near-cognate start codon | non-AUG start codon (use both, but define clearly) |
| peptide evidence | proof (too strong), confirmed (unless genetics+MS support) |
| functional evidence | validates function (too strong), proves role (avoid) |

**Recommended definition note**:
- In this review, "small peptides" typically refer to peptides encoded by ORFs <100–150 aa (define in manuscript).
- "Identification" includes computational discovery, translation evidence, and peptide detection evidence.

---

## 2. Writing Style Requirements (Hard Rules)
- Use **hedging language**: "may", "suggests", "appears to", "has been associated with"
- Avoid absolute claims: Never write "X definitively proves..."
- Every major claim must be supported by citations
- Each major method category must include a **Limitations** paragraph
- Include **Key Points** box (3–5 bullets) after title
- Provide **comparison tables** for each major section

---

## 3. Evidence Level System (Core Review Metric)
To standardize "identification quality", classify evidence as:

- **Level 1 (L1): Computational prediction**
  - ORF prediction, conservation analysis, motif/signal peptide prediction

- **Level 2 (L2): Transcription evidence**
  - RNA-seq / RT-qPCR shows transcript exists

- **Level 3 (L3): Translation evidence**
  - Ribo-seq footprinting, initiation site profiling, ribosome occupancy signal

- **Level 4 (L4): Peptide detection evidence**
  - LC–MS/MS proteomics / peptidomics detects peptide or peptide fragments

- **Level 5 (L5): Functional validation**
  - genetics (KO/CRISPR/OE), synthetic peptide treatment, receptor binding or rescue

When summarizing each peptide candidate or class, always indicate the highest evidence level.

---

## 4. Literature Sources & Search Plan

### 4.1 PubMed / PubMed Central (biological discovery core)
Search keywords (example):
- "plant small peptide identification"
- "plant smORF peptide ribosome profiling"
- "Arabidopsis non-AUG translation smORF"
- "lncRNA encoded peptide plant"
- "peptidomics plant secreted peptide"

Filters:
- Year: 2016–2029
- Article types: Review + Research articles (keep both)
- Focus: Methods / discovery / validation pipelines

### 4.2 Google Scholar (high recall)
Use citation chasing:
- Start from key methods papers, follow "Cited by"
- Use combinations:
  - plant AND smORF AND micropeptide AND ribo-seq
  - plant peptidomics AND secreted peptide AND mass spectrometry
  - non-AUG translation plant ORF

### 4.3 Plant-focused journals (for targeted search)
- The Plant Cell
- Plant Physiology
- New Phytologist
- Nature Plants
- PNAS / Nature Communications (cross-field high-impact)

### 4.4 Databases & Omics Resources (make a dedicated table)
Candidate resource types to include:
- Plant genome annotation portals
- Proteomics repositories
- Ribo-seq repositories / browsers
- small peptide specialized databases (if available)

---

## 5. Major Method Categories to Cover

| Category | Subcategories | Notes | Status |
|----------|---------------|------|--------|
| Computational discovery | ORF prediction / conservation / motif & signal peptide | high recall, low certainty | [ ] |
| Translation-based evidence | Ribo-seq / initiation profiling | central for smORF | [ ] |
| Proteomics detection | LC–MS/MS / peptidomics / targeted MS | direct peptide evidence | [ ] |
| Functional validation | genetics / synthetic peptide / receptor assays | gold standard | [ ] |
| Integration & benchmarking | multi-omics scoring / reproducibility | standards needed | [ ] |

---

## 6. Required Tables & Figures

### Required Tables
- Table 1: Definition & taxonomy of plant small peptides (canonical vs non-canonical)
- Table 2: Method comparison table for identification pipelines (core)
- Table 3: Evidence level scoring rubric and examples
- Table 4: Representative peptides/candidates and evidence levels (case studies)
- Table 5: Databases & resources for peptide discovery

### Figure Placeholders
- Figure 1: Overview pipeline from genome/transcriptome → translation evidence → proteomics → functional validation
- Figure 2: Evidence pyramid (L1–L5) for discovery confidence
- Figure 3: Timeline of major breakthroughs (2016–2029)
- Figure 4: Example workflow for non-canonical smORF discovery

---

## 7. Quality Checklist

### Structure
- [ ] Key Points section (3-5 bullets)
- [ ] Clear scope definition
- [ ] Table per major section
- [ ] Figure placeholders with captions

### Content
- [ ] Coding-region identification covered
- [ ] Non-coding / non-canonical ORF identification covered
- [ ] Start codon diversity discussed (AUG and near-cognate)
- [ ] Evidence-level framework applied consistently
- [ ] Limitations written for each method category

### Language
- [ ] Hedging language used appropriately
- [ ] Terminology consistent
- [ ] Avoid overstating causality

### References
- [ ] 80–120 references typical
- [ ] Include both Arabidopsis and crop studies
- [ ] Include 2016–2029 coverage
