# Research Skills

A collection of Claude Code skills for academic research workflows.

## Skills

| Skill | Description | Trigger |
|-------|-------------|---------|
| [medical-imaging-review](./medical-imaging-review/) | Write comprehensive literature reviews for medical imaging AI | `/medical-imaging-review`, "review paper", "survey", "综述" |
| [paper-slide-deck](./paper-slide-deck/) | Generate professional slides from academic papers with auto figure extraction | `/paper-slide-deck paper.pdf` |

## Installation

Copy the desired skill folder to your Claude Code skills directory:

```bash
# For medical-imaging-review
cp -r medical-imaging-review ~/.claude/skills/

# For paper-slide-deck
cp -r paper-slide-deck ~/.claude/skills/
```

Or copy to project-local skills:

```bash
cp -r medical-imaging-review .agents/skills/
cp -r paper-slide-deck .agents/skills/
```

---

## Medical Imaging Review Skill

A systematic workflow for writing survey papers, systematic reviews, and literature analyses on medical imaging AI topics.

### Features

- **Structured 7-phase workflow** for literature review writing
- **Domain-specific templates** covering multiple medical imaging domains
- **Standardized writing style** with hedging language and citation patterns
- **Quality checklists** ensuring completeness
- **Zotero integration** for reference management

### Supported Domains

- Coronary Artery Analysis (CCTA)
- Lung Imaging (CT/X-ray)
- Brain Imaging (MRI/CT)
- Cardiac Imaging (MRI/CT/Echo)
- Pathology (Whole Slide Images)
- Retinal Imaging (Fundus/OCT)

### Files

| File | Description |
|------|-------------|
| `SKILL.md` | Main skill definition and quick reference |
| `WORKFLOW.md` | Detailed 7-phase workflow guide |
| `TEMPLATES.md` | Project file templates |
| `DOMAINS.md` | Domain-specific method categories and datasets |

---

## Paper Slide Deck Skill

Transform academic papers into professional slide decks with automatic figure extraction and AI-generated visuals.

### Features

- **Auto figure detection** from PDF papers
- **Smart figure-to-slide mapping** based on caption analysis
- **17 visual styles** (academic-paper, sketch-notes, minimal, etc.)
- **Gemini API integration** for AI slide generation
- **PPTX/PDF export** with merge scripts

### Workflow

1. Analyze paper and detect figures/tables
2. Generate outline with auto IMAGE_SOURCE mapping
3. Extract figures from PDF (or AI-generate)
4. Apply academic templates
5. Merge to PPTX/PDF

### Files

| Path | Description |
|------|-------------|
| `SKILL.md` | Main skill definition and workflow |
| `references/` | Analysis framework, templates, style definitions |
| `scripts/` | Python/TypeScript automation scripts |

### Scripts

| Script | Purpose |
|--------|---------|
| `generate-slides.py` | Gemini API image generation |
| `detect-figures.ts` | PDF figure/table detection |
| `extract-figure.ts` | PDF page extraction |
| `apply-template.ts` | Academic figure container template |
| `merge-to-pptx.ts` | PPTX generation |
| `merge-to-pdf.ts` | PDF generation |

---

## License

MIT License
