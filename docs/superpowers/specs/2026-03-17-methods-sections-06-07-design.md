# Design: Methods Documentation for Sections 06 & 07 + Cleanup

**Date**: 2026-03-17
**Scope**: Write clean manuscript-ready supplementary methods for spatial biopsies (06) and organoid single-cell (07), clean existing methods docs (01-05), and consolidate internal review notes into a separate document.

## Deliverables

### 1. New Methods Files

#### `docs/methods/06_spatial_biopsies_methods.md`

Clean supplementary methods covering:

**Data Source & Deconvolution**
- 6 10x Visium biopsies from patient HCC22-088 (Carleton et al.)
- Cell-type deconvolution performed using CITEgeist (Lee Lab) into 9 cell types: Cancer Cells, Immunosuppressive Macrophages (CD163+), Inflammatory Macrophages (CD14+HLA-DR+), T Cell-Interacting Macrophages (CD11c+), CD4 T Cells, CD8 T Cells, B Cells, Endothelial Cells, Fibroblasts
- Sample dimensions: 316-2,251 spots per sample, 2,141-12,251 genes

**Spatial Cytokine Visualization**
- Queried cytokines: IL4, IL10, IL13, TGFB1, TGFB2, TGFB3
- Cell-type-resolved spatial plots using scanpy (viridis colormap)
- Layer-specific expression extracted from CITEgeist pass1 layers per macrophage subtype

**Macrophage Pathway Enrichment**
- Top 200 expressed genes per macrophage subtype per sample
- Enrichment: gseapy enrichr (MSigDB_Hallmark_2020, KEGG_2021_Human, Reactome_2022), FDR<0.05
- GSEA prerank: gene scores = mean expression, 1,000 permutations, seed=12345, min pathway size=5, max=1,000
- IL/TGF pathway filtering via regex
- Summary dotplot: CD163+ macrophage pathway significance (-log10 adjusted p-value) across 6 samples

**Software**
- Python (scanpy, gseapy, seaborn, matplotlib)
- Conda environment: erp_brca_aging
- Seed: 12345

---

#### `docs/methods/07_organoid_scrnaseq_methods.md`

Clean supplementary methods covering:

**Experimental Design**
- PDO-296: patient-derived organoid from study cohort, ER+ breast cancer
- 7 conditions: Vehicle, E1, E2, E1+fulvestrant, E2+fulvestrant, E1+HSD17B7i, E2+HSD17B7i
- 2 pools on 10x Chromium Flex platform

**Cell Ranger & Reference**
- Custom GRCh38-2024-A reference with `filter-probes=false` in Cell Ranger multi config
- Rationale: HSD17B7 is marked as excluded in the default 10x Genomics v1.1.0 probe set; since HSD17B7 is the drug target under investigation, probe filtering was disabled to retain expression measurements for this gene

**Quality Control**
- Filters: min 500 genes, min 1,000 counts, <15% mitochondrial
- Doublet detection: Scrublet per sample (expected_doublet_rate=0.06, n_prin_comps=30, min_gene_variability_pctl=85), samples with <50 cells skipped
- Seed: 42 (preserved for consistency with original analysis)

**Preprocessing & Dimensionality Reduction**
- Gene filtering: expressed in >=3 cells
- Normalization: CPM (target_sum=10,000), log1p
- HVG selection: 3,000 genes (seurat_v3 method on counts layer)
- Scaling: unit variance, max_value=10
- PCA: 50 components on HVGs
- Cell cycle scoring: Tirosh et al. 2016 (43 S-phase genes, 47 G2M genes)
- Neighbors: k=15, n_pcs=30
- UMAP: default parameters
- Clustering: Leiden, resolution=0.5

**Pseudobulk Differential Expression**
- Aggregation: sum of raw counts per sample
- Gene filter: >=10 total counts across samples
- Method: PyDESeq2 (Wald test, dispersion estimation)
- Significance: padj<0.05, |log2FC|>0.5
- Comparisons: E1 vs E2, E1 vs E1+HSD17B7i, E1+HSD17B7i vs E2+HSD17B7i

**Pathway & Transcription Factor Analysis**
- GSEA prerank: MSigDB_Hallmark_2020, rank=log2FC x -log10(pval), 1,000 permutations, seed=42, min_size=15, max_size=500, significance FDR<0.25
- Enrichr over-representation analysis (script 09): cell-level Wilcoxon rank-sum DE (E1 vs E1+HSD17B7i), separate analysis for up/down-regulated genes; 5 libraries (MSigDB_Hallmark_2020, KEGG_2021_Human, Reactome_2022, GO_Biological_Process_2023, WikiPathway_2023_Human), padj<0.05
- PROGENy: top 300 footprints, multivariate linear model (decoupler 2.x), organism=human
- DoRothEA: levels A/B/C, multivariate linear model (decoupler 2.x), organism=human
- Target TFs: ESR1, E2F1, E2F4, MYC, TP53
- Note: pseudobulk DE (PyDESeq2, script 04) and cell-level Wilcoxon DE (script 09) serve different purposes — pseudobulk for the 3 primary comparisons, cell-level for pathway enrichment input

**Single-Cell Pathway Scoring**
- 8 custom gene sets: Estrogen Response Early (23 genes), Estrogen Response Late (20 genes), E2F Targets (24 genes), G2M Checkpoint (20 genes), MYC Targets (22 genes), Proliferation (14 genes), Steroid Biosynthesis (12 genes), ER Targets Direct (15 genes)
- 4 ER subprograms: Proliferative (11), Transcriptional (12), Metabolic (10), Signaling Crosstalk (10)
- Scoring: scanpy sc.tl.score_genes()
- Statistical tests: Mann-Whitney U (medians), Kolmogorov-Smirnov (distributions), Cohen's d (effect size), BH-FDR correction
- Comparisons: E1 vs E2, E1 vs E1+HSD17B7i, E1+HSD17B7i vs E2+HSD17B7i

**HSD17B7 Inhibitor Mechanism Analysis**
- E1/E2-dominant gene signature classification: Wilcoxon rank-sum DE (E1 vs Vehicle, E2 vs Vehicle), significance padj<0.05, |log2FC|>0.25, classify as E1-only, E2-only, or shared responsive
- Signature projection across all conditions to test cycling hypothesis (H1: E2+HSD17B7i gains E1-like character)
- Steroidogenic gene panel: HSD17B family (12 members), CYP19A1, STS, SULT1E1, HSD3B1/2, estrogen receptors (ESR1, ESR2, GPER1)
- ER subprogram dissection: score each subprogram per condition

**Heterogeneity Analysis**
- UMAP KDE density contours per condition
- Neighborhood mixing: per-cell entropy of neighbor treatment composition
- Within-condition dispersion: silhouette scores
- Cluster-treatment association: Fisher's exact test on contingency table
- Biological axis identification: PCA on pathway score columns
- Estrogen response continuum: composite score (mean of ER Early + ER Late), KDE density overlays per condition, sliding-window composition analysis, pairwise KDE overlap matrix, UMAP colored by continuum score

**Proliferation & Quiescence**
- Cycling classification: 2-component GMM on proliferation score (BIC comparison for bimodality), threshold = mean of component means, fallback to 75th percentile
- Quiescence scoring: 8-gene set (CDKN1A, CDKN1B, BTG1, BTG2, TOB1, GAS1, CDKN2A, CCNG2) via sc.tl.score_genes()
- Cycling sub-classification: S vs G2M based on score comparison
- Statistics: Mann-Whitney U (score distributions), Fisher's exact (cycling fractions), Spearman correlation (proliferation-ER coupling)

**Treatment Label Convention**
- h5ad files store treatment values as E1+ICI / E2+ICI (historical); all figures display E1+fulv / E2+fulv (fulvestrant, to avoid confusion with immune checkpoint inhibitors)

**Software**
- Python (scanpy, PyDESeq2, decoupler 2.x, gseapy, scikit-learn GMM, scipy)
- Conda environment: erp_brca_aging
- Seed: 42 throughout (exception to project default of 12345, preserved for reproducibility with original analysis)
- Figures: dual PNG (150 DPI) + PDF (vector)

---

### 2. Clean Existing Methods Files (01-05)

For each of `01_human_bulk_rnaseq_methods.md` through `05_rat_bulk_rnaseq_methods.md`:

**Remove:**
- All `[WARNING]`, `[CRITICAL]`, `[RESOLVED]` tags and associated discussion (formatted as blockquotes `> **[WARNING]** ...` — remove entire blockquote blocks, not just tag text)
- Inline issue commentary (e.g., "doublet rate not calibrated", "applymap deprecated")
- Software bug workarounds and compatibility notes
- Planned-but-unimplemented features (MICA stubs, MCPcounter)

**Keep:**
- All parameter values, statistical methods, thresholds
- Tool versions and software lists
- Analytical workflow descriptions
- Reproducibility notes (seeds, determinism sources)

### 3. Clean `00b_supplementary_methods.md`

Convert to a clean table of contents pointing to sections 01-07 (add entries for 06 and 07 with correct filenames). Remove the issue tracker entirely.

### 4. Consolidated Issues Document

**`docs/code_review/code_review_issues.md`**

Organized by section (01-07), preserving:
- All C1-C4 critical issues with resolution status
- All W1-W19 warnings
- New issues from sections 06-07 continuing existing numbering (W20+ for warnings, C5+ for critical), e.g., W20: gseapy `organism='Human'` parameter in 06 (deprecated/ignored in recent versions)
- Each issue retains: severity, description, resolution status, and any fix applied

### 5. Main Text Additions

Add concise paragraphs to `00_main_text_methods.md` for:

**Spatial Biopsy Analysis (Section 06):**
Dense paragraph covering: HCC22-088 biopsies, CITEgeist deconvolution, 3 macrophage subtypes, cytokine spatial mapping, pathway enrichment (Hallmark/KEGG/Reactome + GSEA prerank), IL/TGF focus, seed 12345.

**Organoid Single-Cell RNA-seq (Section 07):**
Dense paragraph covering: PDO-296 from study cohort, 7 conditions (E1/E2 x Vehicle/fulvestrant/HSD17B7i), 10x Flex with custom reference (filter-probes=false for HSD17B7), Scrublet QC, CPM+log1p, Leiden clustering, PyDESeq2 pseudobulk DE, GSEA/PROGENy/DoRothEA pathway analysis, custom gene set scoring (Mann-Whitney U, Cohen's d, BH-FDR), GMM-based cycling classification, inhibitor mechanism analysis, seed 42.

## File Structure After Implementation

```
docs/
├── methods/
│   ├── 00_main_text_methods.md          (updated: +2 paragraphs)
│   ├── 00b_supplementary_methods.md     (cleaned: TOC only)
│   ├── 01_human_bulk_rnaseq_methods.md  (cleaned: no warnings)
│   ├── 02_rat_snrnaseq_methods.md       (cleaned: no warnings)
│   ├── 03_rat_wes_methods.md            (cleaned: no warnings)
│   ├── 04_human_scrnaseq_methods.md     (cleaned: no warnings)
│   ├── 05_rat_bulk_rnaseq_methods.md    (cleaned: no warnings)
│   ├── 06_spatial_biopsies_methods.md   (NEW)
│   └── 07_organoid_scrnaseq_methods.md  (NEW)
├── code_review/
│   ├── 2026-02-17-integration-issues-log.md  (existing)
│   └── code_review_issues.md                  (NEW: consolidated)
└── superpowers/
    └── specs/
        └── 2026-03-17-methods-sections-06-07-design.md  (this doc)
```
