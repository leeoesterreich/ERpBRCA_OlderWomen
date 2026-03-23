# Computational Methods: Human Bulk RNA-seq Analysis (`01_human_bulk_rnaseq`)

## 1. Data Acquisition

Human ER+ breast cancer bulk RNA-seq data were obtained as a pre-compiled feature count matrix (`HumanERpAge_39404g168s_FeatureCount.txt`) containing 39,404 genes across 168 samples, accompanied by a sample annotation file (`HumanERpAge_BulkRNAseq_SampleInformation.txt`). A parallel log2-transformed TPM matrix (`HumanERpAge_39404g168s_TPMlog2.txt`) was provided as a fallback for analyses not requiring raw counts. Sample metadata included chronological age, tissue type (Tumor or Tumor-Adjacent), and age range classification. Gene annotations were obtained from NCBI (`Homo_sapiens.gene_info.txt`) for protein-coding gene filtering.

For the MICA concordance analysis (subdirectory `mica/`), external validation cohorts were planned from TCGA, METABRIC, and SCAN-B datasets (`ilc_clean_data.RData`, `METABRIC_Clinical_Info.csv`, `TCGA_Key.csv`).

### Age Group Definitions

Samples were classified into age groups as follows:

- **Primary analysis (main pipeline):** Age groups were derived from the sample annotation column `AgeRange` (categorical values: Young, Middle, Elderly). Numeric age boundaries are Young: 35–44 years, Middle: 54–69 years, Elderly: 70–85 years; the correlation analysis (script `04_correlations.R`) restricted to Young and Elderly tumor samples only, excluding Middle-aged samples.
- **MICA sub-analysis:** Young: 35-45 years; Middle-Aged: 55-69 years; Elderly: >= 70 years. A postmenopausal subdivision was also defined: Early: 55-59, Middle: 60-69, Elderly: >= 70.

## 2. Preprocessing

Preprocessing was performed in `01_preprocess.R` with two operational modes depending on input data availability.

### Count Mode (Primary)

When raw feature counts were available, the pipeline performed:

1. **Gene annotation filtering:** Genes were filtered to retain only protein-coding genes from the NCBI gene annotation file, excluding genes matching the pattern `^LOC\d+` (uncharacterized loci). This was applied using `dplyr::filter(type_of_gene == "protein-coding", !grepl("^LOC\\d+", Symbol))`.

2. **Sample name harmonization:** Sample identifiers were stripped of the `_LEE*` suffix using `gsub("_LEE(.*)", "", colnames(count_data))` to match annotation file naming.

3. **Expression normalization:** Pre-computed log2(TPM) values (`HumanERpAge_39404g168s_TPMlog2.txt`) were used as the expression matrix. The matrix was deduplicated by gene symbol (the row with the highest mean expression across samples was retained).

Note: The pipeline supports DESeq2 VST normalization from raw counts as an alternative path, but the current analysis uses the log2(TPM) matrix as the primary input.

### Quality Control

QC visualizations were generated as a multi-page PDF (`qc_plots.pdf`):
- **PCA plot:** Principal component analysis with `prcomp(t(expr_matrix), center = TRUE, scale. = TRUE)` (TPM mode) or `plotPCA()` from DESeq2 (count mode), colored by age range and shaped by tissue type.
- **Sample similarity:** Pearson correlation heatmap (TPM mode) or Euclidean distance heatmap (count mode) with hierarchical clustering, annotated by age range and tissue type, generated using `pheatmap::pheatmap()`.
- **Library size distribution:** Bar plot of per-sample library sizes in millions (count mode only).

### Duplicate Gene Handling

When multiple rows shared the same gene symbol, the row with the highest mean expression across samples was retained.

### Random Seed

`set.seed(12345)` was set at the start of all scripts for reproducibility.

## 3. Core Analyses

### 3.1 Gene Set Variation Analysis (GSVA)

GSVA was performed in `02_run_gsva.R` to quantify per-sample estrogen pathway activity.

**Gene sets assembled:**

| Gene Set | Source | Category |
|---|---|---|
| HALLMARK_ESTROGEN_RESPONSE_EARLY | MSigDB Hallmark (H) | `msigdbr(species = "Homo sapiens", category = "H")` |
| HALLMARK_ESTROGEN_RESPONSE_LATE | MSigDB Hallmark (H) | `msigdbr(species = "Homo sapiens", category = "H")` |
| E1UpRegGene | Custom (407 genes; Liguori et al., Cell Metabolism 2020, Table S1: genes uniquely upregulated by E1 not E2, FC>2, q<0.05) | `SuppleTable1_UpregulatedByE1_NotE2_407g.txt` |
| REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION | MSigDB C2:REACTOME | `msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")` |
| WP_ESTROGEN_SIGNALING_PATHWAY | MSigDB C2:WIKIPATHWAYS | `msigdbr(species = "Homo sapiens", category = "C2", subcategory = "WIKIPATHWAYS")` |
| GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY | MSigDB C5:BP | `msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")` |
| GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS | MSigDB C5:BP | `msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")` |

For the E1-upregulated gene set, gene symbols containing hyphens were truncated at the first hyphen (`gsub("-.*", "", e1_genes$Gene)`), and duplicates were removed.

**GSVA parameters:**
- Function: `gsva(gsvaParam(...))` (GSVA 2.x API)
- `kcdf = "Gaussian"` (appropriate for continuous normalized data)
- `maxDiff = TRUE` (Gaussian-like enrichment statistic)
- Input: log2(TPM) expression matrix

**Visualization:** Row-scaled heatmap of GSVA enrichment scores across all samples, annotated by age range (Young = green `#4DAF4A`, Middle = blue `#377EB8`, Elderly = red `#E41A1C`) and tissue type (Tumor = purple `#984EA3`, TumorAdj = orange `#FF7F00`). Generated with `pheatmap()` with `scale = "row"` and `show_colnames = FALSE`.

### 3.2 PROGENy Pathway Activity

PROGENy pathway activity inference was performed in `03_run_progeny.R`.

**Input data:** log2(TPM) expression matrix (`vst_normalized_matrix.rds`, named for historical compatibility) produced by `01_preprocess.R`, consistent with the GSVA input. The matrix was filtered to protein-coding genes (same filter as preprocessing), with sample names matched to annotation via inner join.

**PROGENy parameters:**
- Function: `progeny::progeny()`
- `scale = FALSE`
- `organism = "Human"`
- `top = 100` (number of footprint genes per pathway)
- `perm = 1000` (permutations for significance assessment)

### 3.3 Correlation Analysis

Spearman rank correlations between gene expression and pathway activity scores were computed in `04_correlations.R`.

**Sample selection:** Only tumor samples from the Young and Elderly age groups were included. Middle-aged samples and tumor-adjacent samples were excluded (`filter(Group == "Tumor", AgeRange %in% c("Young", "Elderly"))`).

**Genes of interest:** PAK4, HSD17B7, GREB1, PGR, ESR1, TFF1, CYP19A1, HSD17B2, SAA1, and chronological Age (9 genes + 1 clinical variable).

**Pathway variables:** PROGENy Estrogen pathway activity score, plus all GSVA enrichment scores (7 estrogen-related pathway scores, see Section 3.1).

**Correlation method:**
- `cor.test(x, y, method = "spearman", exact = FALSE)` for each gene-pathway pair
- Minimum sample size threshold: 5 non-NA paired observations required per test
- Gene expression values were extracted from the raw log2(TPM) file (not the VST matrix)

**Multiple testing correction:**
- Benjamini-Hochberg FDR applied globally across all correlation tests using `p.adjust(method = "BH")`
- Significance thresholds: nominal p < 0.05 and FDR q < 0.05

**Output:** Both original (uncorrected) and FDR-corrected results were saved to separate directories (`results/original/` and `results/corrected/`). A combined RDS file was also saved for visualization.

The correlation analysis (`04_correlations.R`) filters the TPM matrix to protein-coding genes using the same NCBI gene annotation and filtering logic as `03_run_progeny.R`, ensuring consistent gene content across all downstream analyses.

## 4. Downstream Analyses

### 4.1 MICA Concordance Analysis

The MICA (Mutual Information-based Concordance Analysis) sub-pipeline validated estrogen pathway findings across independent cohorts (TCGA, METABRIC, SCAN-B). GSVA was used to compute pathway activity scores across cohorts with two age stratification schemes, and concordance was assessed using 500 permutations (`N_PERM = 500`) at a significance threshold of `P_THRESHOLD = 0.05` with `N_PARALLEL = 10` parallel threads. CPM normalization was applied via `edgeR::cpm()` for count data.

### 4.2 HSD17B Gene Family Analysis (Planned)

Figure 4C was planned to display a heatmap of median expression for 14 genes (ESR1, GREB1, PGR, SAA1, RAB19, KRT37, TRPM8, CYP19A1, HSD17B1, HSD17B7, HSD17B12, HSD17B2, HSD17B10, HSD17B14) across age groups, aggregated by median and visualized using `ComplexHeatmap::Heatmap()`.

## 5. Visualization

### Implemented Visualizations

| Output File | Type | Script | Description |
|---|---|---|---|
| `qc_plots.pdf` | Multi-page PDF | `01_preprocess.R` | PCA, sample correlation/distance heatmap, library size bars |
| `gsva_heatmap.pdf` | PDF | `02_run_gsva.R` | Row-scaled heatmap of GSVA scores (14 x 8 inches) |
| `correlation_bubbleplot_original.pdf` | PDF + PNG | `05_visualize.R` | Bubble plot: size = -log10(p), color = Spearman rho (12 x 7 inches) |
| `correlation_bubbleplot_fdr.pdf` | PDF + PNG | `05_visualize.R` | Bubble plot with FDR-corrected p-values (12 x 7 inches) |
| `correlation_bubbleplot_comparison.pdf` | PDF + PNG | `05_visualize.R` | Side-by-side original vs FDR using patchwork (20 x 8 inches) |

**Bubble plot aesthetics:**
- Color scale: blue-dodgerblue-yellow-orange-red gradient via `colorRampPalette()`, mapped to Spearman rho with limits [-1, 1]
- Size: `-log10(p-value)` or `-log10(FDR)`, range 1-10
- Gene order (y-axis): TFF1, SAA1, PGR, PAK4, HSD17B7, HSD17B2, GREB1, ESR1, CYP19A1, Age
- RAB19 was explicitly filtered out of visualization (`filter(!grepl("RAB19", GeneSymb))`)
- Pathway names were abbreviated: `HALLMARK_` to `HM_`, `REACTOME_` to `RC_`, `GOBP_` to `GO_`

### Planned Visualizations (Not Implemented)

| Output File | Type | Script | Description |
|---|---|---|---|
| `fig4c_hsd17b_heatmap.pdf` | PDF + PNG | `mica/04_fig4_estrogen.R` | HSD17B gene heatmap across cohorts |
| `fig4h_estrogen_boxplots.pdf` | PDF + PNG | `mica/04_fig4_estrogen.R` | EstroGene pathway boxplots by age |
| `fig5c_inflammatory_heatmap.pdf` | PDF + PNG | `mica/05_fig5_immune.R` | Inflammatory pathway heatmap |
| `fig5e_immune_boxplots.pdf` | PDF + PNG | `mica/05_fig5_immune.R` | Immune cell abundance by age |

## 6. Software Versions

From `environment.yml` (conda environment `erp_brca_aging`):

| Package | Version | Purpose |
|---|---|---|
| R | 4.4.1 | Base language |
| DESeq2 | (Bioconductor, version not pinned) | Count normalization and VST |
| GSVA | (Bioconductor, version not pinned) | Gene set variation analysis |
| PROGENy | (Bioconductor, version not pinned) | Pathway activity inference |
| msigdbr | (Bioconductor, version not pinned) | MSigDB gene set retrieval |
| tidyverse | (not pinned) | Data manipulation |
| data.table | (not pinned) | Fast I/O |
| pheatmap | (not pinned) | Heatmaps |
| ggplot2 | (not pinned) | Plotting |
| patchwork | (not pinned) | Multi-panel figures |
| edgeR | (not pinned) | CPM normalization (MICA sub-pipeline) |
| ComplexHeatmap | (not listed in environment.yml) | MICA figure generation |
| circlize | (not listed in environment.yml) | Color functions for ComplexHeatmap |
| MICA | (not listed in environment.yml) | Concordance analysis |

### SLURM Execution Parameters

| Parameter | Main Pipeline (`run_analysis.sbatch`) | MICA Pipeline (`mica/run_all.sbatch`) |
|---|---|---|
| Partition | htc | htc |
| Time limit | 4 hours | 4 hours |
| Memory | 32 GB | 32 GB |
| CPUs | 4 | 10 |
| Conda env | `erp_bulk_rnaseq` | `/ix1/.../envs/erp_brca_aging` |

## 7. Reproducibility Notes

- **Random seeds:** `set.seed(12345)` was set at the start of every R script in the pipeline.
- **Execution order:** Scripts were executed sequentially (`01` through `05`) via shell scripts, with each script depending on outputs from prior steps. The sbatch script enforced sequential execution with `set -eo pipefail`.
- **Completion markers:** Upon successful completion, the sbatch script wrote a marker file to `.pipeline_markers/01_human_bulk_rnaseq.complete`.
- **Intermediate checkpoints:** All intermediate results (DESeq2 object, VST matrix, GSVA scores, PROGENy activity, correlation results) were saved as RDS files. Key results were additionally exported as CSV for cross-platform accessibility.
- **Hardcoded paths:** The sbatch script contains an absolute path to `SCRIPT_DIR`. All R scripts derive paths relative to the script location using `commandArgs()`.

---

## Main Text Summary

Human ER+ breast cancer bulk RNA-seq data (39,404 genes, 168 samples) were preprocessed by filtering to protein-coding genes using NCBI gene annotations and normalizing with DESeq2 variance stabilizing transformation (design: ~ AgeRange + Group, `blind = FALSE`). Estrogen pathway activity was quantified using GSVA (v2.x, `kcdf = "Gaussian"`, `maxDiff = TRUE`) across seven curated gene sets from MSigDB Hallmark, Reactome, WikiPathways, GO Biological Process, and a custom E1-responsive gene set (407 genes). PROGENy pathway activity was inferred using the top 100 footprint genes per pathway with 1,000 permutations. Spearman rank correlations between nine estrogen-related genes (PAK4, HSD17B7, GREB1, PGR, ESR1, TFF1, CYP19A1, HSD17B2, SAA1) plus chronological age and all pathway activity scores were computed for Young and Elderly tumor samples, with Benjamini-Hochberg false discovery rate correction applied across all tests. MICA concordance analysis was performed across TCGA, METABRIC, and SCAN-B cohorts using GSVA with two age stratification schemes and 500 permutations. All analyses used R 4.4.1 with seed 12345 for reproducibility.
