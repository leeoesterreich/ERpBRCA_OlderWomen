# Computational Methods: Human Single-Cell RNA-seq Analysis (Section 04)

## 1. Data Acquisition

Single-cell RNA-seq data were obtained from two sources. The primary dataset was the Xu et al. 2024 Primary Breast Tumor Atlas, comprising treatment-naive HR+ patients from 9 studies. Raw count matrices (`matrix.mtx`), gene annotations (`genes.tsv`), cell barcodes (`barcodes.tsv`), and per-cell metadata (`metadata.csv`) were loaded from a local institutional copy. An earlier version of the pipeline used Wu et al. (2021) data from GEO accession GSE176078 (10 ER+ patients, ~31,000 cells); script `00a_download_geo.sh` downloads this archive via `wget` from `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz`. The final pipeline used the Xu et al. atlas exclusively.

From the Xu atlas, cells were filtered to those annotated as `Clinical_Subtype == "HR+"` and `Treatment_Status == "Naive"`, yielding approximately 36 patients and ~115,000 cells. Patient 0319 was subsequently excluded as ER-low/PR+ only (`01_load_subset_data.R`, line 54).

## 2. Preprocessing and Quality Control

A Seurat (v5.x) object was created from the sparse count matrix using `CreateSeuratObject()` (`00b_preprocess_seurat.R`). Gene names were deduplicated via `make.unique()`. Pre-computed `nCount_RNA` and `nFeature_RNA` values from the atlas metadata were discarded in favor of Seurat-recomputed values from raw integer counts (line 124).

Quality control filtering applied the following thresholds:
- Minimum genes per cell (`nFeature_RNA`): 200
- Maximum genes per cell (`nFeature_RNA`): 6,000
- Minimum UMI counts per cell (`nCount_RNA`): 400
- Maximum mitochondrial gene percentage (`percent.mito`): 15%

Mitochondrial content was computed using `PercentageFeatureSet()` with the pattern `^MT-`.

## 3. Normalization and Dimensionality Reduction

Normalization was performed using SCTransform (`SCTransform()`) with the `glmGamPoi` method for computational efficiency on the ~115K-cell dataset (`00b_preprocess_seurat.R`, line 170). Principal component analysis (PCA) was computed on the top 30 components (`npcs = 30`).

Batch correction across the 9 contributing studies was performed using Harmony (`RunHarmony()`, `group.by.vars = "Dataset"`, `assay.use = "SCT"`). UMAP embeddings were computed on the Harmony-corrected space using dimensions 1-30 (`RunUMAP()`, `reduction = "harmony"`, `dims = 1:30`).

For downstream analyses requiring log-normalized data (e.g., violin plots, FindMarkers), `NormalizeData()` and `ScaleData()` were run on the RNA assay. `ScaleData()` was applied to variable features plus a curated set of key markers (CCL2, TGFB1, CD163, MRC1, EPCAM, KRT19, CD68, CD3D, MS4A1, PECAM1, CTLA4, PDCD1) rather than all genes, to avoid the ~50 GB memory requirement of a full 59K x 115K scale.data matrix (line 219-230).

Random seed was set to 12345 in all R scripts via `set.seed(12345)`.

## 4. Cell Type Annotation

Cell type labels were transferred directly from the Xu et al. 2024 atlas metadata field `Cell_Type_Annotation` and mapped to pipeline-standard names via a hardcoded dictionary (`03_cell_type_annotation.R`, lines 43-62). The mapping included:

| Xu Annotation | Pipeline Name |
|---|---|
| Cancer Epithelial Cells | CancerEpithelial |
| CD4+ T Cells | TcellsCD4 |
| CD8+ T Cells | TcellsCD8 |
| Regulatory T Cells | Tregs |
| NK Cells | NKcells |
| B Cells | Bcells |
| Plasma Cells | Plasmablasts |
| Macrophages | Macrophage |
| Monocytes | Monocyte |
| Dendritic Cells | DCs |
| Fibroblasts | CAFs |
| Endothelial Cells | Endothelial |
| Perivascular-like (PVL) Cells | PVL |
| Epithelial Cells | NormalEpithelial |
| Myoepithelial Cells | Myoepithelial |
| MDSCs | MDSCs |
| Mast Cells | MastCells |
| Neutrophils | Neutrophils |

Unmapped cell types were labeled "Unknown". No de novo clustering or marker-based annotation was performed; the atlas-provided labels were used as given.

## 5. Age Group Assignment

Patients were assigned to three age groups based on parsed ages from metadata (`01_load_subset_data.R`):
- **Young**: age <= 50
- **MidAge**: 51-80
- **Elderly**: > 80

Age strings in range format (e.g., "60-65") were converted to midpoints. Cells with unparseable or missing ages were excluded.

## 6. Cell Type Fraction Analysis

Per-patient cell type proportions were computed by dividing cell counts per type by total cells per patient (`04_cell_fractions.R`). Statistical comparison of Young vs Elderly fractions used the Wilcoxon rank-sum test (`wilcox.test()`, `exact = FALSE`) for each cell type independently. P-values were adjusted across all cell types using the Benjamini-Hochberg (BH) FDR method (`p.adjust()`, `method = "BH"`). Significance was determined at FDR < 0.05.

## 7. Gene Expression Violin Plots

Violin plots were generated per cell type across age groups for five gene sets (`05_gene_expression_violin.R`):
1. **Immune checkpoint genes**: CTLA4, PDCD1, PDCD1LG2, CD274, LAG3, HAVCR2 (plotted from log-normalized `data` slot)
2. **M2 macrophage markers**: CCL2, CCL3, CCL4, TNF, TGFB1, CD163, MRC1 (plotted from z-scored `scale.data` slot)
3. **CellPhoneDB interaction genes**: CCR1, CCR3, CCR4, CD74, APP, CXCL10, DPP4, EGFR, AGER, GPR75, CCL5, IL2RG, IL15, S100A11 (log-normalized)
4. **Cytokines/chemokines**: 28 genes including IL1A, IL1B, IL6, CXCL8, IL10, IFNG, CXCL1-3, and others (log-normalized)
5. **Hedgehog pathway genes**: SHH, KLF4, GLI1, KIF7, SMO (log-normalized, grouped by cell type rather than age)

A `safe_vlnplot()` wrapper function was used to handle plot combination using `wrap_plots()`. Cell types with fewer than 10 cells were skipped. M2 markers were additionally plotted per patient (CaseID), with patients ordered by ascending age.

## 8. GSVA Pathway Analysis

Gene Set Variation Analysis (GSVA) was performed in two modes (`06_run_gsva.R`):

**Gene sets:** Estrogen-related pathways were assembled from MSigDB via the `msigdbr` package (species = "Homo sapiens"), spanning four collections:
- Hallmark: ESTROGEN_RESPONSE_EARLY, ESTROGEN_RESPONSE_LATE
- WikiPathways: ESTROGEN_METABOLISM_WP697, WP5276, ESTROGEN_SIGNALING
- Reactome: 6 estrogen-related sets (biosynthesis, dependent gene expression, dependent nuclear events, stimulated signaling PRKCZ, extra-nuclear signaling, RUNX1-regulated ER transcription)
- GO Biological Process: RESPONSE_TO_ESTROGEN, ESTROGEN_RECEPTOR_SIGNALING_PATHWAY, POSITIVE_REGULATION_OF_INTRACELLULAR_ESTROGEN, CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS
- LI EstroGene early/late E2 response up (loaded from local GMT files)

**Pseudo-bulk mode:** Raw counts were aggregated per patient via `rowSums()`, converted to CPM, and log2-transformed (`log2(CPM + 1)`). GSVA was run with `gsvaParam()` using `kcdf = "Gaussian"` and `maxDiff = TRUE`.

**Single-cell mode (original approach):** Raw counts were filtered: genes with `rowSums > 10`, removal of genes containing "." in name (proxy for non-protein-coding), and cells with `colSums > 1000` (targeting ~18,063 genes x ~28,732 cells per original code). GSVA was run with `kcdf = "Poisson"` and `maxDiff = TRUE`. Per-cell GSVA scores were aggregated to per-patient means.

Pathway activity was stratified by HSD17B7 expression status (median split into HSD17B7+ and HSD17B7- groups). Heatmaps were generated using ComplexHeatmap with a diverging blue-white-red color scale, clamped to [-4, 4].

## 9. PROGENy Pathway Activity

PROGENy pathway activity scores were computed on pseudo-bulk profiles (`07_run_progeny.R`). Per-patient pseudo-bulk was created by averaging log-normalized expression across cells per patient (`rowMeans()` of the `data` slot). PROGENy was run with `top = 100` footprint genes per pathway, `perm = 1000` permutations, `scale = FALSE`, and `organism = "Human"`. Heatmaps were annotated by age group and visualized with `pheatmap()` using column-wise z-scoring (`scale = "column"`).

## 10. WCSEA Pathway Analysis

Weighted Concept Signature Enrichment Analysis (WCSEA) was attempted using the `indepthPathway` package (`08_run_wcsea.R`). Differentially expressed genes between Elderly and Young (all cell types pooled) were identified using Seurat's `FindMarkers()` with `min.pct = 0.1` and `logfc.threshold = 0.25`. FDR-adjusted p-values (Seurat's `p_val_adj`) were used at a threshold of 0.05. If the `indepthPathway` package was available, WCSEA was run with `pathway_db = "GO_BP"`, `min_genes = 10`, and `max_genes = 500`, with BH-FDR correction applied to output p-values. If unavailable, significant marker genes were saved for manual analysis.

## 11. CellPhoneDB Input Preparation

Input files for CellPhoneDB were prepared in `09_cellphonedb_prep.R` using three filtering strategies:

**Young/Elderly subsets:** Genes with `rowSums > 20` were retained. Hyphens in cell barcodes were replaced with underscores. Special characters in cell type names were stripped.

**Whole cohort:** A three-step filter matching the original code was applied: (1) `rowSums > 10`, (2) removal of genes containing "." in name, (3) cells with `colSums > 1000`. This targeted ~18,063 genes x ~28,732 cells.

Output files consisted of tab-separated count matrices (gene x cell) and metadata tables (Cell, cell_type) for each group.

## 12. Macrophage Subsetting

Macrophage and monocyte cells were extracted from the full annotated Seurat object by filtering on `CellTypeAnnotSH %in% c("Macrophage", "Monocyte")` (`10_load_macrophage_seurat.R`, line 92). The resulting object was saved as `macrophage_seurat.rds`.

## 13. Macrophage Pseudobulk Differential Expression

Differential gene expression between elderly and young macrophages was performed using a pseudobulk DESeq2 approach to avoid single-cell pseudoreplication (`11_macrophage_deg_analysis.R`). Cells annotated as `CellTypeAnnotSH == "Macrophage"` were subset to Elderly and Young age groups. Patients with fewer than 10 macrophages were excluded (`MIN_CELLS_PER_PSEUDOBULK = 10`).

Raw counts were aggregated per patient using `AggregateExpression()` (Seurat). Low-count genes were filtered, requiring >= 10 counts in >= 3 patients. DESeq2 was run with the design formula `~ AgeGroup` and results extracted with `contrast = c("AgeGroup", "Elderly", "Young")` at `alpha = 0.05`. DEGs were defined at FDR < 0.05. Volcano plots used thresholds of |log2FoldChange| > 0.5 and FDR < 0.05 for coloring.

**Output column names:** `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj` (DESeq2 convention).

## 14. Macrophage Pathway Enrichment

Pathway enrichment on macrophage DEGs was performed using a hypergeometric test (one-sided Fisher's exact test) (`12_macrophage_pathway_enrichment.R`). Gene set collections included MSigDB Hallmark (50 pathways) and BioCarta (`category = "C2"`, `subcategory = "CP:BIOCARTA"`). Enrichment was run separately for upregulated and downregulated gene sets. P-values were adjusted using BH-FDR. Pathways were considered significant at FDR < 0.05 (with a relaxed display threshold of FDR < 0.1 for visualization).

This script uses DESeq2-style column names (`log2FoldChange`, `padj`) from the DEG input file produced by `11_macrophage_deg_analysis.R`.

## 15. ENRICHR Pathway Analysis

ENRICHR pathway analysis was performed via the GSEApy Python package (`12_enrichr_pathway_analysis.py`). Six ENRICHR libraries were queried per gene list:
- MSigDB_Hallmark_2020
- KEGG_2021_Human
- Reactome_Pathways_2024
- GO_Biological_Process_2023
- WikiPathway_2023_Human
- BioCarta_2016

Each library was queried independently to ensure adjusted p-values were computed within each library rather than across all libraries combined. Analyses were run on per-cell-type pseudobulk DEGs from script 14 and macrophage-specific DEGs from script 11. For macrophage-specific DEGs, a relaxed nominal p-value threshold (p < 0.01) was used as a fallback when no genes passed FDR < 0.05.

Heatmaps of signed -log10(FDR) scores were generated across cell types, with significance markers (* FDR < 0.05, ** FDR < 0.01). Row clustering used Ward's method with Euclidean distance. A maximum of 30 pathways were displayed per heatmap, selected by balanced representation of upregulated and downregulated pathways.

## 16. Macrophage Cell-Cell Communication Visualization

Script `13_macrophage_cellphonedb.R` generated a bar chart visualization of communication-related DEGs from the macrophage pseudobulk DE results, filtering for genes matching chemokine (CCL, CXCL), cytokine (IL, TNF, TGF, IFNG), and immune checkpoint (CD80, CD86, PDCD1, CD274, CTLA4, LAG3) patterns. This is a DEG-based visualization, not a CellPhoneDB or CellChat interaction analysis. The script loaded the full annotated Seurat object (`seurat_annotated.rds`) as a fallback but primarily operated on `macrophage_seurat.rds`.

This script uses DESeq2-style column names (`log2FoldChange`, `padj`) from the DEG input file.

## 17. Multi-Cell-Type Pseudobulk Differential Expression

Pseudobulk DEG analysis was performed for each cell type independently (`14_multicelltype_deg.R`). For each cell type, cells from Elderly and Young patients were aggregated per patient using `AggregateExpression()`. Patients required >= 10 cells per cell type (`MIN_CELLS_PER_PSEUDOBULK = 10`) and >= 3 patients per age group (`MIN_PATIENTS_PER_GROUP = 3`). Gene filtering required >= 10 counts in >= 3 patients. DESeq2 was run with `design = ~ AgeGroup` and `contrast = c("AgeGroup", "Elderly", "Young")` at `alpha = 0.05`. Per-cell-type DEG tables were saved with DESeq2-style column names (`gene`, `log2FoldChange`, `padj`).

## 18. Multi-Cell-Type Pathway Enrichment

Per-patient pseudo-bulk GSVA was performed within each cell type, followed by statistical comparison of pathway activity between age groups (`15_multicelltype_pathway.R`). Two modes were available: "curated" (default, using 25 pre-selected Hallmark + BioCarta pathways) and "divergent" (data-driven selection of most divergent pathways).

**Input data:** Full annotated Seurat object (`seurat_annotated.rds`) containing all cell types, enabling proper cross-cell-type pathway analysis.

**Pseudo-bulk construction:** For each cell type present in the object, per-patient expression profiles were created by averaging SCTransform-normalized data (`GetAssayData()`, layer = "data") across cells, requiring >= 10 cells per patient per cell type (`min_patients_per_group = 2` in script 15, compared to `MIN_PATIENTS_PER_GROUP = 3` in the DESeq2-based script 14). If the SCT assay was unavailable, RNA counts were normalized to CPM + log2.

**GSVA parameters:** `gsvaParam()` with `kcdf = "Gaussian"` and `maxDiff = TRUE`. Gene sets included all MSigDB Hallmark pathways, BioCarta (subcategory "CP:BIOCARTA"), and LI_ESTROGENE early/late E2 response gene sets from local GMT files.

**Statistical testing:** For each pathway within each cell type, Elderly vs Young patient GSVA scores were compared using a Welch t-test (`t.test(..., var.equal = FALSE)`). FDR correction was applied using two strategies: (1) across all pathway-cell type combinations (`padj_all`), and (2) restricted to 25 curated pathways only (`padj`, reducing multiple testing burden from ~4000+ to ~312 tests). Significance was determined at FDR < 0.05.

**Curated pathway set (Hallmark):** ESTROGEN_RESPONSE_EARLY, ESTROGEN_RESPONSE_LATE, INFLAMMATORY_RESPONSE, TNFA_SIGNALING_VIA_NFKB, TGF_BETA_SIGNALING, IL6_JAK_STAT3_SIGNALING, IL2_STAT5_SIGNALING, INTERFERON_GAMMA_RESPONSE, INTERFERON_ALPHA_RESPONSE, EPITHELIAL_MESENCHYMAL_TRANSITION, ANGIOGENESIS, HYPOXIA, APOPTOSIS, plus LI_ESTROGENE_EARLY_E2_RESPONSE_UP and LI_ESTROGENE_LATE_E2_RESPONSE_UP. The unified set includes IL2_STAT5, IFN_ALPHA, and EMT; COMPLEMENT is excluded.

**Curated pathway set (BioCarta):** INFLAM_PATHWAY, IL6_PATHWAY, IL2_PATHWAY, NFKB_PATHWAY, TNFR1_PATHWAY, DEATH_PATHWAY, FAS_PATHWAY, CASPASE_PATHWAY, P53_PATHWAY, CELLCYCLE_PATHWAY, G1_PATHWAY, G2_PATHWAY.

**Visualization:** Dual heatmaps (Hallmark and BioCarta) were generated using ComplexHeatmap. GSVA scores were z-scored per pathway across all cell type-age group combinations. Rows were grouped by cell lineage category (Lymphocyte, Myeloid, Epithelial, Stromal). Columns were hierarchically clustered. Significance stars (* FDR < 0.05, ** FDR < 0.01) were overlaid. Color scale: blue-white-orange-red diverging palette (Hallmark: [-3, 3]; BioCarta: [-4, 4]).

## 19. GSVA Polarity Diagnostic

A diagnostic script (`15b_pathway_polarity_diagnostic.R`) was implemented to investigate potential GSVA score polarity reversals between normalization methods. This script compared GSVA scores computed from CPM + log2-normalized pseudo-bulk versus SCTransform-normalized data for five key Hallmark pathways (ESTROGEN_RESPONSE_EARLY, ESTROGEN_RESPONSE_LATE, INFLAMMATORY_RESPONSE, INTERFERON_GAMMA_RESPONSE, TNFA_SIGNALING_VIA_NFKB). Z-score normalization effects on polarity were also assessed.

## 20. CellChat Cell-Cell Communication Analysis

CellChat analysis was performed to infer cell-cell communication between age groups (`16_cellchat_analysis.R`). CellChat objects were created separately for Elderly and Young subsets from the full annotated Seurat object (`seurat_annotated.rds`). The CellChatDB.human database was used, subset to "Secreted Signaling" interactions.

CellChat preprocessing steps included:
1. `subsetData()` - subset to signaling genes
2. `identifyOverExpressedGenes()` - identify overexpressed ligands/receptors
3. `identifyOverExpressedInteractions()` - identify overexpressed L-R pairs
4. `computeCommunProb()` with `type = "triMean"` (truncated mean) and `min.cells = 10`
5. `computeCommunProbPathway()` - aggregate to pathway level
6. `aggregateNet()` - aggregate communication network

Young and Elderly CellChat objects were merged using `mergeCellChat()` for comparative visualization. Bubble plots were generated using `netVisual_bubble()` with `color.heatmap = "Spectral"`.

## 21. CellPhoneDB Statistical Analysis

CellPhoneDB was run via the Python API in two versions:

**CellPhoneDB v5** (`16_run_cellphonedb.py`): Used `cpdb_statistical_analysis_method.call()` with parameters: `counts_data = "gene_name"`, `threshold = 0.1` (minimum fraction of cells expressing a gene), `iterations = 1000` (permutations), `threads = 4`. No random seed was set.

**CellPhoneDB v4** (`16b_run_cellphonedb_v4.py`): Used the same API with `counts_data = "hgnc_symbol"`, `threshold = 0.1`, `threads = 4`, `debug_seed = 42`, `result_precision = 3`.

Both versions were run separately for Elderly and Young age groups using pre-prepared count matrices and metadata files from script 09.

### 21.1 Independent GSEA (Script 09)

Script 09 (`09_pathway_enrichment_analysis.py`) performed an independent GSEA using cell-level Wilcoxon rank-sum z-scores for gene ranking (distinct from the pseudobulk-based log2FC x -log10(p) ranking in script 05), querying three libraries: MSigDB_Hallmark_2020, KEGG_2021_Human, and Reactome_2022.

## 22. CellPhoneDB Dot Plot Visualization

A comparative dot plot of CellPhoneDB results was generated in R (`17_cellphonedb_dotplot.R`). Two modes were available:
- **top50** (default): Despite the mode name, selects the top 40 ligand-receptor pairs (via `head(40)`) by minimum p-value across all macrophage-immune cell interactions
- **curated**: A predefined list of L-R pairs with normalized name matching to handle ordering differences between CellPhoneDB versions

The plot focused on Macrophage-as-sender interactions with six target cell types: B cells, Cycling T cells, NK cells, NKT cells, CD4+ T cells, and CD8+ T cells. Dot size encoded `-log10(p-value)` (capped at 3), and dot color encoded `log2(mean expression + 1)` using an RdBu diverging palette (range: -10 to +5). Only interactions with BH-FDR corrected p < 0.05 were displayed. X-axis labels were placed on top (manuscript convention). Age groups were labeled "Younger" and "Older".

## 23. Software Versions

From `environment.yml` (conda environment `erp_brca_aging`):

| Software | Version |
|---|---|
| R | 4.3.3 |
| Seurat | >= 5.0 |
| Harmony | (via conda) |
| glmGamPoi | (Bioconductor) |
| DESeq2 | (Bioconductor) |
| GSVA | (Bioconductor) |
| PROGENy | (Bioconductor) |
| msigdbr | (Bioconductor) |
| ComplexHeatmap | (via R) |
| pheatmap | (via R) |
| CellChat | (pre-installed, version not pinned in environment.yml) |
| indepthPathway | (optional, may not be installed) |
| Python | 3.10 |
| GSEApy | (via pip, version not pinned) |
| CellPhoneDB | v4/v5 (version not pinned in environment.yml) |
| SingleR | (Bioconductor, loaded but not used in current scripts) |
| speckle | (for compositional analysis, loaded but not used in current scripts) |
| pandas | (via conda) |
| numpy | (via conda) |
| matplotlib | (via conda) |
| seaborn | (via conda) |
| scipy | (via conda, used for hierarchical clustering in ENRICHR heatmaps) |

## 24. Reproducibility Notes

- All R scripts set `set.seed(12345)` for reproducibility of stochastic operations (PCA, UMAP, permutation tests).
- `12_enrichr_pathway_analysis.py` sets `np.random.seed(12345)`, though ENRICHR API calls are inherently server-side and non-deterministic.
- `16_run_cellphonedb.py` does not set a random seed; `16b_run_cellphonedb_v4.py` uses `debug_seed=42`.
- Seurat 5.x column name mangling (prepending "g" to numeric IDs, replacing "_" with "-") was explicitly handled in pseudobulk aggregation scripts (11, 14) via a reverse-mapping function.
- Intermediate Seurat objects were saved as RDS checkpoints at each major processing step (preprocessing, annotation, macrophage subset) to enable re-entry without recomputation.
- All figures were saved in dual format (PDF vector + PNG raster at 300 DPI).
- The pipeline was executed on an HPC cluster via SLURM batch scripts with varying resource allocations (8 GB for downloads up to 64+ GB for GSVA on full matrices).

---

## Main Text Summary

Single-cell RNA-seq data from the Xu et al. (2024) Primary Breast Tumor Atlas were analyzed to characterize age-dependent transcriptional changes in hormone receptor-positive breast cancer. Treatment-naive HR+ patients (n ~ 36) were stratified into Young (age <= 50), MidAge (51-80), and Elderly (> 80) groups. Following quality control (200-6,000 genes, >= 400 UMIs, < 15% mitochondrial), SCTransform normalization, Harmony batch correction across 9 source studies, and atlas-derived cell type annotation (18 cell types), the dataset comprised approximately 115,000 cells. Cell type proportions were compared across age groups using Wilcoxon rank-sum tests with BH-FDR correction. GSVA pathway activity analysis was performed on estrogen-related gene sets from Hallmark, Reactome, WikiPathways, and GO collections, using both pseudo-bulk (CPM + log2, Gaussian kernel) and single-cell (Poisson kernel) approaches. PROGENy footprint-based pathway activity and WCSEA enrichment supplemented GSVA findings. Macrophage-focused pseudobulk differential expression (DESeq2) identified age-associated genes, with downstream pathway enrichment (hypergeometric test, Hallmark/BioCarta) and ENRICHR analysis (6 databases, per-library FDR correction). Cell-cell communication was assessed using both CellPhoneDB (v4/v5, 1,000-iteration permutation test) and CellChat (Secreted Signaling database, truncated mean method), comparing macrophage-immune cell interactions between age groups. Multi-cell-type GSVA with Welch t-tests on per-patient pseudo-bulk scores (curated-pathway FDR correction) identified age-differential pathway activity patterns across immune, epithelial, and stromal compartments.
