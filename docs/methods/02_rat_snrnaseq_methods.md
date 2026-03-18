# Computational Methods: Rat Mammary Tumor snRNA-seq (02_rat_snrnaseq)

## 1. Data Acquisition

Single-nucleus RNA sequencing (snRNA-seq) data were obtained from rat mammary tumors deposited under GEO accession GSE276758. The dataset comprised six samples (GSM8505435--GSM8505440, designated Lee_021924_Nuclei1 through Lee_021924_Nuclei6), with three biological replicates per age group: Aged (Nuclei1--3) and Young (Nuclei4--6). Raw 10X Chromium data were read using `Seurat::Read10X()` with `gene.column = 2`.

## 2. Preprocessing

### 2.1 Quality Control Filtering

Quality control was performed per sample prior to merging. Mitochondrial gene content was calculated using `PercentageFeatureSet()` with the pattern `"^MT-|^Mt-"`, and ribosomal protein content was calculated with the pattern `"^RP[SL]|^Rp[sl]"`. Nuclei were retained if they satisfied all of the following criteria:

- `nFeature_RNA > 200` (minimum detected genes)
- `nFeature_RNA < 6000` (maximum detected genes)
- `nCount_RNA > 400` (minimum UMI counts)
- `percent.mt < 15` (maximum mitochondrial content, %)

### 2.2 Doublet Detection and Removal

Doublet detection was performed per sample using DoubletFinder (v1.0). Each sample was independently normalized with `SCTransform()` (method = `"glmGamPoi"`, regressing `percent.mt`) and subjected to PCA (30 components) prior to doublet detection. The optimal `pK` parameter was determined via `paramSweep()` over PCs 1:30 with `sct = TRUE`, followed by `summarizeSweep()` (without ground truth, `GT = FALSE`) and `find.pK()`, selecting the `pK` value that maximized the BCmvn metric. The expected doublet rate was set to 5% (`doublet_rate = 0.05`), and `pN` was set to the default value of 0.25. Predicted doublets were removed, retaining only singlet-classified nuclei. After QC filtering and doublet removal, samples were merged using `Seurat::merge()` with sample-specific cell ID prefixes.

## 3. Normalization and Integration

### 3.1 SCTransform Normalization

The merged Seurat object was re-normalized using `SCTransform()` with the following parameters:

- `method = "glmGamPoi"` (fast negative binomial regression via the glmGamPoi package)
- `vars.to.regress = "percent.mt"` (regressing mitochondrial content)
- `seed.use = 12345`

### 3.2 Principal Component Analysis

PCA was computed on the SCTransform-normalized data using `RunPCA()` with `npcs = 50` and `seed.use = 12345`.

### 3.3 Harmony Integration

Batch correction across samples was performed using Harmony (`RunHarmony()`) with the following parameters:

- `group.by.vars = "orig.ident"` (integration by sample identity)
- `reduction.use = "pca"` (input PCA embeddings)
- `assay.use = "SCT"`
- `plot_convergence = FALSE`

Age group metadata (Young vs. Aged) was mapped to each nucleus from the external annotation file `Rat_scRNAseq_AgeGroup.txt` using exact matching on sample identifiers.

## 4. Clustering

### 4.1 UMAP Embedding

UMAP dimensionality reduction was computed on the Harmony-corrected embeddings using `RunUMAP()`:

- `reduction = "harmony"`
- `dims = 1:30`
- `seed.use = 12345`

### 4.2 Shared Nearest Neighbor Graph and Clustering

Neighbor finding was performed with `FindNeighbors()` using Harmony embeddings (dims 1:30). Two clustering workflows were implemented:

**Primary pipeline (03_cluster_annotate.R):** Clustering was performed at multiple resolutions (0.2, 0.4, 0.6, 0.8, 1.0) using `FindClusters()` with `random.seed = 12345` and the default Louvain algorithm. Resolution 0.4 was selected as the default for downstream analysis. The default `k.param` (20) was used for `FindNeighbors()`.

**scType pipeline (03b/03c):** Clustering was performed at resolution 1.5 with `k.param = 15` for `FindNeighbors()`, yielding approximately 40 clusters to match the original analysis granularity.

## 5. Cell Type Annotation

Three annotation strategies were implemented:

### 5.1 Marker-Based Cluster Scoring (03_cluster_annotate.R)

Average expression of canonical marker genes was computed per cluster using `AverageExpression()` on log-normalized RNA assay data (`NormalizeData()` with default parameters). For each cluster, mean expression across markers for each of eight cell types was calculated, and the cluster was assigned to the cell type with the highest mean score. The marker gene sets were:

| Cell Type | Markers |
|-----------|---------|
| Luminal | Krt8, Krt18, Krt19, Epcam, Cd24a |
| Basal | Krt14, Krt5, Krt17, Trp63, Acta2 |
| T cells | Cd3d, Cd3e, Cd3g, Cd4, Cd8a |
| B cells | Cd19, Cd79a, Cd79b, Ms4a1 |
| Macrophages | Cd68, Adgre1, Csf1r, Mrc1 |
| Fibroblasts | Col1a1, Col1a2, Col3a1, Dcn, Pdgfra |
| Endothelial | Pecam1, Cdh5, Vwf, Kdr |
| Adipocytes | Adipoq, Pparg, Fabp4, Lep |

Cluster-level markers were independently identified using `FindAllMarkers()` on the RNA assay with `test.use = "wilcox"`, `only.pos = TRUE`, `min.pct = 0.25`, and `logfc.threshold = 0.5`. The top 10 markers per cluster (ranked by `avg_log2FC`) were saved for verification.

### 5.2 scType Cell-by-Cell Annotation (03b_sctype_annotate.R)

Cell-level annotation was performed using a custom implementation of the scType algorithm (Ianevski et al., 2022). A custom rat mammary tissue marker database was defined with 14 cell types: CancerEpithelial, Myoepithelial, Fibroblast, Endothelial, Monocyte, Macrophage, M1_Macrophage, M2_Macrophage, DendriticCell, NKcell, CD4Tcell, CD8Tcell, Treg, and NaiveTcell. Each cell type was defined by both positive markers (expected to be expressed) and negative markers (expected to be absent). For each cell, a score was computed per cell type as:

`score = sum(scaled_expression[positive_markers]) / sqrt(n_positive) - sum(scaled_expression[negative_markers]) / sqrt(n_negative)`

Scaled expression data from the SCT assay (`scale.data` layer) were used as input. Each cell was assigned to the cell type with the maximum positive score; cells with all non-positive scores were classified as "Unknown." Confidence scores were computed as the difference between the top two scores. Detailed cell types were mapped to broader categories: Monocyte/Macrophage/M1/M2 to "Myeloid"; NKcell/CD4Tcell/CD8Tcell/Treg/NaiveTcell to "NKTcell."

### 5.3 scType Cluster-Level Annotation (03c_sctype_cluster.R)

Cell-level scType scores (computed as in Section 5.2) were aggregated by cluster: for each cluster, the mean score per cell type was calculated by summing cell-level scores and dividing by the number of cells in the cluster. Each cluster was assigned to the cell type with the highest mean score. Results were compared against the original manual annotations from the source analysis, which assigned 40 clusters at resolution 1.5 based on visual inspection of dot plots.

### 5.4 Method Comparison

A systematic comparison of annotation methods (`method_comparison.R`, `spot_check_markers.R`) validated immune cell annotations (Ptprc/CD45 expression in cells called immune) and identified potential misclassification of endothelial cells that expressed epithelial markers (Epcam, Krt18). The comparison recommended a hybrid approach: scType for immune cell identification and the original manual strategy for epithelial subtype annotation.

## 6. Differential Expression

Differential expression testing between Aged and Young groups was performed per cell type using a pseudobulk DESeq2 approach. Raw counts are aggregated per sample using `AggregateExpression()`, with DESeq2 `~ group` design and Aged vs Young contrast. This properly treats biological replicates (n=3 per group) as the unit of analysis. Cell types with fewer than 10 cells in either group were excluded. P-values were corrected for multiple testing using the Benjamini-Hochberg method (`p.adjust(method = "BH")`) applied within each cell type. Genes were classified as significant at FDR < 0.05 and further annotated by direction (upregulated or downregulated in Aged).

## 7. Differential Abundance

Differential abundance of cell types between Young and Aged groups was tested using the propeller method from the speckle package (`speckle::propeller()`). Propeller performs an empirical Bayes moderated test on arcsin-square-root-transformed cell type proportions, properly accounting for the compositional nature of single-cell data and using biological replicates (samples) as the unit of analysis.

Two levels of cell type granularity were tested:

1. **Main cell types** (`CellTypeByMarker_RatsnRNAseq`): broad categories (e.g., CancerEpithelial, Myeloid, NKTcell, Fibroblast, Endothelial)
2. **Subtypes** (`CellTypeMacroTcell_RatsnRNAseq`): finer-grained types from scType annotation

P-values from propeller were additionally corrected using the Benjamini-Hochberg method. Results at both levels were combined for reporting.

## 8. Visualization

The following visualizations were generated:

| Plot | Description | Output File |
|------|-------------|-------------|
| QC violin plots | nFeature_RNA, nCount_RNA, percent.mt per sample | `outputs/qc_plots.pdf` |
| Feature scatter | nCount_RNA vs. nFeature_RNA | `outputs/qc_plots.pdf` |
| Elbow plot | PCA variance explained (50 PCs) | `outputs/integration_plots.pdf` |
| UMAP by sample/age group | Integration quality assessment | `outputs/integration_plots.pdf` |
| UMAP by cluster | At resolutions 0.2, 0.4, 0.6, 0.8, 1.0 | `outputs/clustering_plots.pdf` |
| UMAP by cell type | Annotated cell type labels | `outputs/clustering_plots.pdf` |
| Marker heatmap | Scaled average expression per cluster (pheatmap) | `outputs/marker_heatmap.pdf` |
| scType UMAP | Detailed and main cell type annotations | `figures/by_analysis/rat_snrnaseq/sctype_annotation_umap.pdf` |
| scType cluster UMAP | Cluster-level scType annotations | `figures/by_analysis/rat_snrnaseq/sctype_cluster_annotation.pdf` |
| Volcano plots | Per-cell-type DE results, top 10 FDR genes labeled (ggrepel) | `figures/by_analysis/rat_snrnaseq/DE_volcano_plots.pdf` |
| DA stacked bar | Cell type proportions per sample, faceted by age group | `figures/by_analysis/rat_snrnaseq/DA_proportion_plots.pdf` |
| DA boxplot | Cell type proportions (%) by age group with jittered points | `figures/by_analysis/rat_snrnaseq/DA_proportion_plots.pdf` |
| DA log2FC bar | Propeller log2 fold change by cell type, FDR-colored | `figures/by_analysis/rat_snrnaseq/DA_proportion_plots.pdf` |

All figure panels were saved in both PDF (vector) and PNG (300 DPI) formats, with SVG versions for selected panels.

## 9. Software Versions

The analysis was performed using the `erp_brca_aging` conda environment with the following key software:

| Package | Version | Purpose |
|---------|---------|---------|
| R | 4.4.1 | Base language |
| Seurat | >= 5.0 | snRNA-seq analysis framework |
| SeuratObject | (bundled) | Seurat data structures |
| harmony | (latest conda) | Batch correction / integration |
| glmGamPoi | (Bioconductor) | Fast SCTransform backend |
| DoubletFinder | (latest) | Doublet detection |
| speckle | (latest conda) | Differential abundance (propeller) |
| pheatmap | (latest) | Heatmap visualization |
| ggplot2 | (tidyverse bundle) | Plotting |
| ggrepel | (latest conda) | Non-overlapping text labels |
| patchwork | (latest conda) | Multi-panel figure composition |
| HGNChelper | (latest) | Gene symbol validation (scType) |
| openxlsx | (latest) | Excel file I/O (scType) |
| dplyr / tidyr | (tidyverse) | Data wrangling |

## 10. Reproducibility Notes

- **Random seed:** All scripts set `set.seed(12345)` at the top. Seed was explicitly passed to `SCTransform()` (`seed.use = 12345`), `RunPCA()` (`seed.use = 12345`), `RunUMAP()` (`seed.use = 12345`), and `FindClusters()` (`random.seed = 12345`).
- **Parallelization:** `future.globals.maxSize` was set to 4 GB (`4 * 1024^3`). No explicit `plan()` call was made, so processing ran sequentially by default.
- **Checkpoints:** Intermediate Seurat objects were saved as RDS files at each major step: `seurat_qc_filtered.rds` (post-QC), `seurat_integrated.rds` (post-integration), `seurat_annotated.rds` (post-clustering/annotation), with additional checkpoints for scType variants.
- **SLURM parameters:** The pipeline was executed via SLURM (`run_analysis.sbatch`) with 64 GB memory, 8 CPUs, 12-hour wall time, on the `htc` partition.

---

## Main Text Summary

Single-nucleus RNA sequencing data from six NMU-induced rat mammary tumors (three Aged, three Young; GEO: GSE276758) were processed using Seurat v5. Nuclei were filtered by gene count (200--6,000), UMI count (>400), and mitochondrial content (<15%), and doublets were removed using DoubletFinder (5% expected rate, per-sample pK optimization). Data were normalized with SCTransform (glmGamPoi, regressing mitochondrial percentage) and integrated across samples using Harmony. UMAP embedding and shared nearest neighbor clustering (Louvain algorithm, resolution 0.4) were computed on the first 30 Harmony-corrected components. Cell types were annotated by scoring clusters against canonical mammary tissue marker gene sets; a parallel scType-based annotation confirmed immune cell identification. Differential expression between age groups was tested per cell type using pseudobulk DESeq2 (BH-corrected FDR < 0.05), and differential abundance was assessed using propeller (speckle package) with Benjamini-Hochberg correction. All analyses used random seed 12345.
