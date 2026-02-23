# Wu scRNA-seq Preprocessing Pipeline Design

**Date:** 2026-02-23
**Status:** Approved
**Author:** Claude (with user approval)

## Overview

Make the Wu et al. (GSE176078) single-cell RNA-seq analysis fully reproducible by adding data download and preprocessing scripts that were missing from the original codebase.

## Problem

The existing `analysis/04_human_scrnaseq/` pipeline expects a pre-processed Seurat object (`SeuratObj_GSE176078_ERpos_NewMeta_AfterQCSCT.rds`) that:
1. Does not exist in the repository
2. Has no documentation of how it was generated
3. Cannot be reproduced without the missing preprocessing code

## Solution

Add scripts `00a_download_geo.sh` and `00b_preprocess_seurat.R` that download raw data from GEO and process it following the exact methods from Sanghoon's publication.

## Methods (from publication)

> We compiled single-cell RNA sequencing (scRNA-seq) datasets from ten ER+ breast tumor samples, focusing exclusively on cells from pre-treatment patients. The processed scRNA-seq data, originally published by Wu et al. were obtained from the NCBI Gene Expression Omnibus (GEO) under accession number GSE176078. To analyze the data, we imported the matrix.mtx, features.tsv, and barcodes.tsv files into a Seurat object using the Read10X function in the Seurat R package (version 5.0.3). Our cell selection criteria included expression of fewer than 6,000 genes, with a minimum of 200 expressed genes and 400 UMI counts. Additionally, we filtered out cells with more than 15% of reads mapped to mitochondrial gene expression. To address technical variations and batch effects, we applied SCTransform to each Seurat object. Principal component analysis (PCA) was performed on the filtered feature-by-barcode matrix, and Uniform Manifold Approximation and Projection (UMAP) embeddings were based on the first 30 principal components. Overall, our processed dataset comprised 29,733 genes across 30,959 cells. For cell type annotation, we followed the major and subset cell types provided by Wu et al.

## Pipeline Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| MIN_FEATURES | 200 | Methods |
| MAX_FEATURES | 6000 | Methods |
| MIN_COUNTS | 400 | Methods |
| MAX_MT_PERCENT | 15% | Methods |
| PCA_DIMS | 30 | Methods |
| Seurat version | 5.0.3 | Methods |
| SCTransform method | glmGamPoi | Sanghoon's Rat workflow |

## Sample IDs (10 ER+ samples)

| Age Group | Sample IDs |
|-----------|------------|
| Young | CID3941, CID4530N, CID4535 |
| MidAge | CID4463, CID4040, CID4471, CID4461 |
| Elderly | CID3948, CID4067, CID4290A |

## File Structure

### New Scripts

```
analysis/04_human_scrnaseq/
├── 00a_download_geo.sh           # Download count matrices from GEO
├── 00b_preprocess_seurat.R       # QC, SCT, merge, PCA, UMAP
├── 01_load_subset_data.R         # (existing, minor edits)
├── 03_cell_type_annotation.R     # (existing, use Wu annotations)
...
```

### Data Directory

```
data/human_scrnaseq/
├── raw/                          # Downloaded from GEO
│   ├── CID3941/
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   └── .../
├── metadata/                     # Wu et al. cell type annotations
├── SeuratObj_GSE176078_ERpos_AfterQCSCT.rds   # Output of 00b
└── ClinicalData_Wu_scRNAseq_26p.txt
```

## Script Designs

### 00a_download_geo.sh

Downloads Wu et al. scRNA-seq count matrices from GEO:
- Source: `GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz`
- Extracts to `data/human_scrnaseq/raw/`
- SLURM job with 2hr time limit, 8GB memory

### 00b_preprocess_seurat.R

```r
# Pseudocode
for each sample_id in [10 ER+ samples]:
    counts <- Read10X(raw/{sample_id}/)
    obj <- CreateSeuratObject(counts)
    obj <- QC_filter(obj)  # nFeature 200-6000, nCount >400, MT <15%
    obj <- SCTransform(obj, method="glmGamPoi", vars.to.regress="percent.mt")

merged <- merge(all_objects)
merged <- RunPCA(merged, npcs=30)
merged <- RunUMAP(merged, dims=1:30)
# Add Wu et al. cell type annotations from metadata

saveRDS(merged, "SeuratObj_GSE176078_ERpos_AfterQCSCT.rds")
```

Expected output: ~30,959 cells × 29,733 genes

## Changes to Existing Scripts

| Script | Change |
|--------|--------|
| `01_load_subset_data.R` | Remove redundant PCA/UMAP calls |
| `02_harmony_integrate.R` | **DELETE** — not used per methods |
| `03_cell_type_annotation.R` | Load Wu annotations from metadata instead of re-annotating |
| `run_analysis.sh` | Add 00a, 00b to pipeline sequence |

## Data Provenance Documentation

Add to `data/human_scrnaseq/README.md`:

```markdown
## Wu et al. 2021 scRNA-seq Data

**Source:** GSE176078 (NCBI GEO)
**Reference:** Wu et al. 2021, Nature Genetics

### Processing Steps
1. Download: matrix.mtx, features.tsv, barcodes.tsv per sample
2. QC: nFeature 200-6000, nCount >400, percent.mt <15%
3. Normalization: SCTransform per sample (glmGamPoi)
4. Integration: Merge (no Harmony)
5. Dim reduction: PCA (30 dims) → UMAP

### Final Dataset
- 29,733 genes × 30,959 cells
- 10 ER+ samples (3 Young, 4 MidAge, 3 Elderly)
- Cell types: Original Wu et al. annotations
```

## Dependencies

Already in `environment.yml`:
- Seurat >= 5.0.3
- glmGamPoi (for SCTransform)

## Notes

- No Harmony integration per methods — SCTransform handles batch effects
- Cell type annotations come from Wu et al., not re-derived
- Clinical data file `ClinicalData_Wu_scRNAseq_26p.txt` already exists
