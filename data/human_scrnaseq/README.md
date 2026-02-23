# Wu et al. 2021 scRNA-seq Data

**Source:** GSE176078 (NCBI GEO)
**Reference:** Wu SZ, et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet. 2021;53(9):1334-1347.

## Data Acquisition

### Option 1: Run download script (recommended)

```bash
# Download raw count matrices
sbatch analysis/04_human_scrnaseq/00a_download_geo.sh

# Or run interactively
bash analysis/04_human_scrnaseq/00a_download_geo.sh
```

### Option 2: Manual download

1. Go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
2. Download `GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz`
3. Extract to `data/human_scrnaseq/raw/`

## Processing Steps

Per published methods:

1. **Load:** `Read10X()` to import matrix.mtx, features.tsv, barcodes.tsv
2. **QC filtering:**
   - nFeature_RNA > 200 and < 6000
   - nCount_RNA > 400
   - percent.mt < 15%
3. **Normalization:** SCTransform per sample (glmGamPoi, regress percent.mt)
4. **Merge:** Combine all samples
5. **Dimensionality reduction:** PCA (30 dims) → UMAP

## Samples Used (10 ER+ samples)

| Age Group | Sample IDs |
|-----------|------------|
| Young (≤50) | CID3941, CID4530N, CID4535 |
| MidAge (51-80) | CID4463, CID4040, CID4471, CID4461 |
| Elderly (>80) | CID3948, CID4067, CID4290A |

## Output Files

| File | Description |
|------|-------------|
| `SeuratObj_GSE176078_ERpos_AfterQCSCT.rds` | Processed Seurat object (all 10 samples) |
| `ClinicalData_Wu_scRNAseq_26p.txt` | Clinical metadata |

## Expected Dimensions

- ~30,959 cells
- ~29,733 genes

## Cell Type Annotations

Cell type annotations follow the original Wu et al. classification:
- CellTypeMajor: Main cell categories
- CellTypeMinor: Subcategories for immune cells
