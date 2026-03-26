# Human scRNA-seq Data

## Primary Dataset: Xu et al. 2024 Breast Tumor Atlas

The main scRNA-seq reanalysis uses the Xu et al. 2024 primary breast tumor atlas, subsetted to HR-positive, treatment-naive patients (~115K cells across 36 patients from 9 studies).

**Reference:** Xu Y, et al. An integrated single-cell atlas of human breast cancers. Cancer Cell. 2024.

The atlas data must be obtained directly from the original authors. See `analysis/04_human_scrnaseq/00b_preprocess_seurat.R` for the expected input format (matrix.mtx, genes.tsv, barcodes.tsv, metadata.csv).

## Clinical Metadata

`ClinicalData_Wu_scRNAseq_26p.txt` contains clinical annotations used for age-group stratification.
