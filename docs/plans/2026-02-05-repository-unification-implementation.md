# Repository Unification Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Unify four GitHub repositories into ERpBRCA_OlderWomen with biostatistical corrections and full reproducibility.

**Architecture:** Extract embedded R code from README files into standalone scripts, organize into analysis directories with conda environments, implement FDR correction and DE/DA testing, create comparison framework for original vs corrected results.

**Tech Stack:** R 4.4.1, Seurat, DESeq2, GSVA, PROGENy, Harmony, conda, bash

---

## Phase 1: Repository Structure Setup

### Task 1.1: Create Directory Structure

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/.gitkeep`
- Create: `analysis/02_rat_snrnaseq/.gitkeep`
- Create: `analysis/03_rat_wes/.gitkeep`
- Create: `data/human_bulk_rnaseq/.gitkeep`
- Create: `data/rat_snrnaseq/.gitkeep`
- Create: `data/rat_wes/.gitkeep`
- Create: `figures/by_analysis/.gitkeep`
- Create: `figures/manuscript/.gitkeep`
- Create: `results/original/.gitkeep`
- Create: `results/corrected/.gitkeep`
- Create: `results/comparison/.gitkeep`
- Create: `scripts/.gitkeep`
- Create: `envs/.gitkeep`

**Step 1: Create all directories**

```bash
mkdir -p analysis/01_human_bulk_rnaseq
mkdir -p analysis/02_rat_snrnaseq
mkdir -p analysis/03_rat_wes
mkdir -p data/human_bulk_rnaseq/raw
mkdir -p data/human_bulk_rnaseq/external
mkdir -p data/rat_snrnaseq/raw
mkdir -p data/rat_snrnaseq/metadata
mkdir -p data/rat_wes/raw
mkdir -p figures/by_analysis/human_bulk_rnaseq
mkdir -p figures/by_analysis/rat_snrnaseq
mkdir -p figures/by_analysis/rat_wes
mkdir -p figures/manuscript/Fig1
mkdir -p figures/manuscript/Fig2
mkdir -p figures/manuscript/Fig3
mkdir -p figures/manuscript/FigS1
mkdir -p figures/manuscript/FigS2
mkdir -p results/original/human_bulk_rnaseq
mkdir -p results/original/rat_snrnaseq
mkdir -p results/corrected/human_bulk_rnaseq
mkdir -p results/corrected/rat_snrnaseq
mkdir -p results/comparison/figures
mkdir -p scripts
mkdir -p envs
```

**Step 2: Add .gitkeep files to empty directories**

```bash
find . -type d -empty -not -path "./.git/*" -exec touch {}/.gitkeep \;
```

**Step 3: Commit**

```bash
git add -A
git commit -m "chore: create unified repository directory structure"
```

---

### Task 1.2: Create Master Conda Environment

**Files:**
- Create: `environment.yml`

**Step 1: Write environment file**

```yaml
name: erp_brca_aging
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  # R base
  - r-base=4.4.1
  - r-essentials

  # R core packages
  - r-tidyverse
  - r-data.table
  - r-ggplot2
  - r-patchwork
  - r-rmarkdown
  - r-knitr
  - r-ggrepel
  - r-gridextra
  - r-rcolorbrewer
  - r-pheatmap

  # Bioconductor packages
  - bioconductor-deseq2
  - bioconductor-gsva
  - bioconductor-progeny
  - bioconductor-msigdbr
  - bioconductor-qusage

  # Seurat ecosystem
  - r-seurat>=5.0
  - r-seuratobject
  - r-harmony
  - bioconductor-glmgampoi
  - bioconductor-singlecellexperiment
  - bioconductor-singler
  - r-dittoseq

  # Statistical packages
  - r-speckle           # For differential abundance
  - r-lme4

  # Python for WES and utilities
  - python=3.10
  - pandas
  - numpy
  - jupyter
  - matplotlib
  - seaborn

  # System tools
  - pandoc
  - wget
  - curl

  # PDF generation
  - r-pagedown
  - chromium
```

**Step 2: Commit**

```bash
git add environment.yml
git commit -m "chore: add master conda environment file"
```

---

### Task 1.3: Create Analysis-Specific Conda Environments

**Files:**
- Create: `envs/bulk_rnaseq.yml`
- Create: `envs/snrnaseq.yml`
- Create: `envs/wes.yml`

**Step 1: Write bulk RNA-seq environment**

```yaml
# envs/bulk_rnaseq.yml
name: erp_bulk_rnaseq
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - r-base=4.4.1
  - r-tidyverse
  - r-data.table
  - r-ggplot2
  - r-patchwork
  - r-rmarkdown
  - bioconductor-deseq2
  - bioconductor-gsva
  - bioconductor-progeny
  - bioconductor-msigdbr
  - bioconductor-qusage
  - r-ggrepel
  - r-pheatmap
  - r-rcolorbrewer
```

**Step 2: Write snRNA-seq environment**

```yaml
# envs/snrnaseq.yml
name: erp_snrnaseq
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - r-base=4.4.1
  - r-tidyverse
  - r-data.table
  - r-ggplot2
  - r-patchwork
  - r-seurat>=5.0
  - r-seuratobject
  - r-harmony
  - bioconductor-glmgampoi
  - bioconductor-singler
  - r-dittoseq
  - r-speckle
  - r-doubletfinder
  - r-ggrepel
  - r-gridextra
```

**Step 3: Copy and adapt WES environment from RatAgingWES**

```bash
cp ../RatAgingWES/aging_wes.yml envs/wes.yml
```

**Step 4: Commit**

```bash
git add envs/
git commit -m "chore: add analysis-specific conda environments"
```

---

### Task 1.4: Create Data Download Scripts

**Files:**
- Create: `data/download_all.sh`
- Create: `data/human_bulk_rnaseq/download.sh`
- Create: `data/rat_snrnaseq/download.sh`
- Create: `data/rat_wes/download.sh`
- Create: `data/README.md`

**Step 1: Write master download script**

```bash
#!/bin/bash
# data/download_all.sh
# Master script to download all data from GEO
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Downloading all GEO data ==="

echo "[1/3] Downloading human bulk RNA-seq (GSE276755)..."
bash "$SCRIPT_DIR/human_bulk_rnaseq/download.sh"

echo "[2/3] Downloading rat snRNA-seq (GSE276758)..."
bash "$SCRIPT_DIR/rat_snrnaseq/download.sh"

echo "[3/3] Downloading rat WES (GSE276759)..."
bash "$SCRIPT_DIR/rat_wes/download.sh"

echo "=== All downloads complete ==="
```

**Step 2: Write human bulk RNA-seq download script**

```bash
#!/bin/bash
# data/human_bulk_rnaseq/download.sh
# Downloads GSE276755 human bulk RNA-seq data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"
EXT_DIR="$SCRIPT_DIR/external"

mkdir -p "$RAW_DIR" "$EXT_DIR"

echo "Downloading GSE276755 data..."

# Check if files already exist
if [[ -f "$RAW_DIR/HumanERpAge_39404g168s_FeatureCount.txt" ]]; then
    echo "Data already downloaded, skipping..."
    exit 0
fi

# Download from GEO supplementary files
# Note: Update these URLs with actual GEO links when available
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE276nnn/GSE276755/suppl"

wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_39404g168s_FeatureCount.txt.gz" || echo "Warning: Could not download count file"
wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_39404g168s_TPMlog2.txt.gz" || echo "Warning: Could not download TPM file"
wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_BulkRNAseq_SampleInformation.txt.gz" || echo "Warning: Could not download sample info"

# Decompress if downloaded
gunzip -k "$RAW_DIR"/*.gz 2>/dev/null || true

# Download gene annotation from NCBI
echo "Downloading gene annotation..."
wget -q -P "$EXT_DIR" "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
gunzip -k "$EXT_DIR/Homo_sapiens.gene_info.gz" 2>/dev/null || true
mv "$EXT_DIR/Homo_sapiens.gene_info" "$EXT_DIR/Homo_sapiens.gene_info.txt" 2>/dev/null || true

echo "Human bulk RNA-seq download complete."
```

**Step 3: Write rat snRNA-seq download script**

```bash
#!/bin/bash
# data/rat_snrnaseq/download.sh
# Downloads GSE276758 rat snRNA-seq data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"
META_DIR="$SCRIPT_DIR/metadata"

mkdir -p "$RAW_DIR" "$META_DIR"

echo "Downloading GSE276758 data..."

# Check if files already exist
if [[ -d "$RAW_DIR/Lee_021924_Nuclei1" ]]; then
    echo "Data already downloaded, skipping..."
    exit 0
fi

# Note: 10X data from GEO requires specific download structure
# Update with actual GEO supplementary file URLs
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE276nnn/GSE276758/suppl"

# Download the tar archive of 10X matrices
wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276758_RAW.tar" || echo "Warning: Could not download 10X data"

# Extract if downloaded
if [[ -f "$RAW_DIR/GSE276758_RAW.tar" ]]; then
    tar -xf "$RAW_DIR/GSE276758_RAW.tar" -C "$RAW_DIR"
    rm "$RAW_DIR/GSE276758_RAW.tar"
fi

# Create age group annotation file
cat > "$META_DIR/Rat_scRNAseq_AgeGroup.txt" << 'EOF'
CaseID	AgeGroup
Lee_021924_Nuclei1	Aged
Lee_021924_Nuclei2	Aged
Lee_021924_Nuclei3	Aged
Lee_021924_Nuclei4	Young
Lee_021924_Nuclei5	Young
Lee_021924_Nuclei6	Young
EOF

echo "Rat snRNA-seq download complete."
```

**Step 4: Write rat WES download script**

```bash
#!/bin/bash
# data/rat_wes/download.sh
# Downloads GSE276759 rat WES data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"

mkdir -p "$RAW_DIR"

echo "Downloading GSE276759 data..."

# Check if files already exist
if [[ -f "$RAW_DIR/.downloaded" ]]; then
    echo "Data already downloaded, skipping..."
    exit 0
fi

# Note: WES data may need SRA download
# Update with actual GEO/SRA links
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE276nnn/GSE276759/suppl"

wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276759_processed_data.tar.gz" || echo "Warning: Could not download WES data"

# Extract if downloaded
if [[ -f "$RAW_DIR/GSE276759_processed_data.tar.gz" ]]; then
    tar -xzf "$RAW_DIR/GSE276759_processed_data.tar.gz" -C "$RAW_DIR"
fi

touch "$RAW_DIR/.downloaded"
echo "Rat WES download complete."
```

**Step 5: Write data README**

```markdown
# Data Directory

This directory contains all input data for the analyses.

## GEO Accessions

| Dataset | GEO Accession | Description |
|---------|---------------|-------------|
| Human Bulk RNA-seq | GSE276755 | 168 samples (83 tumor, 85 adjacent normal) |
| Rat snRNA-seq | GSE276758 | 6 samples (3 Young, 3 Aged) |
| Rat WES | GSE276759 | Whole exome sequencing of rat tumors |

## Download Instructions

Run the master download script:

```bash
./download_all.sh
```

Or download individual datasets:

```bash
./human_bulk_rnaseq/download.sh
./rat_snrnaseq/download.sh
./rat_wes/download.sh
```

## Directory Structure

```
data/
├── human_bulk_rnaseq/
│   ├── raw/                    # Count and TPM matrices
│   └── external/               # Gene annotation files
├── rat_snrnaseq/
│   ├── raw/                    # 10X filtered matrices per sample
│   └── metadata/               # Sample annotations
└── rat_wes/
    └── raw/                    # VEP and CNVKit outputs
```
```

**Step 6: Make scripts executable and commit**

```bash
chmod +x data/download_all.sh
chmod +x data/human_bulk_rnaseq/download.sh
chmod +x data/rat_snrnaseq/download.sh
chmod +x data/rat_wes/download.sh
git add data/
git commit -m "feat: add GEO data download scripts"
```

---

## Phase 2: Extract and Refactor Code

### Task 2.1: Extract Human Bulk RNA-seq Preprocessing Script

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/01_preprocess.R`

**Step 1: Write preprocessing script**

```r
#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/01_preprocess.R
# Preprocess human bulk RNA-seq data: read counts, filter genes, normalize
#
# Inputs:
#   - data/human_bulk_rnaseq/raw/HumanERpAge_39404g168s_FeatureCount.txt
#   - data/human_bulk_rnaseq/raw/HumanERpAge_BulkRNAseq_SampleInformation.txt
#   - data/human_bulk_rnaseq/external/Homo_sapiens.gene_info.txt
#
# Outputs:
#   - analysis/01_human_bulk_rnaseq/outputs/dds_norm_AllAgeGroup.rds
#   - analysis/01_human_bulk_rnaseq/outputs/vst_normalized_matrix.rds
#   - analysis/01_human_bulk_rnaseq/outputs/qc_plots.pdf

# Set random seed for reproducibility
set.seed(12345)

# Load libraries
suppressPackageStartupMessages({
  library(tidyr)
  library(DESeq2)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
})

# Define paths relative to project root
project_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")
output_dir <- file.path(dirname(sys.frame(1)$ofile), "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
count_file <- file.path(data_dir, "raw/HumanERpAge_39404g168s_FeatureCount.txt")
sample_file <- file.path(data_dir, "raw/HumanERpAge_BulkRNAseq_SampleInformation.txt")
gene_annot_file <- file.path(data_dir, "external/Homo_sapiens.gene_info.txt")

cat("=== Human Bulk RNA-seq Preprocessing ===\n")
cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Read Gene Annotation
# -----------------------------------------------------------------------------
cat("Step 1: Reading gene annotation...\n")
gene_annot <- fread(gene_annot_file, header = TRUE, stringsAsFactors = FALSE)
gene_annot_filtered <- gene_annot %>%
  dplyr::select(Symbol, type_of_gene) %>%
  dplyr::filter(type_of_gene == "protein-coding", !grepl("^LOC\\d+", Symbol))
cat("  Protein-coding genes:", nrow(gene_annot_filtered), "\n")

# -----------------------------------------------------------------------------
# Step 2: Read Count Data
# -----------------------------------------------------------------------------
cat("Step 2: Reading count data...\n")
count_data <- fread(count_file, header = TRUE, stringsAsFactors = FALSE)
colnames(count_data)[1] <- "GeneSymb"
colnames(count_data) <- gsub("_LEE(.*)", "", colnames(count_data))
cat("  Raw genes:", nrow(count_data), "\n")
cat("  Samples:", ncol(count_data) - 1, "\n")

# Filter to protein-coding genes
count_data_filtered <- count_data %>%
  dplyr::filter(GeneSymb %in% gene_annot_filtered$Symbol)
cat("  Protein-coding genes in data:", nrow(count_data_filtered), "\n")

# -----------------------------------------------------------------------------
# Step 3: Read Sample Annotation
# -----------------------------------------------------------------------------
cat("Step 3: Reading sample annotation...\n")
sample_annot <- fread(sample_file, stringsAsFactors = FALSE, header = TRUE)
colnames(sample_annot)[2] <- "SampleName"
sample_annot <- sample_annot %>%
  dplyr::mutate(
    AgeRange_Group = paste0(AgeRange, "_", Group),
    SampleNameGroup = paste0(SampleName, "_", AgeRange, Group)
  )
cat("  Samples:", nrow(sample_annot), "\n")
cat("  Age groups:", paste(names(table(sample_annot$AgeRange)), collapse = ", "), "\n")
cat("  Sample types:", paste(names(table(sample_annot$Group)), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Step 4: Prepare DESeq2 Input
# -----------------------------------------------------------------------------
cat("Step 4: Preparing DESeq2 input...\n")

# Transpose and join
count_data_t <- count_data_filtered %>%
  tibble::column_to_rownames("GeneSymb") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleName")

count_annot <- dplyr::inner_join(
  sample_annot[, c("SampleName", "SampleNameGroup", "AgeRange", "Group")],
  count_data_t,
  by = "SampleName"
) %>%
  dplyr::select(-SampleName)

# Create count matrix (genes x samples)
count_matrix <- count_annot %>%
  tibble::column_to_rownames("SampleNameGroup") %>%
  dplyr::select(-AgeRange, -Group) %>%
  t() %>%
  as.data.frame()

# Create metadata
metadata <- data.frame(
  SampleNameGroup = colnames(count_matrix),
  row.names = colnames(count_matrix)
)
metadata <- metadata %>%
  dplyr::mutate(
    AgeRange = gsub(".*_(.*)Tumor.*", "\\1", SampleNameGroup),
    Group = ifelse(grepl("TumorAdj", SampleNameGroup), "TumorAdj", "Tumor")
  )

cat("  Count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# -----------------------------------------------------------------------------
# Step 5: Create DESeq2 Object with CORRECTED Design
# -----------------------------------------------------------------------------
cat("Step 5: Creating DESeq2 object...\n")

# BIOSTATISTICAL FIX: Use proper design formula instead of ~1
# This accounts for AgeRange and Group effects during normalization
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_matrix),
  colData = metadata,
  design = ~ AgeRange + Group  # CORRECTED: was design = ~1
)

cat("  DESeq2 object created with design: ~ AgeRange + Group\n")

# -----------------------------------------------------------------------------
# Step 6: Variance Stabilizing Transformation
# -----------------------------------------------------------------------------
cat("Step 6: Applying variance stabilizing transformation...\n")
vst_data <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_matrix <- assay(vst_data) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneSymb")

# Remove duplicate genes
vst_matrix_nodup <- vst_matrix %>%
  dplyr::filter(!duplicated(GeneSymb))
cat("  VST matrix:", nrow(vst_matrix_nodup), "genes\n")

# -----------------------------------------------------------------------------
# Step 7: Generate QC Plots
# -----------------------------------------------------------------------------
cat("Step 7: Generating QC plots...\n")

pdf(file.path(output_dir, "qc_plots.pdf"), width = 12, height = 10)

# PCA plot
pca_data <- plotPCA(vst_data, intgroup = c("AgeRange", "Group"), returnData = TRUE)
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = AgeRange, shape = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  ggtitle("PCA: Colored by Age Range, Shaped by Sample Type")
print(pca_plot)

# Sample distance heatmap
sample_dists <- dist(t(assay(vst_data)))
sample_dist_matrix <- as.matrix(sample_dists)
pheatmap::pheatmap(
  sample_dist_matrix,
  main = "Sample Distance Matrix",
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists
)

# Library size distribution
lib_sizes <- colSums(counts(dds))
lib_df <- data.frame(
  Sample = names(lib_sizes),
  LibrarySize = lib_sizes,
  Group = metadata$Group,
  AgeRange = metadata$AgeRange
)
lib_plot <- ggplot(lib_df, aes(x = reorder(Sample, LibrarySize), y = LibrarySize / 1e6, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  labs(x = "Sample", y = "Library Size (millions)", title = "Library Size Distribution") +
  theme(axis.text.y = element_text(size = 6))
print(lib_plot)

dev.off()
cat("  QC plots saved to:", file.path(output_dir, "qc_plots.pdf"), "\n")

# -----------------------------------------------------------------------------
# Step 8: Save Outputs
# -----------------------------------------------------------------------------
cat("Step 8: Saving outputs...\n")
saveRDS(dds, file.path(output_dir, "dds_norm_AllAgeGroup.rds"))
saveRDS(vst_matrix_nodup, file.path(output_dir, "vst_normalized_matrix.rds"))
saveRDS(sample_annot, file.path(output_dir, "sample_annotation.rds"))

cat("\n=== Preprocessing complete ===\n")
cat("Outputs saved to:", output_dir, "\n")
```

**Step 2: Commit**

```bash
git add analysis/01_human_bulk_rnaseq/01_preprocess.R
git commit -m "feat(bulk-rnaseq): add preprocessing script with corrected DESeq2 design"
```

---

### Task 2.2: Extract Human Bulk RNA-seq GSVA Script

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/02_run_gsva.R`

**Step 1: Write GSVA script**

```r
#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/02_run_gsva.R
# Run Gene Set Variation Analysis on estrogen-related pathways
#
# Inputs:
#   - analysis/01_human_bulk_rnaseq/outputs/vst_normalized_matrix.rds
#   - data/human_bulk_rnaseq/external/SuppleTable1_UpregulatedByE1_NotE2_407g.txt
#
# Outputs:
#   - analysis/01_human_bulk_rnaseq/outputs/gsva_estrogen_pathways.rds
#   - analysis/01_human_bulk_rnaseq/outputs/gsva_heatmap.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# Define paths
script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")

cat("=== GSVA Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Normalized Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading normalized data...\n")
vst_matrix <- readRDS(file.path(output_dir, "vst_normalized_matrix.rds"))
vst_mat <- vst_matrix %>%
  tibble::column_to_rownames("GeneSymb") %>%
  as.matrix()
cat("  Matrix:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Load Gene Sets
# -----------------------------------------------------------------------------
cat("Step 2: Loading gene sets...\n")

# E1-upregulated genes
e1_file <- file.path(data_dir, "external/SuppleTable1_UpregulatedByE1_NotE2_407g.txt")
if (file.exists(e1_file)) {
  e1_genes <- fread(e1_file, header = TRUE, stringsAsFactors = FALSE)
  colnames(e1_genes) <- gsub(" ", "", colnames(e1_genes))
  e1_genes$Gene <- gsub("-.*", "", e1_genes$Gene)
  e1_genes <- e1_genes %>% filter(!duplicated(Gene))
  e1_gene_list <- list(E1UpRegGene = e1_genes$Gene)
  cat("  E1-upregulated genes:", length(e1_gene_list$E1UpRegGene), "\n")
} else {
  e1_gene_list <- list()
  cat("  Warning: E1 gene file not found\n")
}

# Hallmark estrogen pathways
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
hallmark_estrogen <- hallmark_list[c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)]

# Reactome estrogen pathway
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
reactome_list <- split(reactome_sets$gene_symbol, reactome_sets$gs_name)
reactome_estrogen <- reactome_list["REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION"]

# WikiPathways
wiki_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "WIKIPATHWAYS")
wiki_list <- split(wiki_sets$gene_symbol, wiki_sets$gs_name)
wiki_estrogen <- wiki_list["WP_ESTROGEN_SIGNALING_PATHWAY"]

# GO BP estrogen pathways
gobp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gobp_list <- split(gobp_sets$gene_symbol, gobp_sets$gs_name)
gobp_estrogen <- gobp_list[c(
  "GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS"
)]

# Combine all estrogen pathways
estrogen_pathways <- c(hallmark_estrogen, e1_gene_list, reactome_estrogen, wiki_estrogen, gobp_estrogen)
estrogen_pathways <- estrogen_pathways[!sapply(estrogen_pathways, is.null)]
cat("  Total pathways:", length(estrogen_pathways), "\n")

# -----------------------------------------------------------------------------
# Step 3: Run GSVA
# -----------------------------------------------------------------------------
cat("Step 3: Running GSVA...\n")

gsva_result <- gsva(
  gsvaParam(
    vst_mat,
    estrogen_pathways,
    kcdf = "Gaussian",  # Appropriate for continuous (VST) data
    maxDiff = TRUE
  )
)

cat("  GSVA result:", nrow(gsva_result), "pathways x", ncol(gsva_result), "samples\n")

# -----------------------------------------------------------------------------
# Step 4: Generate Heatmap
# -----------------------------------------------------------------------------
cat("Step 4: Generating heatmap...\n")

# Extract sample metadata from column names
sample_info <- data.frame(
  Sample = colnames(gsva_result),
  AgeRange = gsub(".*_(.*)Tumor.*", "\\1", colnames(gsva_result)),
  Group = ifelse(grepl("TumorAdj", colnames(gsva_result)), "TumorAdj", "Tumor")
)
rownames(sample_info) <- sample_info$Sample

# Annotation colors
ann_colors <- list(
  AgeRange = c(Young = "#4DAF4A", Middle = "#377EB8", Elderly = "#E41A1C"),
  Group = c(Tumor = "#984EA3", TumorAdj = "#FF7F00")
)

pdf(file.path(output_dir, "gsva_heatmap.pdf"), width = 14, height = 8)
pheatmap(
  gsva_result,
  annotation_col = sample_info[, c("AgeRange", "Group")],
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  main = "GSVA Estrogen Pathway Scores",
  scale = "row"
)
dev.off()

# -----------------------------------------------------------------------------
# Step 5: Save Outputs
# -----------------------------------------------------------------------------
cat("Step 5: Saving outputs...\n")
saveRDS(gsva_result, file.path(output_dir, "gsva_estrogen_pathways.rds"))
saveRDS(estrogen_pathways, file.path(output_dir, "estrogen_pathway_genesets.rds"))

cat("\n=== GSVA complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/01_human_bulk_rnaseq/02_run_gsva.R
git commit -m "feat(bulk-rnaseq): add GSVA analysis script"
```

---

### Task 2.3: Extract Human Bulk RNA-seq PROGENy Script

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/03_run_progeny.R`

**Step 1: Write PROGENy script**

```r
#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/03_run_progeny.R
# Run PROGENy pathway activity analysis
#
# Inputs:
#   - data/human_bulk_rnaseq/raw/HumanERpAge_39404g168s_TPMlog2.txt
#   - analysis/01_human_bulk_rnaseq/outputs/sample_annotation.rds
#
# Outputs:
#   - analysis/01_human_bulk_rnaseq/outputs/progeny_pathway_activity.rds

# BIOSTATISTICAL FIX: Set random seed for reproducibility
set.seed(12345)

suppressPackageStartupMessages({
  library(progeny)
  library(dplyr)
  library(data.table)
  library(tibble)
})

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")

cat("=== PROGENy Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load TPM Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading TPM data...\n")
tpm_file <- file.path(data_dir, "raw/HumanERpAge_39404g168s_TPMlog2.txt")
tpm_data <- fread(tpm_file, stringsAsFactors = FALSE, header = TRUE)
colnames(tpm_data)[1] <- "GeneSymb"
colnames(tpm_data) <- gsub("_LEE.*", "", colnames(tpm_data))

# Load gene annotation for filtering
gene_annot <- fread(
  file.path(data_dir, "external/Homo_sapiens.gene_info.txt"),
  header = TRUE, stringsAsFactors = FALSE
)
prot_coding <- gene_annot %>%
  filter(type_of_gene == "protein-coding", !grepl("^LOC\\d+", Symbol)) %>%
  pull(Symbol)

tpm_filtered <- tpm_data %>%
  filter(GeneSymb %in% prot_coding)

cat("  TPM matrix:", nrow(tpm_filtered), "genes x", ncol(tpm_filtered) - 1, "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Prepare Expression Matrix
# -----------------------------------------------------------------------------
cat("Step 2: Preparing expression matrix...\n")
sample_annot <- readRDS(file.path(output_dir, "sample_annotation.rds"))

tpm_t <- tpm_filtered %>%
  column_to_rownames("GeneSymb") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleName")

tpm_annot <- inner_join(
  sample_annot[, c("SampleName", "SampleNameGroup")],
  tpm_t,
  by = "SampleName"
) %>%
  select(-SampleName) %>%
  column_to_rownames("SampleNameGroup") %>%
  t() %>%
  as.data.frame()

cat("  Expression matrix:", nrow(tpm_annot), "genes x", ncol(tpm_annot), "samples\n")

# -----------------------------------------------------------------------------
# Step 3: Run PROGENy
# -----------------------------------------------------------------------------
cat("Step 3: Running PROGENy...\n")

# BIOSTATISTICAL FIX: Random seed already set at script start
pathway_activity <- progeny(
  as.matrix(tpm_annot),
  scale = FALSE,
  organism = "Human",
  top = 100,
  perm = 1000
)

cat("  Pathway activity:", nrow(pathway_activity), "samples x", ncol(pathway_activity), "pathways\n")

# -----------------------------------------------------------------------------
# Step 4: Save Outputs
# -----------------------------------------------------------------------------
cat("Step 4: Saving outputs...\n")
saveRDS(pathway_activity, file.path(output_dir, "progeny_pathway_activity.rds"))

cat("\n=== PROGENy complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/01_human_bulk_rnaseq/03_run_progeny.R
git commit -m "feat(bulk-rnaseq): add PROGENy analysis script with random seed"
```

---

### Task 2.4: Extract Human Bulk RNA-seq Correlation Analysis with FDR Correction

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/04_correlations.R`

**Step 1: Write correlation analysis script with FDR correction**

```r
#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/04_correlations.R
# Calculate correlations between gene expression and pathway activity
# BIOSTATISTICAL FIX: Implements Benjamini-Hochberg FDR correction
#
# Inputs:
#   - analysis/01_human_bulk_rnaseq/outputs/progeny_pathway_activity.rds
#   - analysis/01_human_bulk_rnaseq/outputs/gsva_estrogen_pathways.rds
#   - data/human_bulk_rnaseq/raw/HumanERpAge_39404g168s_TPMlog2.txt
#
# Outputs:
#   - results/original/human_bulk_rnaseq/correlations_no_fdr.csv
#   - results/corrected/human_bulk_rnaseq/correlations_with_fdr.csv

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tibble)
})

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")
results_original <- file.path(project_root, "results/original/human_bulk_rnaseq")
results_corrected <- file.path(project_root, "results/corrected/human_bulk_rnaseq")

dir.create(results_original, showWarnings = FALSE, recursive = TRUE)
dir.create(results_corrected, showWarnings = FALSE, recursive = TRUE)

cat("=== Correlation Analysis with FDR Correction ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")
pathway_activity <- readRDS(file.path(output_dir, "progeny_pathway_activity.rds"))
gsva_result <- readRDS(file.path(output_dir, "gsva_estrogen_pathways.rds"))
sample_annot <- readRDS(file.path(output_dir, "sample_annotation.rds"))

# Load TPM data for gene expression
tpm_file <- file.path(data_dir, "raw/HumanERpAge_39404g168s_TPMlog2.txt")
tpm_data <- fread(tpm_file, stringsAsFactors = FALSE, header = TRUE)
colnames(tpm_data)[1] <- "GeneSymb"
colnames(tpm_data) <- gsub("_LEE.*", "", colnames(tpm_data))

# -----------------------------------------------------------------------------
# Step 2: Define Genes and Pathways of Interest
# -----------------------------------------------------------------------------
cat("Step 2: Setting up analysis...\n")

# Genes of interest (from original analysis)
genes_of_interest <- c("PAK4", "HSD17B7", "GREB1", "PGR", "ESR1",
                       "TFF1", "CYP19A1", "HSD17B2", "SAA1", "Age")

# Filter to tumor samples only, excluding middle age
tumor_samples <- sample_annot %>%
  filter(Group == "Tumor", AgeRange %in% c("Young", "Elderly")) %>%
  pull(SampleNameGroup)

cat("  Genes:", length(genes_of_interest) - 1, "+ Age\n")
cat("  Tumor samples (Young + Elderly):", length(tumor_samples), "\n")

# -----------------------------------------------------------------------------
# Step 3: Prepare Combined Data
# -----------------------------------------------------------------------------
cat("Step 3: Preparing combined data...\n")

# Get Estrogen pathway activity from PROGENy
progeny_estrogen <- pathway_activity[tumor_samples, "Estrogen", drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("SampleNameGroup")

# Get GSVA scores
gsva_t <- gsva_result %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleNameGroup") %>%
  filter(SampleNameGroup %in% tumor_samples)

# Get gene expression
tpm_t <- tpm_data %>%
  column_to_rownames("GeneSymb") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleName")

gene_expr <- inner_join(
  sample_annot[, c("SampleName", "SampleNameGroup", "Age")],
  tpm_t,
  by = "SampleName"
) %>%
  filter(SampleNameGroup %in% tumor_samples) %>%
  select(SampleNameGroup, Age, any_of(genes_of_interest[genes_of_interest != "Age"]))

# Combine all data
combined_data <- gene_expr %>%
  inner_join(progeny_estrogen, by = "SampleNameGroup") %>%
  inner_join(gsva_t, by = "SampleNameGroup")

pathway_names <- c("Estrogen", colnames(gsva_t)[-1])
cat("  Combined data:", nrow(combined_data), "samples\n")
cat("  Pathways:", length(pathway_names), "\n")

# -----------------------------------------------------------------------------
# Step 4: Calculate Correlations
# -----------------------------------------------------------------------------
cat("Step 4: Calculating correlations...\n")

correlation_results <- data.frame()

for (gene in genes_of_interest) {
  for (pathway in pathway_names) {

    if (gene == "Age") {
      x_values <- combined_data$Age
    } else {
      if (!gene %in% colnames(combined_data)) next
      x_values <- combined_data[[gene]]
    }

    y_values <- combined_data[[pathway]]

    # Skip if insufficient data
    if (sum(!is.na(x_values) & !is.na(y_values)) < 5) next

    # Spearman correlation
    cor_test <- cor.test(x_values, y_values, method = "spearman", exact = FALSE)

    correlation_results <- rbind(correlation_results, data.frame(
      GeneSymb = gene,
      PathwayName = pathway,
      Spearman_Rho = cor_test$estimate,
      Spearman_pval = cor_test$p.value,
      n_samples = sum(!is.na(x_values) & !is.na(y_values))
    ))
  }
}

cat("  Total correlation tests:", nrow(correlation_results), "\n")

# -----------------------------------------------------------------------------
# Step 5: Apply FDR Correction
# -----------------------------------------------------------------------------
cat("Step 5: Applying FDR correction...\n")

# BIOSTATISTICAL FIX: Benjamini-Hochberg FDR correction
correlation_results$FDR_qval <- p.adjust(correlation_results$Spearman_pval, method = "BH")

# Add significance flags
correlation_results$Sig_nominal <- correlation_results$Spearman_pval < 0.05
correlation_results$Sig_FDR <- correlation_results$FDR_qval < 0.05

# Summary statistics
n_sig_nominal <- sum(correlation_results$Sig_nominal)
n_sig_fdr <- sum(correlation_results$Sig_FDR)
cat("  Significant at p < 0.05:", n_sig_nominal, "\n")
cat("  Significant at FDR < 0.05:", n_sig_fdr, "\n")

# -----------------------------------------------------------------------------
# Step 6: Save Results
# -----------------------------------------------------------------------------
cat("Step 6: Saving results...\n")

# Original (no FDR) - for comparison
original_results <- correlation_results %>%
  select(GeneSymb, PathwayName, Spearman_Rho, Spearman_pval, n_samples, Sig_nominal)

# Corrected (with FDR)
corrected_results <- correlation_results %>%
  select(GeneSymb, PathwayName, Spearman_Rho, Spearman_pval, FDR_qval,
         n_samples, Sig_nominal, Sig_FDR)

write.csv(original_results, file.path(results_original, "correlations_no_fdr.csv"),
          row.names = FALSE)
write.csv(corrected_results, file.path(results_corrected, "correlations_with_fdr.csv"),
          row.names = FALSE)

# Also save to analysis outputs
saveRDS(correlation_results, file.path(output_dir, "correlation_results.rds"))

cat("\n=== Correlation analysis complete ===\n")
cat("Original results:", file.path(results_original, "correlations_no_fdr.csv"), "\n")
cat("Corrected results:", file.path(results_corrected, "correlations_with_fdr.csv"), "\n")
```

**Step 2: Commit**

```bash
git add analysis/01_human_bulk_rnaseq/04_correlations.R
git commit -m "feat(bulk-rnaseq): add correlation analysis with FDR correction

BIOSTATISTICAL FIX: Implements Benjamini-Hochberg FDR correction on all
88 correlation tests. Outputs both original (no FDR) and corrected results
for comparison."
```

---

### Task 2.5: Extract Human Bulk RNA-seq Visualization Script

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/05_visualize.R`

**Step 1: Write visualization script**

```r
#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/05_visualize.R
# Generate bubble plots and other visualizations
#
# Inputs:
#   - analysis/01_human_bulk_rnaseq/outputs/correlation_results.rds
#
# Outputs:
#   - figures/by_analysis/human_bulk_rnaseq/correlation_bubbleplot_original.pdf
#   - figures/by_analysis/human_bulk_rnaseq/correlation_bubbleplot_fdr.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
fig_dir <- file.path(project_root, "figures/by_analysis/human_bulk_rnaseq")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Generating Visualizations ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Correlation Results
# -----------------------------------------------------------------------------
cat("Step 1: Loading correlation results...\n")
corr_results <- readRDS(file.path(output_dir, "correlation_results.rds"))

# Filter out RAB19 and prepare for plotting
# Note: Document rationale - RAB19 was removed in original analysis
plot_data <- corr_results %>%
  filter(!grepl("RAB19", GeneSymb)) %>%
  mutate(
    GeneSymb = factor(GeneSymb, levels = c(
      "TFF1", "SAA1", "PGR", "PAK4", "HSD17B7", "HSD17B2",
      "GREB1", "ESR1", "CYP19A1", "Age"
    )),
    # Truncate long pathway names for display
    PathwayName_short = gsub("HALLMARK_", "HM_", PathwayName),
    PathwayName_short = gsub("REACTOME_", "RC_", PathwayName_short),
    PathwayName_short = gsub("GOBP_", "GO_", PathwayName_short)
  )

# Color palette
my_palette <- colorRampPalette(c("blue", "dodgerblue", "yellow", "orange", "red"))(100)

# -----------------------------------------------------------------------------
# Step 2: Original Bubble Plot (no FDR)
# -----------------------------------------------------------------------------
cat("Step 2: Creating original bubble plot...\n")

p_original <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName_short)) +
  geom_point(aes(size = -log10(Spearman_pval), color = Spearman_Rho)) +
  scale_color_gradientn("Spearman Rho", colors = my_palette, limits = c(-1, 1)) +
  scale_size_continuous("-log10(p)", range = c(1, 10)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black")
  ) +
  coord_flip() +
  ggtitle("Gene-Pathway Correlations (Original, no FDR)")

ggsave(
  file.path(fig_dir, "correlation_bubbleplot_original.pdf"),
  p_original, width = 12, height = 7
)

# -----------------------------------------------------------------------------
# Step 3: FDR-Corrected Bubble Plot
# -----------------------------------------------------------------------------
cat("Step 3: Creating FDR-corrected bubble plot...\n")

p_fdr <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName_short)) +
  geom_point(aes(size = -log10(FDR_qval), color = Spearman_Rho)) +
  scale_color_gradientn("Spearman Rho", colors = my_palette, limits = c(-1, 1)) +
  scale_size_continuous("-log10(FDR)", range = c(1, 10)) +
  # Add significance threshold line
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black")
  ) +
  coord_flip() +
  ggtitle("Gene-Pathway Correlations (FDR-Corrected)")

ggsave(
  file.path(fig_dir, "correlation_bubbleplot_fdr.pdf"),
  p_fdr, width = 12, height = 7
)

# -----------------------------------------------------------------------------
# Step 4: Side-by-side Comparison
# -----------------------------------------------------------------------------
cat("Step 4: Creating side-by-side comparison...\n")

library(patchwork)

p_combined <- p_original + p_fdr +
  plot_annotation(
    title = "Effect of FDR Correction on Gene-Pathway Correlations",
    subtitle = paste0(
      "Significant at p<0.05: ", sum(plot_data$Sig_nominal), " | ",
      "Significant at FDR<0.05: ", sum(plot_data$Sig_FDR)
    )
  )

ggsave(
  file.path(fig_dir, "correlation_bubbleplot_comparison.pdf"),
  p_combined, width = 20, height = 8
)

cat("\n=== Visualization complete ===\n")
cat("Figures saved to:", fig_dir, "\n")
```

**Step 2: Commit**

```bash
git add analysis/01_human_bulk_rnaseq/05_visualize.R
git commit -m "feat(bulk-rnaseq): add visualization script with original and FDR-corrected plots"
```

---

### Task 2.6: Create Human Bulk RNA-seq Run Script

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/run_analysis.sh`

**Step 1: Write run script**

```bash
#!/bin/bash
# analysis/01_human_bulk_rnaseq/run_analysis.sh
# Runs complete human bulk RNA-seq analysis pipeline
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Human Bulk RNA-seq Analysis Pipeline ==="
echo "Working directory: $SCRIPT_DIR"
echo ""

# Create outputs directory
mkdir -p outputs

# Run analysis scripts in order
echo "[1/5] Preprocessing..."
Rscript 01_preprocess.R

echo "[2/5] Running GSVA..."
Rscript 02_run_gsva.R

echo "[3/5] Running PROGENy..."
Rscript 03_run_progeny.R

echo "[4/5] Calculating correlations..."
Rscript 04_correlations.R

echo "[5/5] Generating visualizations..."
Rscript 05_visualize.R

echo ""
echo "=== Human Bulk RNA-seq Analysis Complete ==="
echo "Outputs: $SCRIPT_DIR/outputs/"
echo "Figures: ../../figures/by_analysis/human_bulk_rnaseq/"
```

**Step 2: Make executable and commit**

```bash
chmod +x analysis/01_human_bulk_rnaseq/run_analysis.sh
git add analysis/01_human_bulk_rnaseq/run_analysis.sh
git commit -m "feat(bulk-rnaseq): add master run script"
```

---

## Phase 2 Continued: Rat snRNA-seq Scripts

### Task 2.7: Create Rat snRNA-seq Preprocessing Script with Doublet Detection

**Files:**
- Create: `analysis/02_rat_snrnaseq/01_qc_filter.R`

**Step 1: Write QC and filtering script with doublet detection**

```r
#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/01_qc_filter.R
# QC, filtering, and doublet detection for rat snRNA-seq
# BIOSTATISTICAL FIX: Enable doublet detection (was commented out in original)
#
# Inputs:
#   - data/rat_snrnaseq/raw/Lee_021924_Nuclei*/
#
# Outputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_qc_filtered.rds
#   - analysis/02_rat_snrnaseq/outputs/qc_plots.pdf

# BIOSTATISTICAL FIX: Set random seed for reproducibility
set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(DoubletFinder)  # BIOSTATISTICAL FIX: Enable doublet detection
})

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
data_dir <- file.path(project_root, "data/rat_snrnaseq/raw")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Rat snRNA-seq QC and Filtering ===\n")
cat("Random seed: 12345\n\n")

# -----------------------------------------------------------------------------
# Step 1: Find Sample Directories
# -----------------------------------------------------------------------------
cat("Step 1: Finding sample directories...\n")
sample_dirs <- list.dirs(data_dir, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("Lee_021924_Nuclei", sample_dirs)]
cat("  Found", length(sample_dirs), "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Load and Process Each Sample
# -----------------------------------------------------------------------------
cat("Step 2: Loading and processing samples...\n")

seurat_list <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  cat("  Processing:", sample_name, "\n")

  # Read 10X data
  expr_matrix <- Read10X(data.dir = sample_dir, gene.column = 2)
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = sample_name)
  seurat_obj$orig.ident <- sample_name

  # Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^Mt-")
  seurat_obj[["percent.rbp"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]|^Rp[sl]")

  cat("    Pre-filter cells:", ncol(seurat_obj), "\n")

  # Apply QC filters
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
             nCount_RNA > 400 & percent.mt < 15
  )

  cat("    Post-filter cells:", ncol(seurat_obj), "\n")

  # SCTransform normalization (required for DoubletFinder)
  seurat_obj <- SCTransform(
    seurat_obj,
    method = "glmGamPoi",
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )

  # BIOSTATISTICAL FIX: Run DoubletFinder
  cat("    Running DoubletFinder...\n")

  # PCA required for DoubletFinder
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

  # Find optimal pK
  sweep_res <- paramSweep(seurat_obj, PCs = 1:30, sct = TRUE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  optimal_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # Estimate doublet rate (~5% for 10X)
  n_cells <- ncol(seurat_obj)
  doublet_rate <- 0.05  # Adjust based on loading density
  n_exp_doublets <- round(doublet_rate * n_cells)

  # Run DoubletFinder
  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs = 1:30,
    pN = 0.25,
    pK = optimal_pk,
    nExp = n_exp_doublets,
    sct = TRUE
  )

  # Get classification column name (varies by run)
  df_col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)[1]

  # Remove doublets
  seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[df_col]] == "Singlet"])
  cat("    Singlets retained:", ncol(seurat_obj), "\n")

  seurat_list[[sample_name]] <- seurat_obj
}

# -----------------------------------------------------------------------------
# Step 3: Merge Samples
# -----------------------------------------------------------------------------
cat("Step 3: Merging samples...\n")

seurat_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

cat("  Total cells after merge:", ncol(seurat_merged), "\n")

# -----------------------------------------------------------------------------
# Step 4: Generate QC Plots
# -----------------------------------------------------------------------------
cat("Step 4: Generating QC plots...\n")

pdf(file.path(output_dir, "qc_plots.pdf"), width = 12, height = 10)

# Violin plots of QC metrics
VlnPlot(seurat_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "orig.ident", pt.size = 0)

# Feature scatter
FeatureScatter(seurat_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               group.by = "orig.ident")

# Cells per sample
cells_per_sample <- data.frame(table(seurat_merged$orig.ident))
ggplot(cells_per_sample, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Sample", y = "Number of Cells", title = "Cells per Sample (Post-QC)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# -----------------------------------------------------------------------------
# Step 5: Save Output
# -----------------------------------------------------------------------------
cat("Step 5: Saving output...\n")
saveRDS(seurat_merged, file.path(output_dir, "seurat_qc_filtered.rds"))

cat("\n=== QC and filtering complete ===\n")
cat("Output:", file.path(output_dir, "seurat_qc_filtered.rds"), "\n")
```

**Step 2: Commit**

```bash
git add analysis/02_rat_snrnaseq/01_qc_filter.R
git commit -m "feat(snrnaseq): add QC script with doublet detection enabled

BIOSTATISTICAL FIX: DoubletFinder is now enabled (was commented out).
Random seed set for reproducibility."
```

---

### Task 2.8: Create Rat snRNA-seq Differential Expression Script

**Files:**
- Create: `analysis/02_rat_snrnaseq/08_differential_expression.R`

**Step 1: Write DE analysis script**

```r
#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/08_differential_expression.R
# Differential expression analysis: Young vs Aged per cell type
# BIOSTATISTICAL FIX: This analysis was missing from original
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - results/corrected/rat_snrnaseq/DE_results_by_celltype.csv
#   - figures/by_analysis/rat_snrnaseq/DE_volcano_plots.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
results_dir <- file.path(project_root, "results/corrected/rat_snrnaseq")
fig_dir <- file.path(project_root, "figures/by_analysis/rat_snrnaseq")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Differential Expression Analysis: Young vs Aged ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading annotated Seurat object...\n")
seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# Ensure AgeGroup is set correctly
table(seurat_obj$AgeGroup)

# Set RNA assay as default for DE
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)

# -----------------------------------------------------------------------------
# Step 2: Run DE for Each Cell Type
# -----------------------------------------------------------------------------
cat("Step 2: Running differential expression per cell type...\n")

cell_types <- unique(seurat_obj$CellTypeByMarker_RatsnRNAseq)
all_de_results <- data.frame()

for (ct in cell_types) {
  cat("  Processing:", ct, "\n")

  # Subset to cell type
  seurat_ct <- subset(seurat_obj, CellTypeByMarker_RatsnRNAseq == ct)

  # Check sample sizes
  n_aged <- sum(seurat_ct$AgeGroup == "Aged")
  n_young <- sum(seurat_ct$AgeGroup == "Young")
  cat("    Aged:", n_aged, "| Young:", n_young, "\n")

  # Skip if too few cells
  if (n_aged < 10 || n_young < 10) {
    cat("    Skipping - too few cells\n")
    next
  }

  # Set identity to AgeGroup
  Idents(seurat_ct) <- "AgeGroup"

  # Run FindMarkers (Wilcoxon test with BH correction)
  tryCatch({
    de_results <- FindMarkers(
      seurat_ct,
      ident.1 = "Aged",
      ident.2 = "Young",
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.25
    )

    if (nrow(de_results) > 0) {
      de_results$gene <- rownames(de_results)
      de_results$celltype <- ct
      de_results$FDR <- p.adjust(de_results$p_val, method = "BH")
      all_de_results <- rbind(all_de_results, de_results)
      cat("    DE genes:", nrow(de_results), "\n")
    }
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
  })
}

# -----------------------------------------------------------------------------
# Step 3: Summarize Results
# -----------------------------------------------------------------------------
cat("\nStep 3: Summarizing results...\n")

all_de_results <- all_de_results %>%
  arrange(FDR, p_val) %>%
  mutate(
    Sig_nominal = p_val < 0.05,
    Sig_FDR = FDR < 0.05,
    Direction = ifelse(avg_log2FC > 0, "Up_in_Aged", "Down_in_Aged")
  )

# Summary by cell type
de_summary <- all_de_results %>%
  group_by(celltype) %>%
  summarise(
    total_DE_genes = n(),
    sig_FDR = sum(Sig_FDR),
    up_in_aged = sum(Sig_FDR & Direction == "Up_in_Aged"),
    down_in_aged = sum(Sig_FDR & Direction == "Down_in_Aged")
  )

print(de_summary)

# -----------------------------------------------------------------------------
# Step 4: Generate Volcano Plots
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating volcano plots...\n")

pdf(file.path(fig_dir, "DE_volcano_plots.pdf"), width = 12, height = 10)

for (ct in unique(all_de_results$celltype)) {
  ct_results <- all_de_results %>% filter(celltype == ct)

  # Label top genes
  top_genes <- ct_results %>%
    filter(Sig_FDR) %>%
    slice_min(FDR, n = 10)

  p <- ggplot(ct_results, aes(x = avg_log2FC, y = -log10(FDR))) +
    geom_point(aes(color = Sig_FDR), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
    geom_text_repel(data = top_genes, aes(label = gene), max.overlaps = 20) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    theme_bw() +
    labs(
      title = paste("Differential Expression:", ct),
      subtitle = paste("Aged vs Young | FDR<0.05:", sum(ct_results$Sig_FDR)),
      x = "log2 Fold Change (Aged/Young)",
      y = "-log10(FDR)"
    ) +
    theme(legend.position = "none")

  print(p)
}

dev.off()

# -----------------------------------------------------------------------------
# Step 5: Save Results
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving results...\n")
write.csv(all_de_results, file.path(results_dir, "DE_results_by_celltype.csv"), row.names = FALSE)
write.csv(de_summary, file.path(results_dir, "DE_summary_by_celltype.csv"), row.names = FALSE)

cat("\n=== Differential expression analysis complete ===\n")
cat("Results:", file.path(results_dir, "DE_results_by_celltype.csv"), "\n")
```

**Step 2: Commit**

```bash
git add analysis/02_rat_snrnaseq/08_differential_expression.R
git commit -m "feat(snrnaseq): add differential expression analysis Young vs Aged

BIOSTATISTICAL FIX: Adds formal DE testing that was missing from original.
Uses Wilcoxon test with BH-FDR correction per cell type."
```

---

### Task 2.9: Create Rat snRNA-seq Differential Abundance Script

**Files:**
- Create: `analysis/02_rat_snrnaseq/09_differential_abundance.R`

**Step 1: Write DA analysis script**

```r
#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/09_differential_abundance.R
# Differential abundance analysis of cell types between Young vs Aged
# BIOSTATISTICAL FIX: This statistical test was missing from original
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - results/corrected/rat_snrnaseq/DA_results_celltypes.csv
#   - figures/by_analysis/rat_snrnaseq/DA_proportion_plots.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(speckle)  # For propeller differential abundance
  library(tidyr)
})

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
results_dir <- file.path(project_root, "results/corrected/rat_snrnaseq")
fig_dir <- file.path(project_root, "figures/by_analysis/rat_snrnaseq")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Differential Abundance Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")
seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# Get metadata
meta <- seurat_obj@meta.data %>%
  select(orig.ident, AgeGroup, CellTypeByMarker_RatsnRNAseq, CellTypeMacroTcell_RatsnRNAseq)

cat("  Total cells:", nrow(meta), "\n")
cat("  Samples:", length(unique(meta$orig.ident)), "\n")

# -----------------------------------------------------------------------------
# Step 2: Calculate Cell Type Proportions
# -----------------------------------------------------------------------------
cat("Step 2: Calculating cell type proportions...\n")

# Main cell types
prop_main <- meta %>%
  group_by(orig.ident, AgeGroup, CellTypeByMarker_RatsnRNAseq) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  )

# Detailed subtypes
prop_sub <- meta %>%
  group_by(orig.ident, AgeGroup, CellTypeMacroTcell_RatsnRNAseq) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  )

# -----------------------------------------------------------------------------
# Step 3: Run Propeller DA Test
# -----------------------------------------------------------------------------
cat("Step 3: Running propeller differential abundance test...\n")

# Prepare data for propeller
# Note: propeller requires counts, sample IDs, and cluster IDs

# Main cell types
da_main <- propeller(
  clusters = seurat_obj$CellTypeByMarker_RatsnRNAseq,
  sample = seurat_obj$orig.ident,
  group = seurat_obj$AgeGroup
)

da_main$CellType <- rownames(da_main)
da_main$FDR <- p.adjust(da_main$P.Value, method = "BH")
da_main$Level <- "Main"

cat("\nMain cell type results:\n")
print(da_main %>% select(CellType, PropMean.Aged, PropMean.Young, P.Value, FDR) %>% arrange(P.Value))

# Detailed subtypes (if enough cells per category)
tryCatch({
  da_sub <- propeller(
    clusters = seurat_obj$CellTypeMacroTcell_RatsnRNAseq,
    sample = seurat_obj$orig.ident,
    group = seurat_obj$AgeGroup
  )

  da_sub$CellType <- rownames(da_sub)
  da_sub$FDR <- p.adjust(da_sub$P.Value, method = "BH")
  da_sub$Level <- "Subtype"

  cat("\nSubtype results:\n")
  print(da_sub %>% select(CellType, PropMean.Aged, PropMean.Young, P.Value, FDR) %>% arrange(P.Value))

  da_all <- rbind(da_main, da_sub)
}, error = function(e) {
  cat("  Warning: Subtype analysis failed -", e$message, "\n")
  da_all <- da_main
})

# -----------------------------------------------------------------------------
# Step 4: Generate Proportion Plots
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating proportion plots...\n")

pdf(file.path(fig_dir, "DA_proportion_plots.pdf"), width = 14, height = 10)

# Stacked bar plot by sample
p1 <- ggplot(prop_main, aes(x = orig.ident, y = proportion, fill = CellTypeByMarker_RatsnRNAseq)) +
  geom_bar(stat = "identity") +
  facet_wrap(~AgeGroup, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Cell Type Proportions by Sample",
    x = "Sample", y = "Proportion", fill = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)

# Box plot comparing proportions
p2 <- ggplot(prop_main, aes(x = CellTypeByMarker_RatsnRNAseq, y = proportion * 100, fill = AgeGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.7) +
  theme_bw() +
  labs(
    title = "Cell Type Proportions: Young vs Aged",
    subtitle = "Note: n=3 per group - interpret with caution",
    x = "Cell Type", y = "Percentage", fill = "Age Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p2)

# DA results visualization
da_plot_data <- da_all %>%
  mutate(
    Sig = FDR < 0.05,
    LogFC = log2(PropMean.Aged / PropMean.Young)
  )

p3 <- ggplot(da_plot_data, aes(x = reorder(CellType, LogFC), y = LogFC, fill = Sig)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~Level, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  theme_bw() +
  labs(
    title = "Differential Abundance: Aged vs Young",
    subtitle = "Red = FDR < 0.05",
    x = "Cell Type", y = "log2 Fold Change (Aged/Young)"
  )
print(p3)

dev.off()

# -----------------------------------------------------------------------------
# Step 5: Save Results
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving results...\n")

# Add sample size warning to results
da_all$Warning <- "n=3 per group; interpret with caution"

write.csv(da_all, file.path(results_dir, "DA_results_celltypes.csv"), row.names = FALSE)
write.csv(prop_main, file.path(results_dir, "cell_proportions_by_sample.csv"), row.names = FALSE)

# Also document the limitation
writeLines(
  c(
    "# Differential Abundance Analysis Notes",
    "",
    "## Sample Size Limitation",
    "- Young: n=3 samples",
    "- Aged: n=3 samples",
    "",
    "With only 3 samples per group, statistical power is severely limited.",
    "Effect sizes and biological trends may be more meaningful than p-values.",
    "",
    "## Method",
    "- propeller (speckle package) for differential abundance",
    "- Benjamini-Hochberg FDR correction",
    ""
  ),
  file.path(results_dir, "DA_analysis_notes.md")
)

cat("\n=== Differential abundance analysis complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/02_rat_snrnaseq/09_differential_abundance.R
git commit -m "feat(snrnaseq): add differential abundance analysis

BIOSTATISTICAL FIX: Adds formal cell type proportion testing using propeller.
Documents n=3 per group limitation explicitly."
```

---

## Phase 3: Comparison Framework

### Task 3.1: Create Results Comparison Script

**Files:**
- Create: `scripts/compare_results.R`

**Step 1: Write comparison script**

```r
#!/usr/bin/env Rscript
# scripts/compare_results.R
# Compare original vs corrected results and generate comparison report
#
# Outputs:
#   - results/comparison/correlation_comparison.csv
#   - results/comparison/significant_findings_summary.md
#   - results/comparison/comparison_report.html
#   - results/comparison/comparison_report.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(rmarkdown)
  library(knitr)
})

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, ".."))
results_dir <- file.path(project_root, "results/comparison")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Generating Comparison Report ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Correlation Results
# -----------------------------------------------------------------------------
cat("Step 1: Loading correlation results...\n")

original_file <- file.path(project_root, "results/original/human_bulk_rnaseq/correlations_no_fdr.csv")
corrected_file <- file.path(project_root, "results/corrected/human_bulk_rnaseq/correlations_with_fdr.csv")

if (!file.exists(original_file) || !file.exists(corrected_file)) {
  cat("Warning: Correlation files not found. Run analyses first.\n")
  original <- NULL
  corrected <- NULL
} else {
  original <- read.csv(original_file)
  corrected <- read.csv(corrected_file)
}

# -----------------------------------------------------------------------------
# Step 2: Create Comparison Table
# -----------------------------------------------------------------------------
cat("Step 2: Creating comparison table...\n")

if (!is.null(corrected)) {
  comparison <- corrected %>%
    mutate(
      Status = case_when(
        Sig_nominal & Sig_FDR ~ "Remains Significant",
        Sig_nominal & !Sig_FDR ~ "Lost Significance (FDR)",
        !Sig_nominal ~ "Not Significant"
      )
    ) %>%
    select(GeneSymb, PathwayName, Spearman_Rho, Spearman_pval, FDR_qval, Status)

  # Summary statistics
  n_total <- nrow(comparison)
  n_sig_original <- sum(comparison$Spearman_pval < 0.05)
  n_sig_fdr <- sum(comparison$FDR_qval < 0.05)
  n_lost <- sum(comparison$Status == "Lost Significance (FDR)")

  cat("  Total tests:", n_total, "\n")
  cat("  Significant (p<0.05):", n_sig_original, "\n")
  cat("  Significant (FDR<0.05):", n_sig_fdr, "\n")
  cat("  Lost significance:", n_lost, "\n")

  write.csv(comparison, file.path(results_dir, "correlation_comparison.csv"), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# Step 3: Generate Summary Markdown
# -----------------------------------------------------------------------------
cat("Step 3: Generating summary markdown...\n")

summary_md <- c(
  "# Biostatistical Corrections: Results Comparison",
  "",
  "## Overview",
  "",
  "This document summarizes the impact of biostatistical corrections on the analysis results.",
  "",
  "## Human Bulk RNA-seq Correlations",
  "",
  if (!is.null(corrected)) {
    c(
      paste("- **Total correlation tests:**", n_total),
      paste("- **Significant at p < 0.05 (original):**", n_sig_original),
      paste("- **Significant at FDR < 0.05 (corrected):**", n_sig_fdr),
      paste("- **Lost significance after FDR:**", n_lost),
      "",
      "### Key Findings Status",
      ""
    )
  } else {
    "Results not yet available. Run the bulk RNA-seq analysis first."
  },
  "",
  "## Rat snRNA-seq Analysis",
  "",
  "### New Analyses Added",
  "- Differential expression (Young vs Aged) per cell type",
  "- Differential abundance testing using propeller",
  "- Doublet detection enabled",
  "",
  "### Sample Size Limitation",
  "- n=3 per age group",
  "- Statistical power is limited",
  "- Effect sizes may be more informative than p-values",
  "",
  "## Corrections Applied",
  "",
  "### Bulk RNA-seq",

  "1. DESeq2 design: `~1` -> `~ AgeRange + Group`",
  "2. Multiple testing: None -> Benjamini-Hochberg FDR",
  "3. Random seeds added for PROGENy",
  "4. QC plots added (PCA, library size, sample distances)",
  "",
  "### snRNA-seq",
  "1. DoubletFinder enabled (was commented out)",
  "2. Random seeds added for UMAP, Harmony, clustering",
  "3. Differential expression analysis added",
  "4. Differential abundance analysis added",
  "5. Cell type validation with reference datasets"
)

writeLines(summary_md, file.path(results_dir, "significant_findings_summary.md"))

# -----------------------------------------------------------------------------
# Step 4: Generate Comparison Figures
# -----------------------------------------------------------------------------
cat("Step 4: Generating comparison figures...\n")

fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

if (!is.null(corrected)) {
  # Bubble plot comparison
  pdf(file.path(fig_dir, "bubbleplot_sidebyside.pdf"), width = 16, height = 8)

  plot_data <- corrected %>%
    filter(!grepl("RAB19", GeneSymb)) %>%
    mutate(
      GeneSymb = factor(GeneSymb, levels = c(
        "TFF1", "SAA1", "PGR", "PAK4", "HSD17B7", "HSD17B2",
        "GREB1", "ESR1", "CYP19A1", "Age"
      ))
    )

  my_palette <- colorRampPalette(c("blue", "dodgerblue", "yellow", "orange", "red"))(100)

  p1 <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName)) +
    geom_point(aes(size = -log10(Spearman_pval), color = Spearman_Rho)) +
    scale_color_gradientn("Rho", colors = my_palette, limits = c(-1, 1)) +
    scale_size_continuous("-log10(p)", range = c(1, 8)) +
    coord_flip() + theme_bw() +
    ggtitle("Original (no FDR)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  p2 <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName)) +
    geom_point(aes(size = -log10(FDR_qval), color = Spearman_Rho)) +
    scale_color_gradientn("Rho", colors = my_palette, limits = c(-1, 1)) +
    scale_size_continuous("-log10(FDR)", range = c(1, 8)) +
    coord_flip() + theme_bw() +
    ggtitle("Corrected (with FDR)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  print(p1 + p2 + plot_annotation(title = "Effect of FDR Correction"))
  dev.off()

  # Significance change plot
  pdf(file.path(fig_dir, "significance_change_heatmap.pdf"), width = 10, height = 8)

  p3 <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName, fill = Status)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c(
      "Remains Significant" = "#2ECC71",
      "Lost Significance (FDR)" = "#E74C3C",
      "Not Significant" = "#BDC3C7"
    )) +
    coord_flip() + theme_bw() +
    ggtitle("Significance Status After FDR Correction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p3)
  dev.off()
}

# -----------------------------------------------------------------------------
# Step 5: Generate HTML and PDF Report
# -----------------------------------------------------------------------------
cat("Step 5: Generating HTML and PDF reports...\n")

# Create R Markdown content
rmd_content <- c(
  "---",
  'title: "Biostatistical Corrections: Comparison Report"',
  'author: "Automated Analysis"',
  paste0('date: "', Sys.Date(), '"'),
  "output:",
  "  html_document:",
  "    toc: true",
  "    toc_float: true",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "```",
  "",
  readLines(file.path(results_dir, "significant_findings_summary.md")),
  "",
  "## Figures",
  "",
  "### Correlation Bubble Plots",
  "",
  paste0("![Side-by-side comparison](figures/bubbleplot_sidebyside.pdf){width=100%}"),
  "",
  "### Significance Changes",
  "",
  paste0("![Significance status](figures/significance_change_heatmap.pdf){width=80%}")
)

rmd_file <- file.path(results_dir, "comparison_report.Rmd")
writeLines(rmd_content, rmd_file)

# Render to HTML
tryCatch({
  rmarkdown::render(
    rmd_file,
    output_format = "html_document",
    output_file = "comparison_report.html",
    quiet = TRUE
  )
  cat("  HTML report generated\n")
}, error = function(e) {
  cat("  Warning: HTML generation failed -", e$message, "\n")
})

# Convert to PDF using pagedown
tryCatch({
  pagedown::chrome_print(
    file.path(results_dir, "comparison_report.html"),
    output = file.path(results_dir, "comparison_report.pdf")
  )
  cat("  PDF report generated\n")
}, error = function(e) {
  cat("  Warning: PDF generation failed -", e$message, "\n")
  cat("  You may need to install chromium: conda install -c conda-forge chromium\n")
})

cat("\n=== Comparison report complete ===\n")
cat("Outputs in:", results_dir, "\n")
```

**Step 2: Commit**

```bash
chmod +x scripts/compare_results.R
git add scripts/compare_results.R
git commit -m "feat: add results comparison script with HTML/PDF report generation"
```

---

### Task 3.2: Create Master Run Script

**Files:**
- Create: `scripts/run_all.sh`

**Step 1: Write master run script**

```bash
#!/bin/bash
# scripts/run_all.sh
# Master script to run complete analysis pipeline
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "=============================================="
echo "  ERpBRCA OlderWomen Analysis Pipeline"
echo "=============================================="
echo "Project root: $PROJECT_ROOT"
echo "Started at: $(date)"
echo ""

# Parse arguments
RUN_DOWNLOAD=true
RUN_BULK=true
RUN_SNRNA=true
RUN_WES=true
RUN_COMPARE=true

while [[ $# -gt 0 ]]; do
  case $1 in
    --skip-download) RUN_DOWNLOAD=false; shift ;;
    --bulk-only) RUN_SNRNA=false; RUN_WES=false; shift ;;
    --snrna-only) RUN_BULK=false; RUN_WES=false; shift ;;
    --skip-compare) RUN_COMPARE=false; shift ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Step 1: Download data
if [ "$RUN_DOWNLOAD" = true ]; then
  echo ""
  echo "[1/5] Downloading data from GEO..."
  echo "----------------------------------------------"
  bash data/download_all.sh
else
  echo "[1/5] Skipping data download..."
fi

# Step 2: Human bulk RNA-seq
if [ "$RUN_BULK" = true ]; then
  echo ""
  echo "[2/5] Running human bulk RNA-seq analysis..."
  echo "----------------------------------------------"
  cd analysis/01_human_bulk_rnaseq
  bash run_analysis.sh
  cd "$PROJECT_ROOT"
else
  echo "[2/5] Skipping bulk RNA-seq analysis..."
fi

# Step 3: Rat snRNA-seq
if [ "$RUN_SNRNA" = true ]; then
  echo ""
  echo "[3/5] Running rat snRNA-seq analysis..."
  echo "----------------------------------------------"
  cd analysis/02_rat_snrnaseq
  bash run_analysis.sh
  cd "$PROJECT_ROOT"
else
  echo "[3/5] Skipping snRNA-seq analysis..."
fi

# Step 4: Rat WES
if [ "$RUN_WES" = true ]; then
  echo ""
  echo "[4/5] Running rat WES analysis..."
  echo "----------------------------------------------"
  cd analysis/03_rat_wes
  bash run_analysis.sh
  cd "$PROJECT_ROOT"
else
  echo "[4/5] Skipping WES analysis..."
fi

# Step 5: Generate comparison report
if [ "$RUN_COMPARE" = true ]; then
  echo ""
  echo "[5/5] Generating comparison report..."
  echo "----------------------------------------------"
  Rscript scripts/compare_results.R
else
  echo "[5/5] Skipping comparison report..."
fi

# Step 6: Organize figures
echo ""
echo "[6/6] Organizing figures for manuscript..."
echo "----------------------------------------------"
Rscript scripts/organize_figures.R 2>/dev/null || echo "Figure organization script not yet implemented"

echo ""
echo "=============================================="
echo "  Pipeline Complete"
echo "=============================================="
echo "Finished at: $(date)"
echo ""
echo "Results:"
echo "  - Comparison report: results/comparison/comparison_report.pdf"
echo "  - Figures: figures/manuscript/"
echo "  - Corrected results: results/corrected/"
echo ""
```

**Step 2: Make executable and commit**

```bash
chmod +x scripts/run_all.sh
git add scripts/run_all.sh
git commit -m "feat: add master run script for complete pipeline"
```

---

### Task 3.3: Update Main README

**Files:**
- Modify: `README.md`

**Step 1: Write updated README**

```markdown
# Chronic Inflammation and Hormone Disposition in Breast Cancer of Older Women

[![DOI](https://img.shields.io/badge/DOI-10.1101%2FXXXX-blue)](https://doi.org/10.1101/XXXX)
[![GEO](https://img.shields.io/badge/GEO-GSE276755-green)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276755)

**Carleton et al. 2024**

> Chronic inflammation and hormone disposition promote a tumor-permissive locale for breast cancer in older women

## Overview

This repository contains all analysis code for the study examining how the aged breast tumor microenvironment promotes ER+ breast cancer development through chronic inflammation and estrogen signaling.

### Key Findings

- Aged F344 rats develop tumors faster than younger counterparts
- Tumors upregulate HSD17B7 to convert E1 to E2 locally
- Tumor-associated macrophages integrate E2 and chemokine signaling
- Age-related immune dysfunction associated with chronic inflammation

## Repository Structure

```
├── analysis/
│   ├── 01_human_bulk_rnaseq/    # Human ER+ breast cancer RNA-seq
│   ├── 02_rat_snrnaseq/         # Rat mammary tumor snRNA-seq
│   └── 03_rat_wes/              # Rat whole exome sequencing
├── data/                         # Data download scripts (GEO)
├── figures/                      # Generated figures
├── results/                      # Analysis outputs
├── scripts/                      # Pipeline scripts
└── docs/                         # Documentation
```

## Quick Start

### 1. Setup Environment

```bash
# Clone repository
git clone https://github.com/leeoesterreich/ERpBRCA_OlderWomen.git
cd ERpBRCA_OlderWomen

# Create conda environment
conda env create -f environment.yml
conda activate erp_brca_aging
```

### 2. Download Data

```bash
# Download all data from GEO
bash data/download_all.sh
```

### 3. Run Analysis

```bash
# Run complete pipeline
bash scripts/run_all.sh

# Or run individual analyses
bash analysis/01_human_bulk_rnaseq/run_analysis.sh
bash analysis/02_rat_snrnaseq/run_analysis.sh
bash analysis/03_rat_wes/run_analysis.sh
```

## Data Availability

| Dataset | GEO Accession | Description |
|---------|---------------|-------------|
| Human Bulk RNA-seq | [GSE276755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276755) | 168 samples |
| Rat snRNA-seq | [GSE276758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276758) | 6 samples |
| Rat Bulk RNA-seq | [GSE276757](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276757) | - |
| Rat WES | [GSE276759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276759) | - |

## Methods Summary

### Human Bulk RNA-seq
- GSVA for estrogen pathway scoring
- PROGENy for pathway activity inference
- Spearman correlations with FDR correction

### Rat snRNA-seq
- Seurat v5 + Harmony integration
- DoubletFinder for doublet removal
- Differential expression (Wilcoxon + BH-FDR)
- Differential abundance (propeller)

## Citation

```bibtex
@article{carleton2024chronic,
  title={Chronic inflammation and hormone disposition promote a tumor-permissive locale for breast cancer in older women},
  author={Carleton, Neil and others},
  journal={bioRxiv},
  year={2024},
  doi={10.1101/XXXX}
}
```

## License

This project is licensed under the terms of the [LICENSE](LICENSE) file.

## Contact

- Neil Carleton - [email]
- Lee-Oesterreich Lab - University of Pittsburgh
```

**Step 2: Commit**

```bash
git add README.md
git commit -m "docs: update README with unified repository structure and instructions"
```

---

## Phase 4-5: Remaining Tasks (Summary)

### Task 4.1-4.5: Complete snRNA-seq Scripts

Extract remaining scripts from README:
- `02_normalize_integrate.R`
- `03_cluster_annotate.R`
- `04_subset_myeloid.R`
- `05_subset_nkt.R`
- `06_cell_fractions.R`
- `07_visualize.R`
- `10_statistical_summary.R`
- `run_analysis.sh`

### Task 4.6-4.8: WES Analysis Scripts

- Clean up existing notebook
- Add random seeds
- Create run script

### Task 5.1: Figure Organization Script

- Create `scripts/organize_figures.R`
- Map analysis outputs to manuscript figures

### Task 5.2: Final Documentation

- Write `docs/methods.md` with detailed methods
- Create figure legends document

### Task 5.3: Validation Run

- Run complete pipeline
- Verify all figures generate
- Confirm comparison report is correct

---

## Success Criteria Checklist

- [ ] All analyses run from `./scripts/run_all.sh`
- [ ] Data downloads automatically from GEO
- [ ] Random seeds ensure identical results
- [ ] FDR-corrected results alongside original
- [ ] Comparison report shows what changed
- [ ] All manuscript figures in `figures/manuscript/`
- [ ] Conda environment installs without errors
- [ ] README provides clear setup instructions
