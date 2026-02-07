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
  library(future)
})

# Increase memory limit for parallelization (4GB per worker)
options(future.globals.maxSize = 4 * 1024^3)

# Helper function to check file existence
check_file_exists <- function(filepath, description = "file") {
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: %s not found: %s", description, filepath))
  }
  cat(sprintf("  Found: %s\n", basename(filepath)))
}

# Helper function to check directory existence
check_dir_exists <- function(dirpath, description = "directory") {
  if (!dir.exists(dirpath)) {
    stop(sprintf("ERROR: %s not found: %s", description, dirpath))
  }
  cat(sprintf("  Found directory: %s\n", basename(dirpath)))
}

# Define paths - use commandArgs to get script directory when run via Rscript
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
data_dir <- file.path(project_root, "data/rat_snrnaseq/raw")

cat("=== Rat snRNA-seq QC and Filtering ===\n")
cat("Random seed: 12345\n")
cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Find Sample Directories
# -----------------------------------------------------------------------------
cat("Step 1: Finding sample directories...\n")
check_dir_exists(data_dir, "Raw data directory")

# Create output directory only after verifying input data exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sample_dirs <- list.dirs(data_dir, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("Lee_021924_Nuclei", sample_dirs)]

if (length(sample_dirs) == 0) {
  stop("ERROR: No sample directories matching 'Lee_021924_Nuclei*' found in ", data_dir)
}

cat("  Found", length(sample_dirs), "samples:\n")
for (sd in sample_dirs) {
  cat("    -", basename(sd), "\n")
}

# -----------------------------------------------------------------------------
# Step 2: Load and Process Each Sample
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading and processing samples...\n")

seurat_list <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  cat("  Processing:", sample_name, "\n")

  # Read 10X data with error handling
  expr_matrix <- tryCatch({
    Read10X(data.dir = sample_dir, gene.column = 2)
  }, error = function(e) {
    stop(sprintf("ERROR: Failed to read 10X data from %s: %s", sample_dir, e$message))
  })
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
cat("\nStep 3: Merging samples...\n")

seurat_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

cat("  Total cells after merge:", ncol(seurat_merged), "\n")

# -----------------------------------------------------------------------------
# Step 4: Save Output (before plots to prevent data loss)
# -----------------------------------------------------------------------------
cat("\nStep 4: Saving output...\n")
saveRDS(seurat_merged, file.path(output_dir, "seurat_qc_filtered.rds"))
cat("  Saved:", file.path(output_dir, "seurat_qc_filtered.rds"), "\n")

# -----------------------------------------------------------------------------
# Step 5: Generate QC Plots
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating QC plots...\n")

tryCatch({
  pdf(file.path(output_dir, "qc_plots.pdf"), width = 12, height = 10)

  # Violin plots of QC metrics - generate individually to avoid patchwork issues
  for (feat in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
    print(VlnPlot(seurat_merged, features = feat, group.by = "orig.ident", pt.size = 0))
  }

  # Feature scatter
  print(FeatureScatter(seurat_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                 group.by = "orig.ident"))

  # Cells per sample
  cells_per_sample <- data.frame(table(seurat_merged$orig.ident))
  print(ggplot(cells_per_sample, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    labs(x = "Sample", y = "Number of Cells", title = "Cells per Sample (Post-QC)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

  dev.off()
  cat("  QC plots saved to:", file.path(output_dir, "qc_plots.pdf"), "\n")
}, error = function(e) {
  cat("  Warning: QC plot generation failed:", conditionMessage(e), "\n")
  cat("  Data was saved successfully - plots can be regenerated\n")
  try(dev.off(), silent = TRUE)
})

# Memory cleanup
rm(seurat_list, expr_matrix)
gc()

cat("\n=== QC and filtering complete ===\n")
cat("Output:", file.path(output_dir, "seurat_qc_filtered.rds"), "\n")
