#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/01_load_subset_data.R
# Load and subset Wu et al. scRNA-seq data for age group analysis
#
# Inputs:
#   - data/human_scrnaseq/SeuratObj_GSE176078_ERpos_AfterQCSCT.rds
#   - data/human_scrnaseq/ClinicalData_Wu_scRNAseq_26p.txt
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_young_midage_elderly.rds

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Define paths
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
data_dir <- file.path(project_root, "data/human_scrnaseq")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Load and Subset scRNA-seq Data ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load clinical data
# -----------------------------------------------------------------------------
cat("Step 1: Loading clinical data...\n")

clinical_file <- file.path(data_dir, "ClinicalData_Wu_scRNAseq_26p.txt")
clinical_data <- fread(clinical_file, header = TRUE, stringsAsFactors = FALSE)
colnames(clinical_data) <- gsub(" ", "", colnames(clinical_data))

# Define age groups
clinical_data$AgeGroup <- ifelse(
  clinical_data$Age > 80, "Elderly",
  ifelse(clinical_data$Age <= 50, "Young", "MidAge")
)

cat("  Total patients:", nrow(clinical_data), "\n")
cat("  Age groups:\n")
print(table(clinical_data$AgeGroup))

# -----------------------------------------------------------------------------
# Step 2: Load Seurat object
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading Seurat object...\n")

seurat_file <- file.path(data_dir, "SeuratObj_GSE176078_ERpos_AfterQCSCT.rds")
seurat_obj <- readRDS(seurat_file)

cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

# Note: PCA and UMAP already computed in 00b_preprocess_seurat.R

# -----------------------------------------------------------------------------
# Step 3: Subset to Young/MidAge/Elderly
# -----------------------------------------------------------------------------
cat("\nStep 3: Subsetting to age groups...\n")

# Sample IDs for each group (from Wu et al. clinical data)
young_ids <- c("CID3941", "CID4530N", "CID4535")
midage_ids <- c("CID4463", "CID4040", "CID4471", "CID4461")
elderly_ids <- c("CID3948", "CID4067", "CID4290A")
all_ids <- c(young_ids, midage_ids, elderly_ids)

seurat_subset <- subset(seurat_obj, subset = orig.ident %in% all_ids)

# Add age group labels
seurat_subset$AgeGroup <- case_when(
  seurat_subset$orig.ident %in% young_ids ~ "Young",
  seurat_subset$orig.ident %in% midage_ids ~ "MidAge",
  seurat_subset$orig.ident %in% elderly_ids ~ "Elderly",
  TRUE ~ NA_character_
)

cat("  Subset cells:", ncol(seurat_subset), "\n")
cat("  By age group:\n")
print(table(seurat_subset$AgeGroup))

# -----------------------------------------------------------------------------
# Step 4: Save output
# -----------------------------------------------------------------------------
cat("\nStep 4: Saving output...\n")

saveRDS(seurat_subset, file.path(output_dir, "seurat_young_midage_elderly.rds"))
saveRDS(clinical_data, file.path(output_dir, "clinical_data.rds"))

cat("\n=== Load complete ===\n")
