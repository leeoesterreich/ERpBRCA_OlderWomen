#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/01_load_subset_data.R
# Load Xu et al. 2024 HR+ scRNA-seq data and assign age groups
#
# Age groups: Young (≤50), MidAge (51-80), Elderly (>80)
# Clinical data is extracted directly from Seurat metadata (no external file needed)
#
# Inputs:
#   - data/human_scrnaseq/SeuratObj_Xu2024_HRpos_AfterQCSCT.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_young_midage_elderly.rds
#   - analysis/04_human_scrnaseq/outputs/clinical_data.rds

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
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

cat("=== Load and Subset scRNA-seq Data (Xu et al. 2024) ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Seurat object
# -----------------------------------------------------------------------------
cat("Step 1: Loading Seurat object...\n")

seurat_file <- file.path(data_dir, "SeuratObj_Xu2024_HRpos_AfterQCSCT.rds")
seurat_obj <- readRDS(seurat_file)

# Fix gene names: pre-built RDS may have tab-embedded names ("ESR1\tESR1")
# from genes.tsv where both columns are identical. Extract just the symbol.
gene_names <- rownames(seurat_obj)
if (any(grepl("\t", gene_names))) {
  cat("  WARNING: Tab characters found in gene names — cleaning...\n")
  clean_names <- sapply(strsplit(gene_names, "\t"), function(x) x[length(x)])
  clean_names <- make.unique(clean_names)
  # Rename features in all assays
  for (assay_name in Assays(seurat_obj)) {
    assay_obj <- seurat_obj[[assay_name]]
    rownames(assay_obj) <- clean_names[match(rownames(assay_obj), gene_names)]
    seurat_obj[[assay_name]] <- assay_obj
  }
  cat("  Cleaned", sum(grepl("\t", gene_names)), "gene names\n")
}

cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

# Drop patient 0319 (ER-low/PR+, not ER+ by standard clinical definition)
drop_patients <- c("0319")
drop_cells <- sum(seurat_obj$Patient_ID %in% drop_patients)
if (drop_cells > 0) {
  cat("  Excluding patient(s):", paste(drop_patients, collapse = ", "),
      "(", drop_cells, "cells) — ER-low/PR+ only\n")
  seurat_obj <- subset(seurat_obj,
    cells = colnames(seurat_obj)[!seurat_obj$Patient_ID %in% drop_patients])
  cat("  Cells after exclusion:", ncol(seurat_obj), "\n")
}

# -----------------------------------------------------------------------------
# Step 2: Extract clinical data from Seurat metadata
# -----------------------------------------------------------------------------
cat("\nStep 2: Extracting clinical data...\n")

# Get one row per patient
meta <- seurat_obj@meta.data
patient_meta <- meta %>%
  group_by(Patient_ID) %>%
  summarise(
    Age_raw = first(Age),
    Dataset = first(Dataset),
    Clinical_Subtype = first(Clinical_Subtype),
    Molecular_Subtype = first(Molecular_Subtype),
    Treatment_Status = first(Treatment_Status),
    n_cells = n(),
    .groups = "drop"
  )

# Parse age — handle mixed formats: numeric ("58") and ranges ("60-65")
parse_age <- function(age_str) {
  age_str <- trimws(age_str)
  if (is.na(age_str) || age_str == "" || age_str == "NA") return(NA_real_)
  # Check for range format "XX-YY" → midpoint
  if (grepl("^\\d+-\\d+$", age_str)) {
    parts <- as.numeric(strsplit(age_str, "-")[[1]])
    return(mean(parts))
  }
  # Plain numeric
  val <- suppressWarnings(as.numeric(age_str))
  return(val)
}

patient_meta$Age <- sapply(patient_meta$Age_raw, parse_age)

# CaseID alias — downstream scripts (04, 05) join on CaseID = orig.ident
patient_meta$CaseID <- patient_meta$Patient_ID

cat("  Patients with parsed ages:\n")
print(patient_meta %>% select(Patient_ID, CaseID, Age_raw, Age, Dataset, n_cells))

# Check for unparseable ages
na_ages <- patient_meta %>% filter(is.na(Age))
if (nrow(na_ages) > 0) {
  cat("\n  WARNING: Patients with unparseable ages (will be excluded):\n")
  print(na_ages %>% select(Patient_ID, Age_raw, Dataset))
}

# -----------------------------------------------------------------------------
# Step 3: Assign age groups
# -----------------------------------------------------------------------------
cat("\nStep 3: Assigning age groups...\n")

patient_meta$AgeGroup <- case_when(
  is.na(patient_meta$Age) ~ NA_character_,
  patient_meta$Age <= 50   ~ "Young",
  patient_meta$Age > 80    ~ "Elderly",
  TRUE                     ~ "MidAge"
)

cat("  Age group distribution (patients):\n")
print(table(patient_meta$AgeGroup, useNA = "ifany"))

# Map age group back to all cells (unname to avoid Seurat 5 barcode mismatch)
age_map <- setNames(patient_meta$AgeGroup, patient_meta$Patient_ID)
seurat_obj$AgeGroup <- unname(age_map[seurat_obj$Patient_ID])

# Map numeric age to all cells
age_numeric_map <- setNames(patient_meta$Age, patient_meta$Patient_ID)
seurat_obj$Age <- unname(age_numeric_map[seurat_obj$Patient_ID])

# Remove cells with NA age group (if any)
na_cells <- sum(is.na(seurat_obj$AgeGroup))
if (na_cells > 0) {
  cat("  Removing", na_cells, "cells with missing age\n")
  seurat_obj <- subset(seurat_obj, subset = AgeGroup %in% c("Young", "MidAge", "Elderly"))
}

cat("  Age group distribution (cells):\n")
print(table(seurat_obj$AgeGroup))

cat("  By age group and dataset:\n")
print(table(seurat_obj$AgeGroup, seurat_obj$Dataset))

# -----------------------------------------------------------------------------
# Step 4: Save output
# -----------------------------------------------------------------------------
cat("\nStep 4: Saving output...\n")

saveRDS(seurat_obj, file.path(output_dir, "seurat_young_midage_elderly.rds"))
saveRDS(patient_meta, file.path(output_dir, "clinical_data.rds"))

cat("  Saved seurat_young_midage_elderly.rds:", ncol(seurat_obj), "cells\n")
cat("  Saved clinical_data.rds:", nrow(patient_meta), "patients\n")

cat("\n=== Load complete ===\n")
