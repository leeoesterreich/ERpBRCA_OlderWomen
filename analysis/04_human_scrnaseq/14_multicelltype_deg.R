#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/14_multicelltype_deg.R
# Pseudobulk DEG analysis across ALL cell types (Elderly vs Young)
# Aggregates cells per patient → DESeq2 at the patient level (no pseudoreplication)
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/multicelltype_deg_summary.csv
#   - analysis/04_human_scrnaseq/outputs/multicelltype_deg_full/  (per-celltype DEG tables)
#   - analysis/04_human_scrnaseq/figures/fig7b_celltype_deg_barplot.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(DESeq2)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "outputs")
figures_dir <- file.path(script_dir, "figures")
deg_dir <- file.path(output_dir, "multicelltype_deg_full")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Multi-Cell-Type Pseudobulk DEG Analysis (Figure 7B) ===\n")

# Minimum patients per group and minimum total cells per patient-celltype
MIN_PATIENTS_PER_GROUP <- 3
MIN_CELLS_PER_PSEUDOBULK <- 10

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_file <- file.path(output_dir, "seurat_annotated.rds")
seurat_obj <- readRDS(seurat_file)

cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Cell types:\n")
print(table(seurat_obj$CellTypeAnnotSH))

# Filter to Elderly and Young only
seurat_ey <- subset(seurat_obj, subset = AgeGroup %in% c("Elderly", "Young"))
cat("\n  Elderly + Young cells:", ncol(seurat_ey), "\n")
cat("  Patients: Young =", length(unique(seurat_ey$orig.ident[seurat_ey$AgeGroup == "Young"])),
    ", Elderly =", length(unique(seurat_ey$orig.ident[seurat_ey$AgeGroup == "Elderly"])), "\n")

rm(seurat_obj); gc()

# -----------------------------------------------------------------------------
# Step 2: Pseudobulk DEG per cell type
# -----------------------------------------------------------------------------
cat("\nStep 2: Running pseudobulk DEG analysis per cell type...\n")

cell_types <- sort(unique(seurat_ey$CellTypeAnnotSH))
cell_types <- cell_types[!is.na(cell_types)]

DefaultAssay(seurat_ey) <- "RNA"

deg_summary <- lapply(cell_types, function(ct) {
  cat("  Processing:", ct, "...")

  ct_cells <- subset(seurat_ey, subset = CellTypeAnnotSH == ct)

  # Count cells per patient
  cells_per_patient <- table(ct_cells$orig.ident)
  # Keep only patients with enough cells
  valid_patients <- names(cells_per_patient[cells_per_patient >= MIN_CELLS_PER_PSEUDOBULK])

  if (length(valid_patients) == 0) {
    cat(" skipped (no patients with >=", MIN_CELLS_PER_PSEUDOBULK, "cells)\n")
    return(data.frame(
      cell_type = ct, n_elderly_patients = 0, n_young_patients = 0,
      n_elderly_cells = 0, n_young_cells = 0,
      n_deg_up = NA, n_deg_down = NA, n_deg_total = NA
    ))
  }

  ct_cells <- subset(ct_cells, cells = colnames(ct_cells)[ct_cells$orig.ident %in% valid_patients])

  # Check minimum patients per age group (use data.frame without rownames for dplyr)
  patient_age <- data.frame(
    orig.ident = ct_cells$orig.ident,
    AgeGroup = ct_cells$AgeGroup,
    stringsAsFactors = FALSE
  ) %>% distinct()
  n_elderly_patients <- sum(patient_age$AgeGroup == "Elderly")
  n_young_patients <- sum(patient_age$AgeGroup == "Young")

  if (n_elderly_patients < MIN_PATIENTS_PER_GROUP || n_young_patients < MIN_PATIENTS_PER_GROUP) {
    cat(" skipped (", n_young_patients, "Young /", n_elderly_patients, "Elderly patients)\n")
    return(data.frame(
      cell_type = ct,
      n_elderly_patients = n_elderly_patients,
      n_young_patients = n_young_patients,
      n_elderly_cells = sum(ct_cells$AgeGroup == "Elderly"),
      n_young_cells = sum(ct_cells$AgeGroup == "Young"),
      n_deg_up = NA, n_deg_down = NA, n_deg_total = NA
    ))
  }

  tryCatch({
    # Pseudobulk: aggregate raw counts per patient
    pseudo <- AggregateExpression(
      ct_cells,
      assays = "RNA",
      return.seurat = FALSE,
      group.by = "orig.ident"
    )$RNA

    # Seurat mangles column names: prepends "g" to IDs starting with digits,
    # replaces "_" with "-". Reverse this to recover original Patient_IDs.
    orig_ids <- unique(ct_cells$orig.ident)
    seurat_mangle <- function(x) {
      x <- gsub("_", "-", x)
      ifelse(grepl("^[0-9]", x), paste0("g", x), x)
    }
    pseudo_to_orig <- setNames(orig_ids, seurat_mangle(orig_ids))
    colnames(pseudo) <- pseudo_to_orig[colnames(pseudo)]
    pseudo <- pseudo[, !is.na(colnames(pseudo)), drop = FALSE]

    # Build coldata
    coldata <- patient_age %>%
      filter(orig.ident %in% colnames(pseudo))
    rownames(coldata) <- coldata$orig.ident
    coldata$orig.ident <- NULL
    coldata$AgeGroup <- factor(coldata$AgeGroup, levels = c("Young", "Elderly"))
    pseudo <- pseudo[, rownames(coldata)]

    # Filter low-count genes (at least 10 counts in at least 3 patients)
    keep <- rowSums(pseudo >= 10) >= 3
    pseudo <- pseudo[keep, ]

    # DESeq2
    dds <- DESeqDataSetFromMatrix(
      countData = round(pseudo),  # ensure integer counts
      colData = coldata,
      design = ~ AgeGroup
    )
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("AgeGroup", "Elderly", "Young"),
                   alpha = 0.05)
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      arrange(padj)

    n_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
    n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)

    cat(" done (", n_up, "up,", n_down, "down |",
        n_young_patients, "Y /", n_elderly_patients, "E patients)\n")

    # Save full DEG table
    fwrite(res_df, file.path(deg_dir, paste0(gsub(" ", "_", ct), "_degs.csv")))

    data.frame(
      cell_type = ct,
      n_elderly_patients = n_elderly_patients,
      n_young_patients = n_young_patients,
      n_elderly_cells = sum(ct_cells$AgeGroup == "Elderly"),
      n_young_cells = sum(ct_cells$AgeGroup == "Young"),
      n_deg_up = n_up,
      n_deg_down = n_down,
      n_deg_total = n_up + n_down
    )
  }, error = function(e) {
    cat(" error:", conditionMessage(e), "\n")
    data.frame(
      cell_type = ct,
      n_elderly_patients = n_elderly_patients,
      n_young_patients = n_young_patients,
      n_elderly_cells = sum(ct_cells$AgeGroup == "Elderly"),
      n_young_cells = sum(ct_cells$AgeGroup == "Young"),
      n_deg_up = NA, n_deg_down = NA, n_deg_total = NA
    )
  })
})

deg_summary_df <- bind_rows(deg_summary)
deg_summary_df <- deg_summary_df %>% arrange(desc(n_deg_total))

cat("\n  DEG Summary:\n")
print(deg_summary_df)

fwrite(deg_summary_df, file.path(output_dir, "multicelltype_deg_summary.csv"))
cat("\n  Saved multicelltype_deg_summary.csv\n")

# -----------------------------------------------------------------------------
# Step 3: Generate Figure 7B - Bar plot
# -----------------------------------------------------------------------------
cat("\nStep 3: Generating Figure 7B...\n")

plot_df <- deg_summary_df %>%
  filter(!is.na(n_deg_total), n_deg_total > 0)

if (nrow(plot_df) > 0) {
  plot_long <- plot_df %>%
    select(cell_type, n_deg_up, n_deg_down) %>%
    tidyr::pivot_longer(cols = c(n_deg_up, n_deg_down),
                        names_to = "direction",
                        values_to = "n_degs") %>%
    mutate(direction = ifelse(direction == "n_deg_up", "Up in Elderly", "Down in Elderly"))

  p <- ggplot(plot_long, aes(x = reorder(cell_type, -n_degs), y = n_degs, fill = direction)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Up in Elderly" = "#E41A1C", "Down in Elderly" = "#377EB8")) +
    coord_flip() +
    labs(
      title = "Pseudobulk DEGs by Cell Type (Elderly vs Young)",
      subtitle = paste0("DESeq2 on patient-level pseudobulk (n=",
                        max(deg_summary_df$n_young_patients), " Young, ",
                        max(deg_summary_df$n_elderly_patients), " Elderly)"),
      x = "",
      y = "Number of DEGs (FDR < 0.05)",
      fill = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 10)
    )

  ggsave(file.path(figures_dir, "fig7b_celltype_deg_barplot.png"),
         p, width = 8, height = 8, dpi = 300)
  cat("  Saved fig7b_celltype_deg_barplot.png\n")
} else {
  cat("  No cell types with significant DEGs for plot\n")
}

cat("\n=== Multi-cell-type pseudobulk DEG analysis complete ===\n")
