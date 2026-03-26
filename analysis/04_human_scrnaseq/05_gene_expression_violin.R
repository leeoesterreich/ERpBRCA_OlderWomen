#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/05_gene_expression_violin.R
# Gene expression violin plots by age group and cell type
#
# Matches Sanghoon's original code which generates violin plots using:
#   - slot='data' (log-normalized) for CellPhoneDB genes (Step 7)
#   - slot='scale.data' (z-scored) for macrophage markers CCL2/TGFB1/CD163/MRC1 (Step 8)
#   - Both AgeGroup and per-patient (CaseID sorted by age) groupings
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#   - analysis/04_human_scrnaseq/outputs/clinical_data.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/violin_*.pdf
#   - analysis/04_human_scrnaseq/figures/violin_*.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(data.table)
})

# -----------------------------------------------------------------------------
# Workaround for Seurat 5.x / patchwork S4 deparse bug
# https://github.com/satijalab/seurat/issues/7653
# -----------------------------------------------------------------------------
safe_vlnplot <- function(seurat_obj, features, group.by, slot = "data",
                         pt.size = 0, title = NULL) {
  plots <- lapply(features, function(gene) {
    tryCatch({
      VlnPlot(seurat_obj, features = gene, group.by = group.by,
              slot = slot, pt.size = pt.size) +
        theme(legend.position = "none")
    }, error = function(e) {
      message("  Warning: VlnPlot failed for ", gene, ": ", e$message)
      NULL
    })
  })

  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(NULL)

  combined <- wrap_plots(plots, ncol = min(3, length(plots)))
  if (!is.null(title)) {
    combined <- combined + plot_annotation(title = title)
  }
  return(combined)
}

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
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Gene Expression Violin Plots ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load and prepare data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# NormalizeData + ScaleData should already be done in 00b_preprocess_seurat.R
# Verify and run if missing
DefaultAssay(seurat_obj) <- "RNA"

# Check if data slot is populated
if (all(GetAssayData(seurat_obj, slot = "data")[1:5, 1:5] ==
        GetAssayData(seurat_obj, slot = "counts")[1:5, 1:5])) {
  cat("  Running NormalizeData (not yet done)...\n")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
}

# Check if scale.data exists — use tryCatch for Seurat 5 compatibility
has_scale_data <- tryCatch({
  sd <- GetAssayData(seurat_obj, slot = "scale.data")
  nrow(sd) > 0
}, error = function(e) FALSE)

if (has_scale_data) {
  cat("  scale.data present\n")
} else {
  cat("  Running ScaleData (all genes)...\n")
  all_genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all_genes, verbose = FALSE)
}

# Order age groups
seurat_obj$AgeGroup <- factor(seurat_obj$AgeGroup, levels = c("Young", "MidAge", "Elderly"))

# Load clinical data for patient ordering by age
clinical_file <- file.path(output_dir, "clinical_data.rds")
if (file.exists(clinical_file)) {
  clinical_data <- readRDS(clinical_file)
  # Sort CaseIDs by age (ascending)
  clinical_sorted <- clinical_data %>%
    filter(CaseID %in% unique(seurat_obj$orig.ident)) %>%
    arrange(Age)
  case_id_order <- clinical_sorted$CaseID
  seurat_obj$CaseID <- factor(seurat_obj$orig.ident, levels = case_id_order)
  cat("  Patient order by age:", paste(case_id_order, collapse = ", "), "\n")
} else {
  cat("  WARNING: clinical_data.rds not found, skipping per-patient plots\n")
  case_id_order <- NULL
}

# Set identity to cell type
Idents(seurat_obj) <- seurat_obj$CellTypeAnnot

cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Define gene sets (matching original code)
# -----------------------------------------------------------------------------
cat("\nStep 2: Defining gene sets...\n")

# Immune checkpoint (original Step 5)
immune_checkpoint <- c("CTLA4", "PDCD1", "PDCD1LG2", "CD274", "LAG3", "HAVCR2")

# M2 macrophage markers (original Step 8 — uses scale.data)
macrophage_markers <- c("CCL2", "CCL3", "CCL4", "TNF", "TGFB1", "CD163", "MRC1")

# Cytokines/chemokines (original Step 7 Feature_Mac)
cytokines <- c("IL1A", "IL1B", "IL1RN", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7",
               "CXCL8", "IL9", "IL10", "IL12B", "IL13", "IL15", "IL17A", "IL33",
               "CD40LG", "EGF", "CCL11", "FGF2", "CSF3", "CSF2", "IFNA17", "IFNG",
               "CXCL1", "CXCL2", "CXCL3")

# CellPhoneDB genes (original Step 7 MyFeature)
cellphonedb_genes <- c("CCR1", "CCR3", "CCR4", "CD74", "APP", "CXCL10", "DPP4",
                       "EGFR", "AGER", "GPR75", "CCL5", "IL2RG", "IL15", "S100A11")

# Hedgehog pathway (original Step 4)
hedgehog <- c("SHH", "KLF4", "GLI1", "KIF7", "SMO")

filter_genes <- function(genes) {
  present <- genes[genes %in% rownames(seurat_obj)]
  if (length(present) < length(genes)) {
    missing <- setdiff(genes, present)
    cat("  Missing:", paste(missing, collapse = ", "), "\n")
  }
  present
}

immune_checkpoint <- filter_genes(immune_checkpoint)
macrophage_markers <- filter_genes(macrophage_markers)
cytokines <- filter_genes(cytokines)
cellphonedb_genes <- filter_genes(cellphonedb_genes)
hedgehog <- filter_genes(hedgehog)

# -----------------------------------------------------------------------------
# Step 3: Generate violin plots by cell type
# Original generates plots for ALL 16 cell types
# -----------------------------------------------------------------------------
cat("\nStep 3: Generating plots per cell type...\n")

cell_types <- sort(unique(seurat_obj$CellTypeAnnot))
cat("  Cell types:", length(cell_types), "\n")

for (ct in cell_types) {
  cat("  Processing:", ct, "\n")

  seurat_subset <- subset(seurat_obj, CellTypeAnnot == ct)

  if (ncol(seurat_subset) < 10) {
    cat("    Skipping (too few cells:", ncol(seurat_subset), ")\n")
    next
  }

  # --- Immune checkpoint (slot='data', group by AgeGroup) ---
  if (length(immune_checkpoint) > 0) {
    p <- safe_vlnplot(seurat_subset, immune_checkpoint, "AgeGroup",
                      slot = "data", title = paste(ct, "- Immune Checkpoint"))
    if (!is.null(p)) {
      pdf(file.path(output_dir, paste0("violin_checkpoint_", ct, ".pdf")), width = 12, height = 8)
      print(p)
      dev.off()
      ggsave(file.path(figures_dir, paste0("violin_checkpoint_", ct, ".png")),
             p, width = 12, height = 8, dpi = 300, bg = "white")
    }
  }

  # --- M2 macrophage markers (slot='scale.data', group by AgeGroup) ---
  # Original Step 8: VlnPlot(..., slot = 'scale.data')
  if (length(macrophage_markers) > 0) {
    p <- safe_vlnplot(seurat_subset, macrophage_markers, "AgeGroup",
                      slot = "scale.data",
                      title = paste(ct, "- M2 Markers (scaled)"))
    if (!is.null(p)) {
      pdf(file.path(output_dir, paste0("violin_m2markers_scaled_", ct, ".pdf")),
          width = 12, height = 8)
      print(p)
      dev.off()
      ggsave(file.path(figures_dir, paste0("violin_m2markers_scaled_", ct, ".png")),
             p, width = 12, height = 8, dpi = 300, bg = "white")
    }

    # Also per-patient (slot='scale.data', group by CaseID sorted by age)
    if (!is.null(case_id_order)) {
      p2 <- safe_vlnplot(seurat_subset, macrophage_markers, "CaseID",
                         slot = "scale.data",
                         title = paste(ct, "- M2 Markers by Patient (scaled)"))
      if (!is.null(p2)) {
        pdf(file.path(output_dir, paste0("violin_m2markers_scaled_bypatient_", ct, ".pdf")),
            width = 14, height = 8)
        print(p2)
        dev.off()
        ggsave(file.path(figures_dir, paste0("violin_m2markers_scaled_bypatient_", ct, ".png")),
               p2, width = 14, height = 8, dpi = 300, bg = "white")
      }
    }
  }

  # --- CellPhoneDB genes (slot='data', group by AgeGroup) ---
  if (length(cellphonedb_genes) > 0) {
    p <- safe_vlnplot(seurat_subset, cellphonedb_genes, "AgeGroup",
                      slot = "data",
                      title = paste(ct, "- CellPhoneDB Genes"))
    if (!is.null(p)) {
      pdf(file.path(output_dir, paste0("violin_cellphonedb_", ct, ".pdf")),
          width = 12, height = 11)
      print(p)
      dev.off()
      ggsave(file.path(figures_dir, paste0("violin_cellphonedb_", ct, ".png")),
             p, width = 12, height = 11, dpi = 300, bg = "white")
    }
  }

  # --- Cytokines (slot='data', group by AgeGroup) ---
  if (length(cytokines) > 0) {
    p <- safe_vlnplot(seurat_subset, cytokines, "AgeGroup",
                      slot = "data",
                      title = paste(ct, "- Cytokines/Chemokines"))
    if (!is.null(p)) {
      pdf(file.path(output_dir, paste0("violin_cytokines_", ct, ".pdf")),
          width = 15, height = 15)
      print(p)
      dev.off()
      ggsave(file.path(figures_dir, paste0("violin_cytokines_", ct, ".png")),
             p, width = 15, height = 15, dpi = 300, bg = "white")
    }
  }
}

# --- Hedgehog pathway: all cell types together (matching original Step 4) ---
cat("\n  Hedgehog pathway (all cell types)...\n")
if (length(hedgehog) > 0) {
  p <- safe_vlnplot(seurat_obj, hedgehog, "CellTypeAnnot",
                    slot = "data",
                    title = "Hedgehog Pathway by Cell Type")
  if (!is.null(p)) {
    pdf(file.path(output_dir, "violin_hedgehog_bycelltype.pdf"), width = 12, height = 9)
    print(p)
    dev.off()
    ggsave(file.path(figures_dir, "violin_hedgehog_bycelltype.png"),
           p, width = 12, height = 9, dpi = 300, bg = "white")
  }
}

cat("\n=== Violin plots complete ===\n")
