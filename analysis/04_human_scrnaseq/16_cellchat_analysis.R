#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/16_cellchat_analysis.R
# Cell-cell communication analysis using CellChat (R alternative to CellPhoneDB)
# Generates Figure 7 Panel D - macrophage communication with T cells
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/macrophage_seurat.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/cellchat_elderly.rds
#   - analysis/04_human_scrnaseq/outputs/cellchat_young.rds
#   - analysis/04_human_scrnaseq/figures/fig7d_cellchat_dotplot.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

# Load CellChat (pre-installed in erp_brca_aging env)
library(CellChat)

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

cat("=== CellChat Analysis (Figure 7D) ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_file <- file.path(output_dir, "macrophage_seurat.rds")
seurat_obj <- readRDS(seurat_file)

cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Cell types:\n")
print(table(seurat_obj$CellTypeAnnotSH))

# -----------------------------------------------------------------------------
# Step 2: Function to run CellChat
# -----------------------------------------------------------------------------
run_cellchat <- function(seurat_subset, name) {
  cat(sprintf("\n=== Running CellChat for %s (%d cells) ===\n", name, ncol(seurat_subset)))

  # Prepare input data
  data.input <- GetAssayData(seurat_subset, slot = "data")
  labels <- seurat_subset$CellTypeAnnotSH
  meta <- data.frame(labels = labels, row.names = colnames(seurat_subset))

  # Create CellChat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

  # Set database (use human Secreted Signaling)
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use

  # Preprocessing
  cat("  Identifying overexpressed genes...\n")
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # Compute communication probability
  cat("  Computing communication probability...\n")
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  # Infer cell-cell communication at signaling pathway level
  cat("  Inferring signaling pathways...\n")
  cellchat <- computeCommunProbPathway(cellchat)

  # Aggregate
  cellchat <- aggregateNet(cellchat)

  cat("  Done!\n")
  return(cellchat)
}

# -----------------------------------------------------------------------------
# Step 3: Run CellChat for Elderly and Young
# -----------------------------------------------------------------------------
cat("\nStep 3: Running CellChat analysis...\n")

seurat_elderly <- subset(seurat_obj, subset = AgeGroup == "Elderly")
seurat_young <- subset(seurat_obj, subset = AgeGroup == "Young")

cellchat_elderly <- run_cellchat(seurat_elderly, "Elderly")
cellchat_young <- run_cellchat(seurat_young, "Young")

# Save CellChat objects
saveRDS(cellchat_elderly, file.path(output_dir, "cellchat_elderly.rds"))
saveRDS(cellchat_young, file.path(output_dir, "cellchat_young.rds"))
cat("\n  Saved CellChat objects\n")

# -----------------------------------------------------------------------------
# Step 4: Compare macrophage interactions
# -----------------------------------------------------------------------------
cat("\nStep 4: Comparing macrophage interactions...\n")

# Get interaction counts involving Macrophage
get_macrophage_interactions <- function(cellchat, name) {
  if (!"Macrophage" %in% levels(cellchat@idents)) {
    cat("  Warning: Macrophage not found in", name, "\n")
    return(NULL)
  }

  # Get network
  net <- cellchat@net$count

  # Get interactions where macrophage is source or target
  macro_as_source <- net["Macrophage", ]
  macro_as_target <- net[, "Macrophage"]

  # Focus on T cell interactions
  tcell_types <- grep("Tcells|NK", rownames(net), value = TRUE)

  data.frame(
    target = c(names(macro_as_source), names(macro_as_target)),
    count = c(macro_as_source, macro_as_target),
    direction = c(rep("outgoing", length(macro_as_source)),
                  rep("incoming", length(macro_as_target))),
    group = name
  )
}

macro_elderly <- get_macrophage_interactions(cellchat_elderly, "Elderly")
macro_young <- get_macrophage_interactions(cellchat_young, "Young")

if (!is.null(macro_elderly) && !is.null(macro_young)) {
  macro_compare <- bind_rows(macro_elderly, macro_young)

  # Focus on T cell related cell types
  macro_compare <- macro_compare %>%
    filter(grepl("Tcells|NK", target))

  cat("  Macrophage-T cell interactions:\n")
  print(macro_compare %>%
          group_by(target, group) %>%
          summarize(total = sum(count), .groups = "drop") %>%
          tidyr::pivot_wider(names_from = group, values_from = total))
}

# -----------------------------------------------------------------------------
# Step 5: Generate Figure 7D - Dot Plot
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating Figure 7D...\n")

# Create merged CellChat object for comparison
object.list <- list(Young = cellchat_young, Elderly = cellchat_elderly)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# Identify macrophage and T/NK cell types present in data
all_celltypes <- union(levels(cellchat_young@idents), levels(cellchat_elderly@idents))
macro_types <- grep("Macro|Mono", all_celltypes, value = TRUE)
tcell_types <- grep("Tcells|NK", all_celltypes, value = TRUE)

cat("  Source cell types (Macrophage/Monocyte):", paste(macro_types, collapse = ", "), "\n")
cat("  Target cell types (T/NK cells):", paste(tcell_types, collapse = ", "), "\n")

if (length(macro_types) > 0 && length(tcell_types) > 0) {
  # Generate dot plot - capture ggplot object and save explicitly
  p <- netVisual_bubble(cellchat_merged,
                        sources.use = macro_types,
                        targets.use = tcell_types,
                        comparison = c(1, 2),
                        angle.x = 45,
                        remove.isolate = TRUE,
                        title.name = "Macrophage/Monocyte -> T/NK Cell Communication",
                        return.data = FALSE)

  # Save using ggsave for reliable output
  ggsave(file.path(figures_dir, "fig7d_cellchat_dotplot.png"),
         plot = p, width = 14, height = 12, dpi = 300)
  cat("  Saved fig7d_cellchat_dotplot.png\n")
} else {
  cat("  WARNING: Could not find matching cell types for dot plot\n")
}

cat("\n=== CellChat analysis complete ===\n")
