#!/usr/bin/env Rscript
# 15b_pathway_polarity_diagnostic.R
# Diagnostic script to investigate GSVA polarity reversal in Figure 7B/C
#
# Key questions:
# 1. What are the raw GSVA scores BEFORE z-normalization?
# 2. Does SCTransform vs CPM+log2 affect the direction?
# 3. Which cell types drive the pathway enrichment?
#
# Output: Diagnostic report and comparison plots

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(msigdbr)
  library(GSVA)
  library(ggplot2)
  library(tidyr)
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

cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("GSVA POLARITY DIAGNOSTIC - Figure 7B/C\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Load Seurat object
# -----------------------------------------------------------------------------
cat("Step 1: Loading Seurat object...\n")

seurat_file <- file.path(output_dir, "macrophage_seurat.rds")
seurat_obj <- readRDS(seurat_file)

# Filter to Elderly and Young only
seurat_ey <- subset(seurat_obj, subset = AgeGroup %in% c("Elderly", "Young"))
cat("  Cells:", ncol(seurat_ey), "\n")
cat("  Age groups:", table(seurat_ey$AgeGroup), "\n")
cat("  Cell types:", length(unique(seurat_ey$CellTypeAnnotSH)), "\n\n")

# -----------------------------------------------------------------------------
# Step 2: Check available assays and normalization
# -----------------------------------------------------------------------------
cat("Step 2: Checking available assays...\n")

assays_available <- Assays(seurat_ey)
cat("  Available assays:", paste(assays_available, collapse = ", "), "\n")

# Check if SCT assay exists
has_sct <- "SCT" %in% assays_available
cat("  Has SCTransform assay:", has_sct, "\n")

# Check current default assay
cat("  Default assay:", DefaultAssay(seurat_ey), "\n\n")

# -----------------------------------------------------------------------------
# Step 3: Load target pathways (curated from manuscript)
# -----------------------------------------------------------------------------
cat("Step 3: Loading target pathways...\n")

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

# Key pathways to check polarity
target_pathways <- c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

target_list <- hallmark_list[target_pathways]
cat("  Target pathways:", length(target_pathways), "\n\n")

# -----------------------------------------------------------------------------
# Step 4: Create pseudo-bulk profiles with BOTH normalization methods
# -----------------------------------------------------------------------------
cat("Step 4: Creating pseudo-bulk profiles...\n")

# Function to create pseudo-bulk
create_pseudobulk <- function(seurat_obj, assay_name, slot_name) {
  DefaultAssay(seurat_obj) <- assay_name
  counts <- GetAssayData(seurat_obj, layer = slot_name)

  # Create grouping
  seurat_obj$ct_age <- paste0(seurat_obj$CellTypeAnnotSH, "_", seurat_obj$AgeGroup)
  groups <- unique(seurat_obj$ct_age)

  # Aggregate
  pseudobulk <- sapply(groups, function(grp) {
    cells <- colnames(seurat_obj)[seurat_obj$ct_age == grp]
    if (length(cells) > 10) {
      Matrix::rowSums(counts[, cells, drop = FALSE])
    } else {
      rep(NA, nrow(counts))
    }
  })

  # Remove NA columns and rows with all zeros
  pseudobulk <- pseudobulk[, !apply(pseudobulk, 2, function(x) all(is.na(x)))]
  pseudobulk <- pseudobulk[rowSums(pseudobulk, na.rm = TRUE) > 0, ]

  return(pseudobulk)
}

# Method 1: CPM + log2 (current implementation)
cat("  Method 1: CPM + log2 normalization...\n")
pb_counts <- create_pseudobulk(seurat_ey, "RNA", "counts")
pb_cpm <- sweep(pb_counts, 2, colSums(pb_counts, na.rm = TRUE), "/") * 1e6
pb_log2_cpm <- log2(pb_cpm + 1)
cat("    Matrix:", nrow(pb_log2_cpm), "genes x", ncol(pb_log2_cpm), "groups\n")

# Method 2: SCT data if available
if (has_sct) {
  cat("  Method 2: SCTransform data...\n")
  pb_sct <- create_pseudobulk(seurat_ey, "SCT", "data")
  cat("    Matrix:", nrow(pb_sct), "genes x", ncol(pb_sct), "groups\n")
} else {
  cat("  Method 2: SCT not available, skipping\n")
  pb_sct <- NULL
}

cat("\n")

# -----------------------------------------------------------------------------
# Step 5: Run GSVA with BOTH methods (no z-score normalization)
# -----------------------------------------------------------------------------
cat("Step 5: Running GSVA (raw scores, no z-normalization)...\n")

run_gsva_raw <- function(expr_mat, pathways, method_name) {
  cat(sprintf("  Running GSVA for %s...\n", method_name))

  gsva_result <- gsva(
    gsvaParam(
      as.matrix(expr_mat),
      pathways,
      kcdf = "Gaussian",
      maxDiff = TRUE
    )
  )

  return(gsva_result)
}

# Run GSVA with CPM+log2
gsva_cpm <- run_gsva_raw(pb_log2_cpm, target_list, "CPM+log2")

# Run GSVA with SCT if available
if (!is.null(pb_sct)) {
  gsva_sct <- run_gsva_raw(pb_sct, target_list, "SCT")
} else {
  gsva_sct <- NULL
}

cat("\n")

# -----------------------------------------------------------------------------
# Step 6: Analyze polarity - Younger vs Older for each pathway
# -----------------------------------------------------------------------------
cat("Step 6: Analyzing polarity (Younger vs Older)...\n\n")

analyze_polarity <- function(gsva_mat, method_name) {
  cat(sprintf("--- %s Method ---\n", method_name))

  results <- data.frame()

  for (pathway in rownames(gsva_mat)) {
    # Get Younger and Older samples
    younger_cols <- grep("_Young$", colnames(gsva_mat))
    older_cols <- grep("_Elderly$", colnames(gsva_mat))

    younger_mean <- mean(gsva_mat[pathway, younger_cols], na.rm = TRUE)
    older_mean <- mean(gsva_mat[pathway, older_cols], na.rm = TRUE)

    diff <- older_mean - younger_mean
    direction <- ifelse(diff > 0, "Older > Younger", "Younger > Older")

    results <- rbind(results, data.frame(
      pathway = gsub("HALLMARK_", "", pathway),
      younger_mean = round(younger_mean, 4),
      older_mean = round(older_mean, 4),
      diff = round(diff, 4),
      direction = direction
    ))

    cat(sprintf("  %s:\n", gsub("HALLMARK_", "", pathway)))
    cat(sprintf("    Younger mean: %.4f\n", younger_mean))
    cat(sprintf("    Older mean:   %.4f\n", older_mean))
    cat(sprintf("    Direction:    %s (diff=%.4f)\n\n", direction, diff))
  }

  return(results)
}

polarity_cpm <- analyze_polarity(gsva_cpm, "CPM+log2")

if (!is.null(gsva_sct)) {
  polarity_sct <- analyze_polarity(gsva_sct, "SCT")

  # Compare directions
  cat("\n--- Direction Comparison (CPM vs SCT) ---\n")
  comparison <- merge(
    polarity_cpm[, c("pathway", "direction")],
    polarity_sct[, c("pathway", "direction")],
    by = "pathway",
    suffixes = c("_CPM", "_SCT")
  )
  comparison$match <- comparison$direction_CPM == comparison$direction_SCT
  print(comparison)

  cat(sprintf("\nMatching directions: %d/%d\n",
              sum(comparison$match), nrow(comparison)))
}

# -----------------------------------------------------------------------------
# Step 7: Z-score normalization effect
# -----------------------------------------------------------------------------
cat("\n" %>% rep(70) %>% paste(collapse = "-"), "\n")
cat("Step 7: Effect of Z-score normalization...\n\n")

gsva_cpm_z <- t(scale(t(gsva_cpm)))

cat("Before z-score (CPM+log2):\n")
cat(sprintf("  Range: %.4f to %.4f\n", min(gsva_cpm), max(gsva_cpm)))
cat(sprintf("  Mean:  %.4f\n", mean(gsva_cpm)))

cat("\nAfter z-score:\n")
cat(sprintf("  Range: %.4f to %.4f\n", min(gsva_cpm_z), max(gsva_cpm_z)))
cat(sprintf("  Mean:  %.4f\n", mean(gsva_cpm_z)))

# Check if z-score changes polarity
cat("\nPolarity after z-score normalization:\n")
polarity_z <- analyze_polarity(gsva_cpm_z, "CPM+log2 (z-scored)")

# -----------------------------------------------------------------------------
# Step 8: Cell type breakdown
# -----------------------------------------------------------------------------
cat("\n" %>% rep(70) %>% paste(collapse = "-"), "\n")
cat("Step 8: Cell type breakdown for key pathways...\n\n")

# Focus on ESTROGEN_RESPONSE_EARLY as example
pathway_focus <- "HALLMARK_ESTROGEN_RESPONSE_EARLY"
cat(sprintf("Pathway: %s\n\n", pathway_focus))

# Get per-cell-type scores
celltype_scores <- data.frame(
  group = colnames(gsva_cpm),
  score = gsva_cpm[pathway_focus, ]
) %>%
  mutate(
    celltype = gsub("_(Young|Elderly)$", "", group),
    age = ifelse(grepl("_Young$", group), "Younger", "Older")
  )

# Print per-cell-type
cat("Cell type breakdown:\n")
celltype_wide <- celltype_scores %>%
  select(celltype, age, score) %>%
  pivot_wider(names_from = age, values_from = score) %>%
  mutate(diff = Older - Younger) %>%
  arrange(desc(abs(diff)))

print(as.data.frame(celltype_wide))

# -----------------------------------------------------------------------------
# Step 9: Save diagnostic report
# -----------------------------------------------------------------------------
cat("\n" %>% rep(70) %>% paste(collapse = "-"), "\n")
cat("Step 9: Saving diagnostic outputs...\n")

# Save polarity comparison
write.csv(polarity_cpm, file.path(output_dir, "diagnostic_polarity_cpm.csv"), row.names = FALSE)
cat("  Saved diagnostic_polarity_cpm.csv\n")

if (!is.null(gsva_sct)) {
  write.csv(polarity_sct, file.path(output_dir, "diagnostic_polarity_sct.csv"), row.names = FALSE)
  cat("  Saved diagnostic_polarity_sct.csv\n")
}

# Save raw GSVA scores
write.csv(as.data.frame(gsva_cpm), file.path(output_dir, "diagnostic_gsva_raw_cpm.csv"))
cat("  Saved diagnostic_gsva_raw_cpm.csv\n")

# Create visualization
cat("  Creating diagnostic plot...\n")

plot_data <- polarity_cpm %>%
  mutate(pathway = factor(pathway, levels = pathway[order(diff)]))

p <- ggplot(plot_data, aes(x = pathway, y = diff, fill = diff > 0)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    labels = c("TRUE" = "Older > Younger", "FALSE" = "Younger > Older"),
                    name = "Direction") +
  labs(
    title = "GSVA Polarity: Older vs Younger (Raw scores, CPM+log2)",
    subtitle = "Positive = higher in Older, Negative = higher in Younger",
    x = "Pathway",
    y = "Mean difference (Older - Younger)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "diagnostic_gsva_polarity.png"), p,
       width = 10, height = 6, dpi = 300)
cat("  Saved diagnostic_gsva_polarity.png\n")

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")

# -----------------------------------------------------------------------------
# Summary: Key findings
# -----------------------------------------------------------------------------
cat("\n")
cat("KEY FINDINGS:\n")
cat("-" %>% rep(70) %>% paste(collapse = ""), "\n")

cat("\n1. NORMALIZATION METHOD:\n")
cat("   Current script uses: CPM + log2\n")
cat("   Manuscript states:   SCTransform per sample\n")
if (has_sct) {
  cat("   SCT assay available: YES\n")
} else {
  cat("   SCT assay available: NO (may need to re-process)\n")
}

cat("\n2. POLARITY DIRECTION (raw GSVA, no z-score):\n")
for (i in 1:nrow(polarity_cpm)) {
  cat(sprintf("   %s: %s\n", polarity_cpm$pathway[i], polarity_cpm$direction[i]))
}

cat("\n3. RECOMMENDATIONS:\n")
cat("   - If SCT assay exists, use SCT data instead of CPM+log2\n")
cat("   - Compare with manuscript Figure 7B/C to confirm expected direction\n")
cat("   - Consider using raw GSVA without z-normalization to preserve absolute direction\n")
cat("\n")
