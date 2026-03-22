#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/12_macrophage_pathway_enrichment.R
# Pathway enrichment on macrophage DEGs
# Generates Figure 7 Panel C visualization (HALLMARK/BIOCARTA pathways)
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/macrophage_degs.csv
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/macrophage_pathway_enrichment.csv
#   - analysis/04_human_scrnaseq/figures/fig7c_pathway_heatmap.png

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(msigdbr)
  library(tibble)
})

# msigdbr v25 uses collection/subcollection; v10 used category/subcategory.
# This wrapper tries the new API first and falls back to the old one.
safe_msigdbr <- function(species, coll, subcoll = NULL) {
  tryCatch({
    if (is.null(subcoll)) msigdbr(species = species, collection = coll)
    else msigdbr(species = species, collection = coll, subcollection = subcoll)
  }, error = function(e) {
    if (is.null(subcoll)) msigdbr(species = species, category = coll)
    else msigdbr(species = species, category = coll, subcategory = subcoll)
  })
}

# Simple enrichment function using hypergeometric test (Fisher's exact)
run_simple_enrichment <- function(deg_genes, pathway_list, all_genes) {
  results <- lapply(names(pathway_list), function(pw_name) {
    pw_genes <- pathway_list[[pw_name]]
    # Overlap
    overlap <- intersect(deg_genes, pw_genes)
    n_overlap <- length(overlap)
    n_deg <- length(deg_genes)
    n_pw <- length(intersect(pw_genes, all_genes))
    n_bg <- length(all_genes)

    # Fisher's exact test
    a <- n_overlap
    b <- n_deg - n_overlap
    c <- n_pw - n_overlap
    d <- n_bg - n_deg - c

    if (n_pw > 0 && n_overlap >= 1) {
      pval <- phyper(a - 1, n_pw, n_bg - n_pw, n_deg, lower.tail = FALSE)
      odds_ratio <- (a * d) / max((b * c), 1)
    } else {
      pval <- 1
      odds_ratio <- 0
    }

    data.frame(
      pathway = pw_name,
      n_overlap = n_overlap,
      n_deg = n_deg,
      n_pathway = n_pw,
      pval = pval,
      odds_ratio = odds_ratio,
      overlap_genes = paste(overlap, collapse = ",")
    )
  })

  result_df <- bind_rows(results)
  result_df$padj <- p.adjust(result_df$pval, method = "BH")
  result_df <- result_df %>% arrange(pval)
  return(result_df)
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

cat("=== Macrophage Pathway Enrichment ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load DEG results
# -----------------------------------------------------------------------------
cat("Step 1: Loading DEG results...\n")

deg_file <- file.path(output_dir, "macrophage_degs.csv")
if (!file.exists(deg_file)) {
  stop("Run 11_macrophage_deg_analysis.R first")
}

degs <- fread(deg_file)
cat("  Total DEGs:", nrow(degs), "\n")
cat("  Significant (adj.p < 0.05):", sum(degs$padj < 0.05, na.rm = TRUE), "\n")

# Get significant DEGs for enrichment
sig_up <- degs %>% filter(padj < 0.05, log2FoldChange > 0) %>% pull(gene)
sig_down <- degs %>% filter(padj < 0.05, log2FoldChange < 0) %>% pull(gene)
all_tested <- degs$gene

cat("  Upregulated genes:", length(sig_up), "\n")
cat("  Downregulated genes:", length(sig_down), "\n")
cat("  All tested genes:", length(all_tested), "\n")

# -----------------------------------------------------------------------------
# Step 2: Load gene sets (HALLMARK and BIOCARTA per manuscript)
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading gene sets...\n")

# HALLMARK pathways
hallmark_sets <- safe_msigdbr("Homo sapiens", "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
cat("  HALLMARK:", length(hallmark_list), "pathways\n")

# BIOCARTA pathways
biocarta_sets <- safe_msigdbr("Homo sapiens", "C2", "CP:BIOCARTA")
if (nrow(biocarta_sets) > 0) {
  biocarta_list <- split(biocarta_sets$gene_symbol, biocarta_sets$gs_name)
  cat("  BIOCARTA:", length(biocarta_list), "pathways\n")
} else {
  cat("  BIOCARTA: Not found, using HALLMARK only\n")
  biocarta_list <- list()
}

# Combine
all_pathways <- c(hallmark_list, biocarta_list)
cat("  Total pathways:", length(all_pathways), "\n")

# -----------------------------------------------------------------------------
# Step 3: Run hypergeometric enrichment
# -----------------------------------------------------------------------------
cat("\nStep 3: Running pathway enrichment...\n")

# Run enrichment on upregulated genes
cat("  Enrichment for upregulated genes...\n")
enrich_up <- run_simple_enrichment(sig_up, all_pathways, all_tested)
enrich_up$direction <- "Up"

# Run enrichment on downregulated genes
cat("  Enrichment for downregulated genes...\n")
enrich_down <- run_simple_enrichment(sig_down, all_pathways, all_tested)
enrich_down$direction <- "Down"

# Combine results
enrich_all <- bind_rows(enrich_up, enrich_down) %>%
  mutate(
    pathway_clean = gsub("HALLMARK_|BIOCARTA_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean)
  ) %>%
  arrange(pval)

cat("  Significant pathways (adj.p < 0.05):", sum(enrich_all$padj < 0.05), "\n")

# Check for TNFa and TGFb pathways (key for Figure 7C)
tnf_tgf <- enrich_all %>%
  filter(grepl("TNF|TGF|TNFA|TGFB", pathway, ignore.case = TRUE))
cat("\n  TNF/TGF pathways found:\n")
print(tnf_tgf %>% select(pathway, direction, n_overlap, pval, padj))

# Save full results
fwrite(enrich_all, file.path(output_dir, "macrophage_pathway_enrichment.csv"))
cat("\n  Saved macrophage_pathway_enrichment.csv\n")

# -----------------------------------------------------------------------------
# Step 4: Generate Figure 7C - Pathway heatmap
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating pathway heatmap...\n")

# Select top pathways by significance
top_pathways <- enrich_all %>%
  filter(padj < 0.1) %>%  # Slightly relaxed threshold to include more
  arrange(pval) %>%
  head(25)

if (nrow(top_pathways) < 5) {
  # If too few significant, take top 15 by p-value
  top_pathways <- enrich_all %>%
    arrange(pval) %>%
    head(15)
}

# Mark TNF/TGF pathways (red dots in manuscript)
top_pathways$is_tnf_tgf <- grepl("TNF|TGF", top_pathways$pathway, ignore.case = TRUE)

# Create signed score (positive for Up, negative for Down)
top_pathways$signed_score <- ifelse(top_pathways$direction == "Up",
                                     log2(top_pathways$odds_ratio + 1),
                                     -log2(top_pathways$odds_ratio + 1))

# Bar plot (similar to Panel C style)
p <- ggplot(top_pathways, aes(x = reorder(pathway_clean, signed_score), y = signed_score,
                               fill = ifelse(is_tnf_tgf, "TNF/TGF", direction))) +
  geom_bar(stat = "identity") +
  geom_point(data = top_pathways %>% filter(is_tnf_tgf),
             aes(x = reorder(pathway_clean, signed_score), y = signed_score),
             color = "red", size = 3, shape = 16) +
  coord_flip() +
  scale_fill_manual(values = c("TNF/TGF" = "#E41A1C", "Up" = "#377EB8", "Down" = "#4DAF4A"),
                    name = "Direction") +
  labs(
    title = "Pathway Enrichment: Macrophages Elderly vs Young",
    subtitle = "Red dots = TNF/TGF pathways",
    x = "",
    y = "Signed Log2(Odds Ratio)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "fig7c_pathway_barplot.png"),
       p, width = 10, height = 8, dpi = 300)
cat("  Saved fig7c_pathway_barplot.png\n")

# Also create a heatmap version showing enrichment direction
# Create a simple matrix for visualization
if (nrow(top_pathways) > 0 && any(!is.na(top_pathways$signed_score) & is.finite(top_pathways$signed_score))) {
  pathway_mat <- matrix(top_pathways$signed_score, ncol = 1,
                        dimnames = list(top_pathways$pathway_clean, "Elderly vs Young"))

  # Annotation for TNF/TGF
  row_ann <- data.frame(
    Category = ifelse(top_pathways$is_tnf_tgf, "TNF/TGF", "Other"),
    row.names = top_pathways$pathway_clean
  )

  ann_colors <- list(
    Category = c("TNF/TGF" = "#E41A1C", "Other" = "grey70")
  )

  max_val <- max(abs(pathway_mat[is.finite(pathway_mat)]), na.rm = TRUE)
  if (max_val > 0) {
    png(file.path(figures_dir, "fig7c_pathway_heatmap.png"),
        width = 8*300, height = 10*300, res = 300)
    pheatmap(
      pathway_mat,
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      annotation_row = row_ann,
      annotation_colors = ann_colors,
      main = "Pathway Enrichment: Macrophages Elderly vs Young",
      color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
      breaks = seq(-max_val, max_val, length.out = 101),
      fontsize_row = 11,
      fontsize_col = 12,
      fontsize = 12,
      show_colnames = TRUE
    )
    dev.off()
    cat("  Saved fig7c_pathway_heatmap.png\n")
  } else {
    cat("  Warning: No valid pathway scores for heatmap\n")
  }
} else {
  cat("  Warning: No pathways to visualize\n")
}

# -----------------------------------------------------------------------------
# Step 5: Show overlap genes for TNF/TGF pathways
# -----------------------------------------------------------------------------
cat("\nStep 5: Overlap genes for TNF/TGF pathways...\n")

for (i in seq_len(nrow(tnf_tgf))) {
  pw <- tnf_tgf$pathway[i]
  overlap_genes <- tnf_tgf$overlap_genes[i]

  if (!is.na(overlap_genes) && overlap_genes != "") {
    genes <- strsplit(overlap_genes, ",")[[1]]
    cat("\n  ", pw, " (", tnf_tgf$direction[i], "):\n", sep = "")
    cat("    Overlap genes:", paste(head(genes, 10), collapse = ", "), "\n")
  }
}

cat("\n=== Pathway enrichment complete ===\n")
cat("Next: Run 13_macrophage_cellphonedb.R\n")
