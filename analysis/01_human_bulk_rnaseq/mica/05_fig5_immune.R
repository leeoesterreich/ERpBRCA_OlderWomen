#!/usr/bin/env Rscript
# 05_fig5_immune.R — Generate manuscript figures 5C, 5E
# Adapted from Jian_MICA/scripts/03_figure_generation.R (immune sections)
#
# Figure 5C: Inflammatory pathway heatmap
# Figure 5E: MCP immune deconvolution boxplots
#
# Inputs:
#   - data/input/raw_data/aging_bulk_rnaseq_result.RData
#   - data/input/raw_data/bulk_pathway_aging.RData
#   - data/input/raw_data/aging_bulk_result_neil_edited_20230205.xlsx
#   - data/input/raw_data/04_immune_mcp_standard_results.RData
#   - data/input/raw_data/04_immune_mcp_aging_standard.RData
#
# Outputs:
#   - figures/fig5c.pdf
#   - figures/fig5e.pdf

set.seed(12345)

suppressPackageStartupMessages({
    library(edgeR)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(tibble)
    library(stringr)
    library(ComplexHeatmap)
    library(circlize)
    library(patchwork)
})

# --- Paths ---
SCRIPT_DIR <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
)
DATA_DIR   <- file.path(SCRIPT_DIR, "data", "input")
FIG_DIR    <- file.path(SCRIPT_DIR, "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

cat("=== Figure 5: Immune / Inflammation ===\n")

# ============================================================
# Figure 5C: Inflammatory Pathway Heatmap
# ============================================================
cat("\n--- Figure 5C ---\n")
load(file.path(DATA_DIR, "raw_data", "aging_bulk_rnaseq_result.RData"))
load(file.path(DATA_DIR, "raw_data", "bulk_pathway_aging.RData"))
stopifnot("aging_data" %in% ls(), "aging_label" %in% ls())

heatmap_median_5c <- function(study.data.list, study.label.list) {
    median_df <- list()
    for (i in 1:length(study.data.list)) {
        median_df[[i]] <- aggregate(x = study.data.list[[i]],
                                    by = list(study.label.list[[i]]),
                                    FUN = median) %>%
            column_to_rownames("Group.1")
    }
    median_df <- lapply(median_df, t)
    median_df_merge <- do.call("cbind", median_df)
    data_source <- rep(names(study.data.list), each = 3)
    data_source <- dplyr::recode(data_source, "metabric" = "METABRIC", "scanb" = "SCAN-B", "tcga" = "TCGA")

    f1 <- colorRamp2(seq(-0.4, 0.4, length = 3), c("blue", "white", "red"))
    h <- Heatmap(t(median_df_merge), name = "Expression", col = f1,
                 row_split = data_source,
                 cluster_columns = TRUE,
                 cluster_rows = FALSE,
                 column_names_max_height = unit(7, "cm"),
                 column_names_gp = gpar(fontsize = 8))
    return(h)
}

inflammatory <- openxlsx::read.xlsx(
    file.path(DATA_DIR, "raw_data", "aging_bulk_result_neil_edited_20230205.xlsx"), sheet = 1)
pathway_col <- colnames(inflammatory)[1]
infl_col <- colnames(inflammatory)[grep("infl", colnames(inflammatory), ignore.case = TRUE)]
inflammatory_path <- inflammatory[[pathway_col]][which(inflammatory[[infl_col]] == 1)]
cat("Inflammatory pathways selected:", length(inflammatory_path), "\n")

fig5c <- heatmap_median_5c(lapply(aging_data, function(x) x[, inflammatory_path]), aging_label)
pdf(file.path(FIG_DIR, "fig5c.pdf"), width = 10, height = 7)
draw(fig5c)
dev.off()
cat("Saved figures/fig5c.pdf\n")

# ============================================================
# Figure 5E: MCP Immune Deconvolution Boxplots
# ============================================================
cat("\n--- Figure 5E ---\n")
rm(list = setdiff(ls(), c("SCRIPT_DIR", "DATA_DIR", "FIG_DIR")))
load(file.path(DATA_DIR, "raw_data", "04_immune_mcp_standard_results.RData"))
load(file.path(DATA_DIR, "raw_data", "04_immune_mcp_aging_standard.RData"))
stopifnot("aging_selected_gene_meta" %in% ls(), "aging_selected_gene_mcp" %in% ls())

cell_boxplot <- function(cell, display_name = NULL) {
    study.tbl <- data.frame(freq = sapply(aging_selected_gene_meta, function(x) length(x)))
    expr <- data.frame(
        expr = unlist(sapply(aging_selected_gene_mcp, function(x) scale(x[, cell]))),
        age = unlist(sapply(aging_selected_gene_meta, function(x) factor(x, levels = c("Young", "Middle-Aged", "Elderly")))),
        study = rep(rownames(study.tbl), study.tbl$freq)
    )
    expr$study <- dplyr::recode(expr$study, "metabric" = "METABRIC", "scanb" = "SCAN-B", "tcga" = "TCGA")
    fig <- ggplot(expr, aes(x = age, y = expr, fill = age)) +
        geom_boxplot() +
        facet_wrap(~study) +
        theme_bw() +
        theme(legend.position = "none",
              strip.text.x = element_text(face = "bold", size = 15),
              axis.title = element_text(size = 15, face = "bold"),
              axis.text = element_text(size = 8)) +
        scale_fill_manual(values = c("#fee0d2", "#fc9272", "#de2d26")) +
        xlab("") + ylab("MCP-counter Score") +
        coord_cartesian(ylim = c(-1, 3)) +
        ggtitle(ifelse(is.null(display_name), cell, display_name))
    return(fig)
}

fig5e <- cell_boxplot("CD8 T cells", "CD8 T Cells") +
    cell_boxplot("Cytotoxic lymphocytes", "Cytotoxic Lymphocytes") +
    cell_boxplot("Monocytic lineage", "Monocytic Lineage") +
    cell_boxplot("Myeloid dendritic cells", "Myeloid Dendritic Cells") +
    plot_layout(nrow = 2)
ggsave(file.path(FIG_DIR, "fig5e.pdf"), fig5e, width = 12, height = 10)
cat("Saved figures/fig5e.pdf\n")

cat("\nFigure 5 generation complete.\n")
