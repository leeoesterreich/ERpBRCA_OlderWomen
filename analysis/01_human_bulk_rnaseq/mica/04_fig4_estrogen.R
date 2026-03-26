#!/usr/bin/env Rscript
# 04_fig4_estrogen.R — Generate manuscript figures 4C, 4H
# Adapted from Jian_MICA/scripts/03_figure_generation.R (estrogen sections)
#
# Figure 4C: HSD17B gene heatmap across age groups
# Figure 4H: Estrogen pathway boxplots (GO BP, WikiPathways)
#
# Inputs:
#   - data/input/raw_data/bulk_selected_gene_aging_v2.RData
#   - data/input/raw_data/02_selected_path_aging_standard.RData
#
# Outputs:
#   - figures/fig4c.pdf
#   - figures/fig4h.pdf

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

cat("=== Figure 4: Estrogen / HSD17B ===\n")

# ============================================================
# Figure 4C: HSD17B Gene Heatmap
# ============================================================
cat("\n--- Figure 4C ---\n")
load(file.path(DATA_DIR, "raw_data", "bulk_selected_gene_aging_v2.RData"))
stopifnot("aging_selected_gene_data" %in% ls(), "aging_selected_gene_meta" %in% ls())

heatmap_median <- function(study.data.list, study.label.list) {
    median_df <- list()
    for (i in 1:length(study.data.list)) {
        median_df[[i]] <- aggregate(x = study.data.list[[i]],
                                    by = list(study.label.list[[i]]),
                                    FUN = median) %>%
            column_to_rownames("Group.1")
    }
    median_df <- lapply(median_df, t)

    # Order genes by max fold change across datasets
    fold_changes <- lapply(median_df, function(x) {
        apply(x, 1, function(row) max(row) - min(row))
    })
    fold_changes_df <- do.call(cbind, fold_changes)
    rownames(fold_changes_df) <- rownames(median_df[[1]])
    max_fold_change <- apply(fold_changes_df, 1, max)
    gene_order <- names(sort(max_fold_change, decreasing = TRUE))
    median_df <- lapply(median_df, function(x) x[gene_order, ])

    median_df_merge <- do.call("cbind", median_df)
    data_source <- rep(names(study.data.list), each = ncol(median_df[[1]]))
    data_source <- dplyr::recode(data_source, "metabric" = "METABRIC", "scanb" = "SCAN-B", "tcga" = "TCGA")

    h <- Heatmap(t(median_df_merge), name = "Expression",
                 row_split = data_source,
                 cluster_columns = FALSE,
                 cluster_rows = FALSE,
                 show_column_names = TRUE,
                 column_names_side = "bottom",
                 row_names_max_width = unit(9, "cm"),
                 row_names_gp = gpar(fontsize = 8))
    return(h)
}

gene_order_initial <- c("ESR1", "GREB1", "PGR", "SAA1", "RAB19", "KRT37", "TRPM8",
                         "CYP19A1", "HSD17B1", "HSD17B7", "HSD17B12", "HSD17B2", "HSD17B10", "HSD17B14")
aging_selected_gene_data2 <- lapply(aging_selected_gene_data, function(x) x[, gene_order_initial])
fig4c <- heatmap_median(
    lapply(aging_selected_gene_data2, function(x) scale(x)),
    lapply(aging_selected_gene_meta, function(x) factor(x, levels = c("Elderly", "Middle-Aged", "Young")))
)
pdf(file.path(FIG_DIR, "fig4c.pdf"), width = 8, height = 6)
draw(fig4c)
dev.off()
cat("Saved figures/fig4c.pdf\n")

# ============================================================
# Figure 4H: Estrogen Pathway Boxplots
# ============================================================
cat("\n--- Figure 4H ---\n")
rm(list = setdiff(ls(), c("SCRIPT_DIR", "DATA_DIR", "FIG_DIR")))
load(file.path(DATA_DIR, "raw_data", "02_selected_path_aging_standard.RData"))
stopifnot("aging_selected_path_data" %in% ls(), "aging_selected_path_meta" %in% ls())

path_boxplot <- function(path, display_name = NULL) {
    study.tbl <- data.frame(freq = sapply(aging_selected_path_meta, function(x) length(x)))
    expr <- data.frame(
        expr = unlist(sapply(aging_selected_path_data, function(x) scale(x[, path]))),
        age = unlist(sapply(aging_selected_path_meta, function(x) factor(x, levels = c("Young", "Middle-Aged", "Elderly")))),
        study = rep(rownames(study.tbl), study.tbl$freq)
    )
    expr$study <- dplyr::recode(expr$study, "metabric" = "METABRIC", "scanb" = "SCAN-B", "tcga" = "TCGA")
    fig <- ggplot(expr, aes(x = age, y = expr, fill = age)) +
        geom_boxplot() +
        facet_wrap(~study) +
        theme_bw() +
        theme(legend.position = "none",
              strip.text.x = element_text(face = "bold", size = 20),
              axis.title = element_text(size = 20, face = "bold"),
              axis.text = element_text(size = 10)) +
        scale_fill_manual(values = c("#fee0d2", "#fc9272", "#de2d26")) +
        xlab("") + ylab("") +
        coord_cartesian(ylim = c(-2, 2)) +
        ggtitle(ifelse(is.null(display_name), path, display_name))
    return(fig)
}

fig4h <- path_boxplot("GOBP_RESPONSE_TO_ESTROGEN", "Response to Estrogen (GO BP)") +
    path_boxplot("WP_ESTROGEN_RECEPTOR_PATHWAY", "Estrogen Receptor Pathway (WikiPathways)") +
    plot_layout(nrow = 2)
ggsave(file.path(FIG_DIR, "fig4h.pdf"), fig4h, width = 8, height = 6)
cat("Saved figures/fig4h.pdf\n")

cat("\nFigure 4 generation complete.\n")
