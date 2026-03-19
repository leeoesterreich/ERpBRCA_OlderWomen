#!/usr/bin/env Rscript
# 03_run_mica.R — Run MICA on GSVA pathway scores
# Adapted from Jian_MICA/scripts/02_mica_analysis.R
#
# Inputs:
#   - outputs/04_selected_path_aging_standard.RData  (from 01_prep_data.R)
#   - outputs/04_selected_path_aging_postM.RData     (from 01_prep_data.R)
#
# Outputs:
#   - outputs/04_selected_path_aging_result.RData
#   - outputs/mica_scheme1_results.csv
#   - outputs/mica_scheme2_results.csv
#   - figures/estrogene_boxplots_scheme2.pdf

set.seed(12345)

suppressPackageStartupMessages({
    library(parallel)
    library(MICA)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
})

# --- Paths ---
SCRIPT_DIR <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIG_DIR    <- file.path(SCRIPT_DIR, "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

cat("=== Stage 2: MICA Analysis ===\n")

# Load from outputs/ (Stage 1 output)
load(file.path(OUTPUT_DIR, "04_selected_path_aging_standard.RData"))
load(file.path(OUTPUT_DIR, "04_selected_path_aging_postM.RData"))

# Verify data structure
stopifnot(is.list(aging_selected_path1_data), length(aging_selected_path1_data) == 3)
stopifnot(is.list(aging_selected_path2_data), length(aging_selected_path2_data) == 3)
cat("Scheme 1 features:", nrow(aging_selected_path1_data[[1]]), "\n")
cat("Scheme 2 features:", nrow(aging_selected_path2_data[[1]]), "\n")

# --- Run MICA ---
# NOTE: p.threshold=1 means ALL features pass to post-hoc.
# This is as in the original script. Flagged in the audit.
cat("\nRunning MICA on Scheme 1 (Young/Middle-Aged/Elderly)...\n")
aging_selected_path1_result <- mica.full(
    study.data.matrix.list = aging_selected_path1_data,
    study.label.list = aging_selected_path1_meta,
    n.perm = 500,
    p.threshold = 1,
    n.parallel = 10
)

cat("Running MICA on Scheme 2 (Early/Middle/Elderly)...\n")
aging_selected_path2_result <- mica.full(
    study.data.matrix.list = aging_selected_path2_data,
    study.label.list = aging_selected_path2_meta,
    n.perm = 500,
    p.threshold = 1,
    n.parallel = 10
)

save(aging_selected_path1_result, aging_selected_path2_result,
     file = file.path(OUTPUT_DIR, "04_selected_path_aging_result.RData"))

# --- Print results ---
cat("\n=== Scheme 1 Results ===\n")
feature_names1 <- rownames(aging_selected_path1_data[[1]])
results1 <- data.frame(
    feature = feature_names1,
    gmi_stat = round(aging_selected_path1_result$gmi.stat, 5),
    gmi_pval = round(aging_selected_path1_result$gmi.pval, 5),
    minmcc_stat = round(aging_selected_path1_result$minmcc.stat, 5),
    minmcc_pval = round(aging_selected_path1_result$minmcc.pval, 5)
)
print(results1)

cat("\n=== Scheme 2 Results ===\n")
feature_names2 <- rownames(aging_selected_path2_data[[1]])
results2 <- data.frame(
    feature = feature_names2,
    gmi_stat = round(aging_selected_path2_result$gmi.stat, 5),
    gmi_pval = round(aging_selected_path2_result$gmi.pval, 5),
    minmcc_stat = round(aging_selected_path2_result$minmcc.stat, 5),
    minmcc_pval = round(aging_selected_path2_result$minmcc.pval, 5)
)
print(results2)

# Save tables as CSV
write.csv(results1, file.path(OUTPUT_DIR, "mica_scheme1_results.csv"), row.names = FALSE)
write.csv(results2, file.path(OUTPUT_DIR, "mica_scheme2_results.csv"), row.names = FALSE)

# --- Generate EstroGene boxplots (Scheme 2) ---
cat("\nGenerating EstroGene boxplots...\n")
gene_boxplot <- function(gene, data_list, meta_list, levels_order) {
    study.tbl <- data.frame(freq = sapply(meta_list, function(x) length(x)))
    expr <- data.frame(
        expr = unlist(sapply(data_list, function(x) scale(t(x)[, gene]))),
        age = unlist(sapply(meta_list, function(x) factor(x, levels = levels_order))),
        study = rep(rownames(study.tbl), study.tbl$freq)
    )
    expr$study <- dplyr::recode(expr$study, "metabric" = "METABRIC", "scanb" = "SCAN-B", "tcga" = "TCGA")
    fig <- ggplot(expr, aes(x = age, y = expr)) +
        geom_boxplot(fill = "gray80") +
        facet_wrap(~study) +
        theme_bw() +
        theme(legend.position = "none",
              strip.text.x = element_text(face = "bold", size = 20),
              axis.title = element_text(size = 20, face = "bold"),
              axis.text = element_text(size = 10)) +
        xlab("") + ylab("") +
        coord_cartesian(ylim = c(-2, 2)) +
        ggtitle(gene)
    return(fig)
}

fig_estrogene <- gene_boxplot("EstroGene_Early", aging_selected_path2_data, aging_selected_path2_meta, c("Early", "Middle", "Elderly")) +
    gene_boxplot("EstroGene_Mid", aging_selected_path2_data, aging_selected_path2_meta, c("Early", "Middle", "Elderly")) +
    gene_boxplot("EstroGene_Late", aging_selected_path2_data, aging_selected_path2_meta, c("Early", "Middle", "Elderly")) +
    plot_layout(nrow = 3)
ggsave(file.path(FIG_DIR, "estrogene_boxplots_scheme2.pdf"), fig_estrogene,
       width = 20, height = 30, units = "cm")

cat("\nStage 2 complete.\n")
