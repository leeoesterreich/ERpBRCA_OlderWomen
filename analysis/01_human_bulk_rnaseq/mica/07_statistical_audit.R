#!/usr/bin/env Rscript
# 07_statistical_audit.R — Independent statistical audit of MICA results
# Adapted from Jian_MICA/scripts/05_statistical_audit.R
#
# Checks: permutation null, p-value calibration, effect sizes, plotting issues,
# cross-cohort consistency, and MICA method validity.
#
# Inputs:
#   - outputs/04_selected_path_aging_standard.RData
#   - outputs/04_selected_path_aging_postM.RData
#   - outputs/04_selected_path_aging_result.RData
#
# Outputs:
#   - outputs/statistical_audit.json
#   - outputs/audit_summary.md
#   - figures/audit_null_distributions.pdf

set.seed(12345)

suppressPackageStartupMessages({
    library(MICA)
    library(parallel)
    library(jsonlite)
})

# --- Paths ---
SCRIPT_DIR <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIG_DIR    <- file.path(SCRIPT_DIR, "figures")
MICA_INTERNAL <- Sys.getenv("MICA_INTERNAL", unset = file.path(SCRIPT_DIR, "..", "..", "..", "MICA", "R", "00_internal_functions.R"))

cat("=== Stage 5: Statistical Audit ===\n")

findings <- list()
add_finding <- function(severity, category, title, detail) {
    findings[[length(findings) + 1]] <<- list(
        severity = severity,
        category = category,
        title = title,
        detail = detail
    )
    cat(sprintf("[%s] %s: %s\n", toupper(severity), category, title))
}

# ============================================================
# 1. Load data and results
# ============================================================
load(file.path(OUTPUT_DIR, "04_selected_path_aging_standard.RData"))
load(file.path(OUTPUT_DIR, "04_selected_path_aging_postM.RData"))
load(file.path(OUTPUT_DIR, "04_selected_path_aging_result.RData"))

# ============================================================
# 2. MICA Method Audit: Pooled null p-value investigation
# ============================================================
cat("\n--- Audit: Pooled null p-values ---\n")
cat("Re-running permutation to examine null distributions...\n")
n_perm_check <- 200
perm_label_check <- list()
for (p in 1:n_perm_check) {
    perm_label_check[[p]] <- lapply(aging_selected_path2_meta, function(x) sample(x))
}

# Source the internal functions for inspection
source(MICA_INTERNAL)

perm_gmi_matrix <- matrix(nrow = n_perm_check,
                           ncol = nrow(aging_selected_path2_data[[1]]))
for (p in 1:n_perm_check) {
    perm_gmi_matrix[p, ] <- gmi.matr(aging_selected_path2_data, perm_label_check[[p]])
}

null_means <- colMeans(perm_gmi_matrix)
null_sds <- apply(perm_gmi_matrix, 2, sd)
feature_names <- rownames(aging_selected_path2_data[[1]])
cat("Null distribution summary:\n")
for (i in seq_along(feature_names)) {
    cat(sprintf("  %s: mean=%.4f, sd=%.4f\n", feature_names[i], null_means[i], null_sds[i]))
}

sd_ratio <- max(null_sds) / min(null_sds)
if (sd_ratio > 2) {
    add_finding("critical", "MICA_method",
                "Pooled null p-values: feature null distributions differ substantially",
                sprintf("Max/min SD ratio = %.2f. Pooled null is biased.", sd_ratio))
} else {
    add_finding("suggestion", "MICA_method",
                "Pooled null p-values: feature null distributions are reasonably similar",
                sprintf("Max/min SD ratio = %.2f. Pooled null is acceptable but unconventional.", sd_ratio))
}

# ============================================================
# 3. +1 correction check
# ============================================================
cat("\n--- Audit: +1 correction ---\n")
for (scheme in c("path1", "path2")) {
    result <- if (scheme == "path1") aging_selected_path1_result else aging_selected_path2_result
    zero_gmi <- sum(result$gmi.pval == 0)
    zero_minmcc <- sum(result$minmcc.pval == 0)
    if (zero_gmi > 0 || zero_minmcc > 0) {
        add_finding("warning", "MICA_method",
                    sprintf("Exact zero p-values in %s", scheme),
                    sprintf("gMI+ zeros: %d, minMCC zeros: %d. Missing +1 correction.", zero_gmi, zero_minmcc))
    }
}

# ============================================================
# 4. Effect sizes (Cohen's d) for key pathways
# ============================================================
cat("\n--- Audit: Effect sizes ---\n")
cohens_d <- function(x, y) {
    nx <- length(x); ny <- length(y)
    pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
    if (pooled_sd == 0) return(NA)
    (mean(x) - mean(y)) / pooled_sd
}

for (scheme_name in c("Scheme1", "Scheme2")) {
    data_list <- if (scheme_name == "Scheme1") aging_selected_path1_data else aging_selected_path2_data
    meta_list <- if (scheme_name == "Scheme1") aging_selected_path1_meta else aging_selected_path2_meta
    levels_order <- if (scheme_name == "Scheme1") c("Young", "Middle-Aged", "Elderly") else c("Early", "Middle", "Elderly")

    cat(sprintf("\n  %s effect sizes (Elderly vs youngest group):\n", scheme_name))
    for (cohort in names(data_list)) {
        features <- rownames(data_list[[cohort]])
        labels <- meta_list[[cohort]]
        youngest <- levels_order[1]
        oldest <- levels_order[length(levels_order)]

        for (feat in features) {
            vals <- data_list[[cohort]][feat, ]
            x <- vals[labels == youngest]
            y <- vals[labels == oldest]
            d <- cohens_d(as.numeric(x), as.numeric(y))
            cat(sprintf("    %s / %s: d = %.3f (n_young=%d, n_old=%d)\n",
                        cohort, feat, d, length(x), length(y)))
        }
    }
}

# ============================================================
# 5. Class label alignment verification
# ============================================================
cat("\n--- Audit: Class label alignment ---\n")
for (scheme_name in c("Scheme1", "Scheme2")) {
    meta_list <- if (scheme_name == "Scheme1") aging_selected_path1_meta else aging_selected_path2_meta
    labels_per_study <- lapply(meta_list, function(x) sort(unique(x)))
    all_same <- all(sapply(labels_per_study[-1], function(x) identical(x, labels_per_study[[1]])))
    if (all_same) {
        cat(sprintf("  %s: Labels consistent across studies: %s\n",
                    scheme_name, paste(labels_per_study[[1]], collapse = ", ")))
    } else {
        add_finding("critical", "data_integrity",
                    sprintf("Class labels misaligned in %s", scheme_name),
                    paste(capture.output(print(labels_per_study)), collapse = "\n"))
    }
}

# ============================================================
# 6. Sample size adequacy
# ============================================================
cat("\n--- Audit: Sample sizes ---\n")
for (scheme_name in c("Scheme1", "Scheme2")) {
    meta_list <- if (scheme_name == "Scheme1") aging_selected_path1_meta else aging_selected_path2_meta
    for (cohort in names(meta_list)) {
        tbl <- table(meta_list[[cohort]])
        small_groups <- tbl[tbl < 10]
        if (length(small_groups) > 0) {
            add_finding("warning", "sample_size",
                        sprintf("%s/%s has small groups", scheme_name, cohort),
                        sprintf("Groups with n < 10: %s",
                                paste(names(small_groups), small_groups, sep = "=", collapse = ", ")))
        }
    }
}

# ============================================================
# 7. gMI+ numerical stability check
# ============================================================
cat("\n--- Audit: gMI+ numerical stability ---\n")
for (scheme_name in c("Scheme1", "Scheme2")) {
    result <- if (scheme_name == "Scheme1") aging_selected_path1_result else aging_selected_path2_result
    if (any(is.na(result$gmi.stat) | is.infinite(result$gmi.stat))) {
        add_finding("critical", "MICA_method",
                    sprintf("gMI+ contains NA/Inf in %s", scheme_name),
                    paste(result$gmi.stat, collapse = ", "))
    }
    if (any(is.na(result$minmcc.stat) | is.infinite(result$minmcc.stat))) {
        add_finding("critical", "MICA_method",
                    sprintf("minMCC contains NA/Inf in %s", scheme_name),
                    paste(result$minmcc.stat, collapse = ", "))
    }
}

# ============================================================
# 8. Post-hoc pairwise multiple testing audit
# ============================================================
cat("\n--- Audit: Post-hoc pairwise multiple testing ---\n")
for (scheme_name in c("Scheme1", "Scheme2")) {
    result <- if (scheme_name == "Scheme1") aging_selected_path1_result else aging_selected_path2_result
    if (!is.null(result$pairwise.tbl) && length(result$pairwise.tbl) > 0) {
        cat(sprintf("  %s: %d pairwise tables found\n", scheme_name, length(result$pairwise.tbl)))
        for (i in seq_along(result$pairwise.tbl)) {
            tbl <- result$pairwise.tbl[[i]]
            if (!"padj" %in% colnames(tbl) && !"qval" %in% colnames(tbl)) {
                n_pairs <- nrow(tbl)
                n_sig_raw <- sum(tbl$pval < 0.05, na.rm = TRUE)
                add_finding("critical", "MICA_method",
                            sprintf("Post-hoc pairwise p-values uncorrected (feature %d, %s)", i, scheme_name),
                            sprintf("Raw p-values for %d pairs without BH/Bonferroni. %d/%d significant at raw p<0.05.",
                                    n_pairs, n_sig_raw, n_pairs))
            }
        }
    } else {
        cat(sprintf("  %s: No pairwise tables\n", scheme_name))
    }
}

# ============================================================
# 9. Confounder sensitivity flags
# ============================================================
cat("\n--- Audit: Confounder sensitivity ---\n")
add_finding("warning", "confounders",
            "No tumor purity adjustment",
            "GSVA scores computed on bulk expression without tumor purity correction.")
add_finding("warning", "confounders",
            "No platform/batch adjustment across cohorts",
            "TCGA (RNA-seq TPM), METABRIC (microarray, TMM-normalized CPM), SCAN-B (RNA-seq TPM) use different platforms.")
add_finding("warning", "confounders",
            "No menopausal status adjustment",
            "Age groups span pre- and post-menopausal women. Estrogen changes at menopause are a major confounder.")
add_finding("suggestion", "confounders",
            "Subtype composition not verified across age groups",
            "If Luminal A vs B proportions shift with age, subtype-associated differences could mimic age effects.")

# ============================================================
# 10. Permutation null calibration plots
# ============================================================
cat("\n--- Audit: Null calibration plots ---\n")
pdf(file.path(FIG_DIR, "audit_null_distributions.pdf"), width = 12, height = 8)
par(mfrow = c(2, ceiling(ncol(perm_gmi_matrix) / 2)))
for (i in 1:ncol(perm_gmi_matrix)) {
    hist(perm_gmi_matrix[, i], breaks = 30, main = feature_names[i],
         xlab = "gMI+ (null)", col = "lightblue", border = "white")
    abline(v = aging_selected_path2_result$gmi.stat[i], col = "red", lwd = 2)
}
dev.off()
cat("Saved figures/audit_null_distributions.pdf\n")

# ============================================================
# 11. Boxplot issue documentation
# ============================================================
cat("\n--- Audit: Known plotting issues ---\n")
add_finding("warning", "visualization",
            "Boxplot midline shows mean, not median",
            "Original boxplots used aes(middle = mean(expr)). Monorepo versions use standard median.")
add_finding("warning", "visualization",
            "ylim clips outliers silently",
            "Boxplots use coord_cartesian(ylim=...) which zooms without removing data.")
add_finding("warning", "visualization",
            "No confidence intervals or error bars",
            "All boxplots show raw distributions without CIs or significance brackets.")

# ============================================================
# Save findings
# ============================================================
cat("\n\n=== AUDIT SUMMARY ===\n")
cat(sprintf("Total findings: %d\n", length(findings)))
for (sev in c("critical", "warning", "suggestion")) {
    n <- sum(sapply(findings, function(f) f$severity == sev))
    cat(sprintf("  %s: %d\n", toupper(sev), n))
}

write(toJSON(findings, auto_unbox = TRUE, pretty = TRUE),
      file.path(OUTPUT_DIR, "statistical_audit.json"))

sink(file.path(OUTPUT_DIR, "audit_summary.md"))
cat("# MICA Statistical Audit Summary\n\n")
cat(sprintf("Date: %s\n\n", Sys.Date()))
for (sev in c("critical", "warning", "suggestion")) {
    items <- Filter(function(f) f$severity == sev, findings)
    if (length(items) > 0) {
        cat(sprintf("## %s (%d)\n\n", toupper(sev), length(items)))
        for (item in items) {
            cat(sprintf("### %s: %s\n\n", item$category, item$title))
            cat(sprintf("%s\n\n", item$detail))
        }
    }
}
sink()

cat("\nStage 5 complete. Results saved to outputs/\n")
