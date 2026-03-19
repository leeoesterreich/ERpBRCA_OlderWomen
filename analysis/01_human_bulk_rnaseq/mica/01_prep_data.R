#!/usr/bin/env Rscript
# 01_prep_data.R — Load cohort data, create age groups, compute GSVA scores
# Merged from: stub 01_prep_data.R (age group functions) +
#              Jian_MICA/scripts/01_estrogene_gsva.R (GSVA implementation)
#
# Inputs:
#   - data/input/ilc_clean_data.RData
#   - data/input/clinical/METABRIC_Clinical_Info.csv
#   - data/input/clinical/TCGA_Key.csv
#   - data/input/database/primary/EstroGene_Signatures.xlsx
#   - data/input/gmt/LI_ESTROGENE_EARLY_E2_RESPONSE_UP.v2025.1.Hs.gmt
#   - data/input/gmt/LI_ESTROGENE_LATE_E2_RESPONSE_UP.v2025.1.Hs.gmt
#
# Outputs:
#   - outputs/04_selected_path_aging_standard.RData  (Scheme 1 GSVA)
#   - outputs/04_selected_path_aging_postM.RData     (Scheme 2 GSVA)

set.seed(12345)

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(tibble)
    library(stringr)
    library(purrr)
    library(GSVA)
    library(edgeR)
})

# --- Paths ---
SCRIPT_DIR <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
)
DATA_DIR   <- file.path(SCRIPT_DIR, "data", "input")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("=== Stage 1: Data Preparation & GSVA Scoring ===\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")

# --- Age group functions ---
create_age_groups <- function(age) {
    dplyr::case_when(
        age >= 35 & age <= 45 ~ "Young",
        age >= 55 & age <= 69 ~ "Middle-Aged",
        age >= 70             ~ "Elderly",
        TRUE ~ "ND"
    )
}

create_age_groups_postm <- function(age) {
    dplyr::case_when(
        age >= 55 & age < 60 ~ "Early",
        age >= 60 & age < 70 ~ "Middle",
        age >= 70            ~ "Elderly",
        TRUE ~ "ND"
    )
}

# --- Load EstroGene signatures ---
estrogene_geneset <- openxlsx::read.xlsx(
    file.path(DATA_DIR, "database", "primary", "EstroGene_Signatures.xlsx"), sheet = 1)
estrogene_list <- list(
    EstroGene_Early = na.omit(estrogene_geneset$EstroGene_Early),
    EstroGene_Mid   = na.omit(estrogene_geneset$EstroGene_Mid),
    EstroGene_Late  = na.omit(estrogene_geneset$EstroGene_Late)
)
cat("EstroGene signatures loaded:",
    paste(names(estrogene_list), sapply(estrogene_list, length), sep = "=", collapse = ", "), "\n")

# --- Load LI_ESTROGENE GMT files ---
parse_gmt <- function(path) {
    line <- readLines(path, n = 1)
    fields <- strsplit(line, "\t")[[1]]
    fields[-(1:2)]  # skip name and URL
}
estrogene_list[["LI_ESTROGENE_Early_E2_Response_Up"]] <- parse_gmt(
    file.path(DATA_DIR, "gmt", "LI_ESTROGENE_EARLY_E2_RESPONSE_UP.v2025.1.Hs.gmt"))
estrogene_list[["LI_ESTROGENE_Late_E2_Response_Up"]]  <- parse_gmt(
    file.path(DATA_DIR, "gmt", "LI_ESTROGENE_LATE_E2_RESPONSE_UP.v2025.1.Hs.gmt"))
cat("Added LI_ESTROGENE signatures:",
    length(estrogene_list$LI_ESTROGENE_Early_E2_Response_Up), "early,",
    length(estrogene_list$LI_ESTROGENE_Late_E2_Response_Up), "late genes\n")

# --- Load ILC clean data ---
methods::setClass(
    "BCdata",
    slots = c(
        model = "character",
        count = "data.frame",
        tpm = "data.frame",
        microarray = "data.frame",
        clinic = "data.frame",
        annotation = "character"
    )
)
load(file.path(DATA_DIR, "ilc_clean_data.RData"))
stopifnot("scanb" %in% ls(), "metabric" %in% ls(), "tcga" %in% ls())
cat("Loaded ilc_clean_data.RData: scanb, metabric, tcga objects present.\n")

# ============================================================
# SCAN-B
# ============================================================
cat("\n--- Processing SCAN-B ---\n")
scanb_meta <- scanb@clinic %>%
    mutate(Age_Class  = create_age_groups(age),
           Age_Class2 = create_age_groups_postm(age)) %>%
    filter(clinic_grp == "ERpHER2n")

cat("SCAN-B ER+/HER2- patients:", nrow(scanb_meta), "\n")
cat("  Scheme 1 (Young/Middle-Aged/Elderly):", table(scanb_meta$Age_Class), "\n")
cat("  Scheme 2 (Early/Middle/Elderly):", table(scanb_meta$Age_Class2), "\n")

scanb_gsva <- gsva(gsvaParam(as.matrix(scanb@tpm)[, rownames(scanb_meta)], estrogene_list))
scanb_aging <- list(pathway = t(scanb_gsva), meta = scanb_meta)

scanb_selection1 <- scanb_aging$meta %>% filter(Age_Class != "ND") %>% rownames()
scanb_path1 <- list(pathway = scanb_aging$pathway[scanb_selection1, ],
                    meta = scanb_aging$meta[scanb_selection1, "Age_Class"])
scanb_selection2 <- scanb_aging$meta %>% filter(Age_Class2 != "ND") %>% rownames()
scanb_path2 <- list(pathway = scanb_aging$pathway[scanb_selection2, ],
                    meta = scanb_aging$meta[scanb_selection2, "Age_Class2"])
cat("  Scheme 1 n:", length(scanb_selection1), "  Scheme 2 n:", length(scanb_selection2), "\n")

# ============================================================
# METABRIC
# ============================================================
cat("\n--- Processing METABRIC ---\n")
metabric_list <- list()
metabric_list$meta <- data.table::fread(
    file.path(DATA_DIR, "clinical", "METABRIC_Clinical_Info.csv")) %>%
    as.data.frame() %>%
    column_to_rownames(var = "V1") %>%
    mutate(Age_Class  = create_age_groups(age_at_diagnosis),
           Age_Class2 = create_age_groups_postm(age_at_diagnosis))

stopifnot("age_at_diagnosis" %in% colnames(metabric_list$meta))

metabric_aging_selection <- metabric_list$meta %>%
    filter(ER == "Positive" & HER2_IHC_status %in% c(0, 1)) %>%
    rownames()
metabric_aging_selection <- intersect(metabric_aging_selection, colnames(metabric@microarray))
cat("METABRIC ER+/HER2- patients (with expression):", length(metabric_aging_selection), "\n")

metabric_norm <- DGEList(metabric@microarray[-which(apply(metabric@microarray, 1, function(x) sum(is.na(x))) > 0), ])
metabric_norm <- calcNormFactors(metabric_norm, method = "TMM")
metabric_norm <- cpm(metabric_norm, log = TRUE)
metabric_gsva <- gsva(gsvaParam(as.matrix(metabric_norm)[, metabric_aging_selection], estrogene_list))
metabric_aging <- list(pathway = t(metabric_gsva),
                       meta = metabric_list$meta[metabric_aging_selection, ])

metabric_selection1 <- metabric_aging$meta %>% filter(Age_Class != "ND") %>% rownames()
metabric_path1 <- list(pathway = metabric_aging$pathway[metabric_selection1, ],
                       meta = metabric_aging$meta[metabric_selection1, "Age_Class"])
metabric_selection2 <- metabric_aging$meta %>% filter(Age_Class2 != "ND") %>% rownames()
metabric_path2 <- list(pathway = metabric_aging$pathway[metabric_selection2, ],
                       meta = metabric_aging$meta[metabric_selection2, "Age_Class2"])
cat("  Scheme 1 n:", length(metabric_selection1), "  Scheme 2 n:", length(metabric_selection2), "\n")

# ============================================================
# TCGA
# ============================================================
cat("\n--- Processing TCGA ---\n")
tcga_list <- list()
tcga_list$meta <- data.table::fread(
    file.path(DATA_DIR, "clinical", "TCGA_Key.csv")) %>%
    as.data.frame() %>%
    column_to_rownames(var = "V1") %>%
    mutate(Age_Class  = create_age_groups(Age),
           Age_Class2 = create_age_groups_postm(Age))
tcga_aging_selection <- tcga_list$meta %>%
    filter(ER == "Positive" & HER2 == "Negative") %>%
    rownames()
cat("TCGA ER+/HER2- patients:", length(tcga_aging_selection), "\n")

tcga_aging <- list()
tcga_aging$tpm <- tcga@tpm
colnames(tcga_aging$tpm) <- substr(colnames(tcga_aging$tpm), 1, 12)
tcga_aging$pathway <- t(GSVA::gsva(gsvaParam(
    as.matrix(tcga_aging$tpm)[, tcga_aging_selection], estrogene_list)))
tcga_aging$meta <- tcga_list$meta[tcga_aging_selection, ]

tcga_selection1 <- tcga_aging$meta %>% filter(Age_Class != "ND") %>% rownames()
tcga_path1 <- list(pathway = tcga_aging$pathway[tcga_selection1, ],
                   meta = tcga_aging$meta[tcga_selection1, "Age_Class"])
tcga_selection2 <- tcga_aging$meta %>% filter(Age_Class2 != "ND") %>% rownames()
tcga_path2 <- list(pathway = tcga_aging$pathway[tcga_selection2, ],
                   meta = tcga_aging$meta[tcga_selection2, "Age_Class2"])
cat("  Scheme 1 n:", length(tcga_selection1), "  Scheme 2 n:", length(tcga_selection2), "\n")

# ============================================================
# Assemble MICA input data
# ============================================================
cat("\n--- Assembling MICA input ---\n")

# Scheme 1: Young / Middle-Aged / Elderly
aging_selected_path1_data <- list(tcga = t(tcga_path1$pathway),
                                  metabric = t(metabric_path1$pathway),
                                  scanb = t(scanb_path1$pathway))
aging_selected_path1_meta <- list(tcga = tcga_path1$meta,
                                  metabric = metabric_path1$meta,
                                  scanb = scanb_path1$meta)

# Scheme 2: Early / Middle / Elderly
aging_selected_path2_data <- list(tcga = t(tcga_path2$pathway),
                                  metabric = t(metabric_path2$pathway),
                                  scanb = t(scanb_path2$pathway))
aging_selected_path2_meta <- list(tcga = tcga_path2$meta,
                                  metabric = metabric_path2$meta,
                                  scanb = scanb_path2$meta)

# Preflight: sample size check
cat("\nSample sizes per group per cohort:\n")
for (scheme_name in c("Scheme1", "Scheme2")) {
    meta_list <- if (scheme_name == "Scheme1") aging_selected_path1_meta else aging_selected_path2_meta
    cat(sprintf("\n  %s:\n", scheme_name))
    for (cohort in names(meta_list)) {
        tbl <- table(meta_list[[cohort]])
        cat(sprintf("    %s: %s\n", cohort, paste(names(tbl), tbl, sep = "=", collapse = ", ")))
        if (any(tbl < 10)) {
            cat(sprintf("    WARNING: %s has group(s) with n < 10!\n", cohort))
        }
    }
}

# Save GSVA outputs
save(aging_selected_path1_data, aging_selected_path1_meta,
     file = file.path(OUTPUT_DIR, "04_selected_path_aging_standard.RData"))
save(aging_selected_path2_data, aging_selected_path2_meta,
     file = file.path(OUTPUT_DIR, "04_selected_path_aging_postM.RData"))

cat("\nStage 1 complete. Output saved to", OUTPUT_DIR, "\n")
