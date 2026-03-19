#!/usr/bin/env Rscript
# 02_run_gsva.R — GSVA scoring is now integrated into 01_prep_data.R
#
# This script exists for backward compatibility with run_all.sbatch.
# GSVA computation on EstroGene + LI_ESTROGENE signatures is performed
# in 01_prep_data.R which outputs:
#   - outputs/04_selected_path_aging_standard.RData  (Scheme 1)
#   - outputs/04_selected_path_aging_postM.RData     (Scheme 2)
#
# Nothing to do here.

set.seed(12345)

cat("02_run_gsva.R: GSVA scoring is integrated into 01_prep_data.R — skipping.\n")
