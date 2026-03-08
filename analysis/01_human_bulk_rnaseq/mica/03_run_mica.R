#!/usr/bin/env Rscript
# 03_run_mica.R - Run MICA concordance analysis
#
# Inputs:
#   - outputs/02_gsva_*.RData
#
# Outputs:
#   - outputs/03_mica_results.RData
#   - outputs/03_mica_summary.csv

set.seed(12345)

library(MICA)
library(parallel)
library(tidyverse)

# MICA parameters
N_PERM <- 500
P_THRESHOLD <- 0.05
N_PARALLEL <- 10

# TODO: Implement once GSVA outputs are available

message("03_run_mica.R: Waiting for GSVA outputs")
