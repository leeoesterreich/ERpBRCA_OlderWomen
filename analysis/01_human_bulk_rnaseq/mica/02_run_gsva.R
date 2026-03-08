#!/usr/bin/env Rscript
# 02_run_gsva.R - Run GSVA on EstroGene and MSigDB pathways
#
# Inputs:
#   - outputs/01_*_filtered.RData
#   - data/input/EstroGene_Signatures.xlsx
#
# Outputs:
#   - outputs/02_gsva_estrogene.RData
#   - outputs/02_gsva_pathways.RData

set.seed(12345)

library(GSVA)
library(tidyverse)

# TODO: Implement once data files are obtained

message("02_run_gsva.R: Waiting for data files from Jian")
