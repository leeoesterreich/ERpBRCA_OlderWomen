#!/usr/bin/env Rscript
# 01_prep_data.R - Load and filter TCGA/METABRIC/SCAN-B data for MICA
#
# Inputs:
#   - data/input/ilc_clean_data.RData
#   - data/input/METABRIC_Clinical_Info.csv
#   - data/input/TCGA_Key.csv
#
# Outputs:
#   - outputs/01_tcga_filtered.RData
#   - outputs/01_metabric_filtered.RData
#   - outputs/01_scanb_filtered.RData

set.seed(12345)

library(tidyverse)
library(edgeR)

# Age group definitions
create_age_groups <- function(age) {
  case_when(
    age >= 35 & age <= 45 ~ "Young",
    age >= 55 & age <= 69 ~ "Middle-Aged",
    age >= 70 ~ "Elderly",
    TRUE ~ NA_character_
  )
}

create_age_groups_postm <- function(age) {
  case_when(
    age >= 55 & age < 60 ~ "Early",
    age >= 60 & age < 70 ~ "Middle",
    age >= 70 ~ "Elderly",
    TRUE ~ NA_character_
  )
}

# TODO: Implement once data files are obtained
# See Jian_MICA/required_data_files.md for required inputs

message("01_prep_data.R: Waiting for data files from Jian")
message("Required files listed in: Jian_MICA/required_data_files.md")
