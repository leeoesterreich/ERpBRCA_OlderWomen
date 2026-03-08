# MICA Integration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Integrate Jian's MICA analysis code into ERpBRCA_OlderWomen to reproduce manuscript Figures 4C, 4H, 5C, 5E.

**Architecture:** Modular R script pipeline in `analysis/01_human_bulk_rnaseq/mica/` with numbered scripts following project conventions. Data prep → GSVA → MICA → Figure generation.

**Tech Stack:** R 4.3.3, MICA package (jianzou75/MICA), GSVA, ComplexHeatmap, ggplot2, patchwork, edgeR

**Blocker:** Data files needed from Jian (see `Jian_MICA/required_data_files.md`). Tasks 1-3 can proceed now; Tasks 4+ require data.

---

## Phase 1: Directory Structure & Scaffolding (No Data Required)

### Task 1: Create Directory Structure

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/mica/`
- Create: `analysis/01_human_bulk_rnaseq/mica/data/input/.gitkeep`
- Create: `analysis/01_human_bulk_rnaseq/mica/outputs/.gitkeep`
- Create: `analysis/01_human_bulk_rnaseq/mica/figures/.gitkeep`
- Create: `analysis/01_human_bulk_rnaseq/mica/logs/.gitkeep`

**Step 1: Create directories**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
mkdir -p analysis/01_human_bulk_rnaseq/mica/{data/input,outputs,figures,logs}
touch analysis/01_human_bulk_rnaseq/mica/{data/input,outputs,figures,logs}/.gitkeep
```

**Step 2: Verify structure**

```bash
ls -la analysis/01_human_bulk_rnaseq/mica/
```
Expected: data/, outputs/, figures/, logs/ directories

**Step 3: Commit**

```bash
git add analysis/01_human_bulk_rnaseq/mica/
git commit -m "feat(mica): create directory structure for MICA integration"
```

---

### Task 2: Create SLURM Batch Script

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/mica/run_all.sbatch`

**Step 1: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=mica_analysis
#SBATCH --partition=htc
#SBATCH --time=4:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=10
#SBATCH --output=logs/mica_%j.out
#SBATCH --error=logs/mica_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

# Activate environment
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
cd "$SCRIPT_DIR"

echo "=== MICA Analysis Pipeline ==="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"

# Run pipeline sequentially
echo "Step 1: Preparing data..."
Rscript 01_prep_data.R

echo "Step 2: Running GSVA..."
Rscript 02_run_gsva.R

echo "Step 3: Running MICA..."
Rscript 03_run_mica.R

echo "Step 4: Generating Figure 4..."
Rscript 04_fig4_estrogen.R

echo "Step 5: Generating Figure 5..."
Rscript 05_fig5_immune.R

echo "=== Pipeline complete ==="
echo "End time: $(date)"
```

**Step 2: Make executable and commit**

```bash
chmod +x analysis/01_human_bulk_rnaseq/mica/run_all.sbatch
git add analysis/01_human_bulk_rnaseq/mica/run_all.sbatch
git commit -m "feat(mica): add SLURM batch script for pipeline execution"
```

---

### Task 3: Create Script Stubs

**Files:**
- Create: `analysis/01_human_bulk_rnaseq/mica/01_prep_data.R`
- Create: `analysis/01_human_bulk_rnaseq/mica/02_run_gsva.R`
- Create: `analysis/01_human_bulk_rnaseq/mica/03_run_mica.R`
- Create: `analysis/01_human_bulk_rnaseq/mica/04_fig4_estrogen.R`
- Create: `analysis/01_human_bulk_rnaseq/mica/05_fig5_immune.R`

**Step 1: Create 01_prep_data.R stub**

```r
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
```

**Step 2: Create 02_run_gsva.R stub**

```r
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
```

**Step 3: Create 03_run_mica.R stub**

```r
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
```

**Step 4: Create 04_fig4_estrogen.R stub**

```r
#!/usr/bin/env Rscript
# 04_fig4_estrogen.R - Generate Figure 4 panels (estrogen/HSD17B)
#
# Figure 4C: HSD17B gene heatmap across age groups
# Figure 4H: EstroGene pathway boxplots
#
# Inputs:
#   - outputs/01_*_filtered.RData (for gene expression)
#   - outputs/02_gsva_estrogene.RData (for pathway scores)
#
# Outputs:
#   - figures/fig4c_hsd17b_heatmap.{pdf,png}
#   - figures/fig4h_estrogen_boxplots.{pdf,png}

set.seed(12345)

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(tidyverse)

# HSD17B genes for Figure 4C
HSD17B_GENES <- c("ESR1", "GREB1", "PGR", "SAA1", "RAB19", "KRT37", "TRPM8",
                  "CYP19A1", "HSD17B1", "HSD17B7", "HSD17B12", "HSD17B2",
                  "HSD17B10", "HSD17B14")

# Age group colors
AGE_COLORS <- c("Young" = "#F94040", "Middle-Aged" = "#5757F9", "Elderly" = "#610051")

# Heatmap function (from Jian's code)
heatmap_median <- function(study.data.list, study.label.list) {
  median_df <- list()
  for(i in seq_along(study.data.list)) {
    median_df[[i]] <- aggregate(x = study.data.list[[i]],
                                by = list(study.label.list[[i]]),
                                FUN = median) %>%
      column_to_rownames("Group.1")
  }
  median_df <- lapply(median_df, t)
  median_df_merge <- do.call("cbind", median_df)
  data_source <- rep(names(study.data.list), each = ncol(median_df[[1]]))

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

# Boxplot function (from Jian's code)
path_boxplot <- function(data_list, meta_list, path, age_levels = c("Young", "Middle-Aged", "Elderly")) {
  study.tbl <- data.frame(freq = sapply(meta_list, length))
  expr <- data.frame(
    expr = unlist(sapply(data_list, function(x) scale(x[, path]))),
    age = unlist(sapply(meta_list, function(x) factor(x, levels = age_levels))),
    study = rep(rownames(study.tbl), study.tbl$freq)
  )

  ggplot(expr, aes(x = age, y = expr, fill = age)) +
    geom_boxplot(aes(middle = mean(expr))) +
    facet_wrap(~study) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text.x = element_text(face = "bold", size = 20),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    scale_fill_manual(values = AGE_COLORS) +
    xlab("") + ylab("") +
    ylim(c(-2, 2)) +
    ggtitle(path)
}

# TODO: Implement once data is available

message("04_fig4_estrogen.R: Waiting for processed data")
```

**Step 5: Create 05_fig5_immune.R stub**

```r
#!/usr/bin/env Rscript
# 05_fig5_immune.R - Generate Figure 5 panels (inflammation/immune)
#
# Figure 5C: Inflammatory pathway heatmap
# Figure 5E: Immune cell boxplots (CD8, cytotoxic, monocyte, DC)
#
# Inputs:
#   - outputs/02_gsva_pathways.RData
#   - outputs/03_mica_results.RData (optional, for significant pathways)
#
# Outputs:
#   - figures/fig5c_inflammatory_heatmap.{pdf,png}
#   - figures/fig5e_immune_boxplots.{pdf,png}

set.seed(12345)

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(tidyverse)

# Immune cell types for Figure 5E
IMMUNE_CELLS <- c("CD8 T cells", "Cytotoxic lymphocytes",
                  "Monocytic lineage", "Myeloid dendritic cells")

# Age group colors
AGE_COLORS <- c("Young" = "#F94040", "Middle-Aged" = "#5757F9", "Elderly" = "#610051")

# Cell boxplot function (from Jian's code)
cell_boxplot <- function(mcp_list, meta_list, cell, age_levels = c("Young", "Middle-Aged", "Elderly")) {
  study.tbl <- data.frame(freq = sapply(meta_list, length))
  expr <- data.frame(
    expr = unlist(sapply(mcp_list, function(x) scale(x[, cell]))),
    age = unlist(sapply(meta_list, function(x) factor(x, levels = age_levels))),
    study = rep(rownames(study.tbl), study.tbl$freq)
  )

  ggplot(expr, aes(x = age, y = expr, fill = age)) +
    geom_boxplot(aes(middle = mean(expr))) +
    facet_wrap(~study) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text.x = element_text(face = "bold", size = 15),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 8)
    ) +
    scale_fill_manual(values = AGE_COLORS) +
    xlab("") + ylab("") +
    ylim(c(-1, 3)) +
    ggtitle(cell)
}

# TODO: Implement once data is available

message("05_fig5_immune.R: Waiting for processed data")
```

**Step 6: Commit all stubs**

```bash
git add analysis/01_human_bulk_rnaseq/mica/*.R
git commit -m "feat(mica): add script stubs with reusable functions from Jian's code"
```

---

## Phase 2: Data Integration (Requires Data From Jian)

### Task 4: Link/Copy Data Files

**Prerequisite:** Obtain data files from Jian (see `Jian_MICA/required_data_files.md`)

**Files:**
- Symlink/copy to: `analysis/01_human_bulk_rnaseq/mica/data/input/`

**Step 1: Create symlinks to data files**

```bash
cd analysis/01_human_bulk_rnaseq/mica/data/input/
# Example (paths will depend on where Jian provides data):
ln -s /path/to/ilc_clean_data.RData .
ln -s /path/to/EstroGene_Signatures.xlsx .
ln -s /path/to/METABRIC_Clinical_Info.csv .
ln -s /path/to/TCGA_Key.csv .
```

**Step 2: Verify data accessibility**

```bash
ls -la data/input/
```
Expected: All required files present (or symlinks)

---

### Task 5: Implement 01_prep_data.R

**Files:**
- Modify: `analysis/01_human_bulk_rnaseq/mica/01_prep_data.R`

**Step 1: Implement full data loading and filtering**

Replace TODO section with actual implementation based on Jian's `04_EstroGene_MICA[1].R`:

```r
# Load BCdata objects
load("data/input/ilc_clean_data.RData")

# Load clinical data
metabric_meta <- read.csv("data/input/METABRIC_Clinical_Info.csv", row.names = 1)
tcga_meta <- read.csv("data/input/TCGA_Key.csv", row.names = 1)

# Process SCAN-B
scanb_meta <- scanb@clinic %>%
  mutate(Age_Class = create_age_groups(age)) %>%
  mutate(Age_Class_PostM = create_age_groups_postm(age)) %>%
  filter(clinic_grp == "ERpHER2n")

# Process METABRIC
metabric_meta <- metabric_meta %>%
  mutate(Age_Class = create_age_groups(age_at_diagnosis)) %>%
  mutate(Age_Class_PostM = create_age_groups_postm(age_at_diagnosis)) %>%
  filter(ER == "Positive" & HER2_IHC_status %in% c(0, 1))

# Process TCGA
tcga_meta <- tcga_meta %>%
  mutate(Age_Class = create_age_groups(Age)) %>%
  mutate(Age_Class_PostM = create_age_groups_postm(Age)) %>%
  filter(ER == "Positive" & HER2 == "Negative")

# Save filtered data
save(scanb, scanb_meta, file = "outputs/01_scanb_filtered.RData")
save(metabric, metabric_meta, file = "outputs/01_metabric_filtered.RData")
save(tcga, tcga_meta, file = "outputs/01_tcga_filtered.RData")

message("01_prep_data.R: Complete")
```

**Step 2: Test locally (if data available)**

```bash
Rscript analysis/01_human_bulk_rnaseq/mica/01_prep_data.R
```

**Step 3: Commit**

```bash
git add analysis/01_human_bulk_rnaseq/mica/01_prep_data.R
git commit -m "feat(mica): implement data preparation script"
```

---

### Task 6: Implement 02_run_gsva.R

Similar pattern: implement full GSVA logic from Jian's code.

---

### Task 7: Implement 03_run_mica.R

Similar pattern: implement MICA analysis from Jian's code.

---

### Task 8: Implement 04_fig4_estrogen.R

Similar pattern: implement Figure 4 generation.

---

### Task 9: Implement 05_fig5_immune.R

Similar pattern: implement Figure 5 generation.

---

## Phase 3: Validation

### Task 10: Run Full Pipeline

**Step 1: Submit SLURM job**

```bash
cd analysis/01_human_bulk_rnaseq/mica/
sbatch run_all.sbatch
```

**Step 2: Monitor progress**

```bash
squeue -u $USER
tail -f logs/mica_*.out
```

**Step 3: Verify outputs**

```bash
ls -la outputs/
ls -la figures/
```
Expected: All RData files and PDF/PNG figures present

---

### Task 11: Validate Figures Against Manuscript

**Step 1: Visual comparison**

Compare generated figures with manuscript:
- `figures/fig4c_hsd17b_heatmap.png` vs Manuscript Figure 4C
- `figures/fig4h_estrogen_boxplots.png` vs Manuscript Figure 4H
- `figures/fig5c_inflammatory_heatmap.png` vs Manuscript Figure 5C
- `figures/fig5e_immune_boxplots.png` vs Manuscript Figure 5E

**Step 2: Document any discrepancies**

If figures differ, document in validation report.

**Step 3: Final commit**

```bash
git add -A
git commit -m "feat(mica): complete MICA integration with all figures"
```

---

## Summary

| Phase | Tasks | Status |
|-------|-------|--------|
| Phase 1: Scaffolding | 1-3 | Ready to implement |
| Phase 2: Data Integration | 4-9 | Blocked on data from Jian |
| Phase 3: Validation | 10-11 | After Phase 2 |

**Next Step:** Implement Tasks 1-3 (no data required), then wait for Jian to provide data files.
