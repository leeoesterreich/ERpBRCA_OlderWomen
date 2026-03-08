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
