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
