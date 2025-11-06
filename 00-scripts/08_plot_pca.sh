#!/usr/bin/env Rscript
# ==============================================================
# Script: 08_plot_pca.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Plot PCA (PC1–PC3) from PLINK outputs and annotate with metadata.
#   Reads:
#     03-analysis/08-downstream/01-pca/<RUN>_pca.eigenvec
#     03-analysis/08-downstream/01-pca/<RUN>_pca.eigenval
#     01-info_files/samples_clean.txt   (fallback: 01-infofiles/samples_clean.txt)
#   Writes:
#     04-results/<RUN>_pca.svg
#
# Usage:
#   Rscript 00-scripts/08_plot_pca.R r70_noTrentino_noOutgroup
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 00-scripts/10_plot_pca.R <RUN>\nExample: Rscript 00-scripts/10_plot_pca.R r70_noTrentino_noOutgroup", call. = FALSE)
}
RUN <- args[[1]]

# paths to input files
pca_dir   <- file.path("03-analysis", "08-downstream", "01-pca")
eigvec_fp <- file.path(pca_dir, paste0(RUN, "_pca.eigenvec"))
eigval_fp <- file.path(pca_dir, paste0(RUN, "_pca.eigenval"))

meta_fp   <- file.path("01-info_files", "samples_clean.txt")
if (!file.exists(meta_fp)) {
  alt <- file.path("01-infofiles", "samples_clean.txt")
  if (file.exists(alt)) meta_fp <- alt
}

# ensure results folder exists
dir.create("04-results", showWarnings = FALSE)
out_svg <- file.path("04-results", paste0(RUN, "_pca.svg"))

# colors and labels for geographic groups
group_colors <- c(
  "Scandinavia"        = "purple",
  "Baltic"             = "#a6cee3",
  "Balkan-Carpathian"  = "#66a61e",
  "France"             = "#e41a1c",
  "Central_Europe"     = "#a65628",
  "Outgroup"           = "orange",
  "Alpine"             = "yellow",
  "Western_Europe"     = "#f781bf"
)
group_labels <- setNames(gsub("_", " ", names(group_colors)), names(group_colors))

# read PLINK outputs
pca <- read_tsv(eigvec_fp, col_names = TRUE, show_col_types = FALSE)
if (!all(c("FID","IID") %in% names(pca))) {
  nPC <- ncol(pca) - 2
  names(pca) <- c("FID","IID", paste0("PC", seq_len(nPC)))
}
pca$IID <- sub("\\.1$", "", pca$IID)

eigenval <- scan(eigval_fp, quiet = TRUE)
variance <- eigenval / sum(eigenval) * 100

# metadata
meta <- read_delim(meta_fp, delim = "\t", show_col_types = FALSE)

# merge + label construction
pca_annot <- pca %>%
  left_join(meta, by = "IID") %>%
  mutate(
    Label = ifelse(!is.na(Name) & Name != "", Name, IID),
    num       = str_extract(Name, "^[0-9]+"),
    cc        = str_extract(Name, "[A-Za-z]{2}$"),
    locality  = ifelse(is.na(Name), NA_character_,
                       Name %>% sub("^[0-9]+_", "", .) %>% sub("_[A-Za-z]{2}$", "", .)),
    num_clean = ifelse(is.na(num), NA_character_, as.character(as.integer(num))),
    Label = dplyr::case_when(
      cc == "SE" & !is.na(num_clean) ~ paste(num_clean, gsub("_", " ", locality)),
      !is.na(num_clean)              ~ num_clean,
      TRUE                           ~ IID
    )
  )

pca_annot$Group <- factor(pca_annot$Group, levels = names(group_colors))

pca_scale <- scale_color_manual(values = group_colors,
                                breaks = names(group_colors),
                                labels = group_labels)

pca_theme <- theme_classic(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(color = "black"),
        axis.ticks       = element_line(color = "black"),
        legend.title     = element_blank(),
        legend.key       = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))

# PLOTS
# PCA1 vs PC2
pc1_pc2 <- ggplot(pca_annot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group), size = 4) +
  geom_text_repel(aes(label = Label), size = 3.2, max.overlaps = 60,
                  color = "black", show.legend = FALSE) +
  labs(x = paste0("PC1 (", round(variance[1], 1), "%)"),
       y = paste0("PC2 (", round(variance[2], 1), "%)")) +
  pca_scale + pca_theme
# PC2 vs PC3
pc2_pc3 <- ggplot(pca_annot, aes(x = PC2, y = PC3)) +
  geom_point(aes(color = Group), size = 4) +
  geom_text_repel(aes(label = Label), size = 3.2, max.overlaps = 60,
                  color = "black", show.legend = FALSE) +
  labs(x = paste0("PC2 (", round(variance[2], 1), "%)"),
       y = paste0("PC3 (", round(variance[3], 1), "%)")) +
  pca_scale + pca_theme
# PC1 vs PC3
pc3_pc1 <- ggplot(pca_annot, aes(x = PC3, y = PC1)) +
  geom_point(aes(color = Group), size = 5) +
  geom_text_repel(aes(label = Label), size = 3.0,
                  max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2,
                  color = "black", show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = paste0("PC3 (", round(variance[3], 1), "%)"),
       y = paste0("PC1 (", round(variance[1], 1), "%)")) +
  pca_scale + pca_theme

combined_row <- (pc1_pc2 + theme(legend.position = "none")) |
                (pc2_pc3 + theme(legend.position = "none")) |
                 pc3_pc1

# ---- save SVG ----
ggsave(out_svg, combined_row, width = 18, height = 6, dpi = 300)
message("Saved PCA plot: ", out_svg)
