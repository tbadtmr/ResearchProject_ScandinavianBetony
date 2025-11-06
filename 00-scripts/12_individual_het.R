#!/usr/bin/env Rscript
# ==============================================================
# Script: 12_individual_het.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Visualize observed heterozygosity (Ho) per individual and test
#   for group differences across all groups (no Scania split).
#   Performs Kruskal–Wallis and pairwise Wilcoxon tests (Holm).
#
# Inputs:
#   - 04-results/diversity_individual_allGroups.tsv
#
# Outputs:
#   - 04-results/heterozygosity_groups.svg
#   - 04-results/heterozygosity_kruskal.txt
#   - 04-results/heterozygosity_pairwise_wilcox.tsv
#
# Usage:
#   Rscript 00-scripts/12_individual_het.R
# ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(svglite)
})

# Inputs & outputs
in_tsv <- "04-results/diversity_individual_allGroups.tsv"
out_svg <- "04-results/heterozygosity_groups.svg"
out_kw  <- "04-results/heterozygosity_kruskal.txt"
out_pw  <- "04-results/heterozygosity_pairwise_wilcox.tsv"

if (!file.exists(in_tsv)) {
  stop("Input file not found: ", in_tsv,
       "\nRun 12_diversity_stats.sh first to generate it.")
}

# Read data
dat <- read.delim(in_tsv, check.names = FALSE) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Ho = as.numeric(Ho),
    Group = as.character(Group)
  )

# Rename & order groups
group_pretty <- c(
  "Western_Europe"    = "Western European",
  "Central_Europe"    = "Central European",
  "Balkan-Carpathian" = "Balkan–Carpathian",
  "Baltic"            = "Baltic",
  "Scandinavia"       = "Scandinavian"
)
group_order <- c("Western European", "Central European",
                 "Balkan–Carpathian", "Baltic", "Scandinavian")

dat <- dat %>%
  mutate(
    Group = recode(Group, !!!group_pretty, .default = Group),
    Group = factor(Group, levels = group_order)
  )

# Colour palette
group_colors <- c(
  "Western European"   = "#f781bf",
  "Central European"   = "#a65628",
  "Balkan–Carpathian"  = "#66a61e",
  "Baltic"             = "#a6cee3",
  "Scandinavian"       = "purple"
)

# Plot
set.seed(123)
pos <- position_jitter(width = 0.2, height = 0)

p <- ggplot(dat, aes(x = Group, y = Ho)) +
  geom_boxplot(fill = "grey85", color = "grey50", alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(color = Group), position = pos, size = 2.8, alpha = 0.9) +
  scale_color_manual(values = group_colors, guide = "none") +
  labs(x = NULL, y = "Observed heterozygosity (Ho)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(out_svg, p, width = 7.2, height = 5.0, dpi = 300)
message("Saved plot: ", out_svg)

# Statistical tests
kw <- kruskal.test(Ho ~ Group, data = dat)
sink(out_kw); cat("Kruskal–Wallis test for Ho by group\n\n"); print(kw); sink()
message("Saved KW result: ", out_kw)

pw <- pairwise.wilcox.test(dat$Ho, dat$Group, p.adjust.method = "holm", exact = FALSE)
pw_tab <- as.data.frame(as.table(pw$p.value)) %>%
  rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
  filter(!is.na(p_adj))
write.table(pw_tab, out_pw, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved pairwise Wilcoxon (Holm): ", out_pw)
