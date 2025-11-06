#!/usr/bin/env Rscript
# ==============================================================
# Script: 09_plot_admixture.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Produces two visualizations of ADMIXTURE results:
#   (1) Cross-validation error plot (mean ± SD across replicates)
#   (2) Faceted admixture barplot for selected K values (here K = 2–4),
#       showing ancestry proportions per individual grouped by population.
# ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)   
  library(ggrepel)     
  library(svglite)    
})

# USER INPUT: which run to plot
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 00-scripts/12_plot_admixture.R <RUN>\n",
       "Example: Rscript 00-scripts/12_plot_admixture.R r70_noTrentino_noOutgroup", call. = FALSE)
}
RUN <- args[[1]]

# FILE PATHS
# All ADMIXTURE output lives under 03-analysis/08-downstream/02-admixture/<RUN>/
admix_dir <- file.path("03-analysis", "08-downstream", "02-admixture", RUN)
fam_fp   <- file.path(admix_dir, paste0(RUN, "_pruned.fam"))     # sample IDs and order
cv_all   <- file.path(admix_dir, "cv_all.tsv")                   # CV error of each replicate
cv_sum   <- file.path(admix_dir, "cv_summary.tsv")               # summarized CV mean ± SD

# Sample metadata (used to color / group bars)
meta_fp  <- file.path("01-info_files", "samples_clean.txt")
if (!file.exists(meta_fp)) {
  alt <- file.path("01-infofiles", "samples_clean.txt")          
  if (file.exists(alt)) meta_fp <- alt
}

# Output figure paths
dir.create("04-results", showWarnings = FALSE)
out_cv   <- file.path("04-results", paste0(RUN, "_admixture_cv.svg"))
out_qsvg <- file.path("04-results", paste0(RUN, "_admixture_K2-4.svg"))

# AESTHETIC SETTINGS
# Which K values to visualize (matches number of ancestral clusters)
Ks_to_plot <- c(2, 3, 4)

# Order and colors of population groups on the x-axis
group_order <- c("Western_Europe","Central_Europe","Balkan-Carpathian","Baltic","Scandinavia")
group_colors <- c(
  "Scandinavia"        = "purple",
  "Baltic"             = "#a6cee3",
  "Balkan-Carpathian"  = "#66a61e",
  "France"             = "#e41a1c",
  "Central_Europe"     = "#a65628",
  "Western_Europe"     = "#f781bf"
)
group_pretty <- c(
  "Western_Europe"    = "Western Europe",
  "Central_Europe"    = "Central Europe",
  "Balkan-Carpathian" = "Balkan–Carpathian",
  "Baltic"            = "Baltic",
  "Scandinavia"       = "Scandinavia"
)


# 1. CROSS-VALIDATION ERROR PLOT
# --------------------------------------------------------------
# Shows how ADMIXTURE model fit (CV error) changes with K.
# Lower CV error = better fit (but beware overfitting).


cvsum <- read_tsv(cv_sum, show_col_types = FALSE) %>%
  arrange(K)

# Identify which K has the lowest mean CV error
bestK <- cvsum$K[which.min(cvsum$mean)]

# Simple line + error bar plot
p_cv <- ggplot(cvsum, aes(K, mean)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15) +
  geom_point(data = subset(cvsum, K == bestK), color = "red", size = 2.5) +  # highlight best K
  labs(x = "K (number of ancestral clusters)", y = "Cross-validation error") +
  theme_classic()

ggsave(out_cv, p_cv, width = 4.5, height = 3.2, dpi = 300)
message("✓ Saved CV plot: ", out_cv)


# 2. BUILD SAMPLE ORDER / GROUPS
# --------------------------------------------------------------
# Combines .fam file (PLINK sample order) with sample metadata.
# Used to:
#   - keep samples in the same order across Ks,
#   - assign population groups for faceting.


fam <- read_table(fam_fp, col_names = c("FID","IID","V3","V4","V5","V6"), show_col_types = FALSE)
meta <- read_delim(meta_fp, delim = "\t", show_col_types = FALSE)

# Merge IDs and metadata; preserve desired group order
order_df <- fam %>%
  select(IID) %>%
  left_join(meta, by = "IID") %>%
  mutate(
    # Add suffix if two individuals share the same name
    Name_lab = if_else(duplicated(Name) | duplicated(Name, fromLast = TRUE),
                       paste0(Name, " (", IID, ")"), Name),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(Group, Name_lab) %>%
  mutate(Name_ord = factor(Name_lab, levels = unique(Name_lab)))

name_levels  <- levels(order_df$Name_ord)
group_levels <- levels(order_df$Group)


# 3. CHOOSE BEST REPLICATE PER K
# --------------------------------------------------------------
# Reads cv_all.tsv (one line per replicate), finds the replicate
# with the smallest CV error for each K, and uses that .Q file
# for plotting.

cvall <- read_tsv(cv_all, show_col_types = FALSE) # columns: K, replicate, cv_error
best_rep <- cvall %>%
  group_by(K) %>%
  filter(cv_error == min(cv_error, na.rm = TRUE)) %>%
  slice(1) %>%                   # keep first if ties
  ungroup() %>%
  filter(K %in% Ks_to_plot) %>%
  select(K, replicate)


# 4. READ BEST Q FILES AND STACK INTO ONE LONG TABLE
# --------------------------------------------------------------
# Each .Q file is a matrix: rows = individuals, columns = clusters.
# After adding sample info, pivot to long format:
# IID, Cluster, Proportion, K, Group, Name_ord


read_Q_best <- function(K, rep_id) {
  qfile <- file.path(admix_dir, sprintf("%s_pruned.%d.rep%d.Q", RUN, K, rep_id))
  q <- read_table(qfile, col_names = FALSE, show_col_types = FALSE)
  colnames(q) <- paste0("Cluster", seq_len(K))
  q$IID <- fam$IID
  q %>%
    inner_join(order_df %>% select(IID, Name_lab, Name_ord, Group), by = "IID") %>%
    pivot_longer(starts_with("Cluster"), names_to = "Cluster", values_to = "Proportion") %>%
    mutate(
      K = factor(paste0("K = ", K), levels = paste0("K = ", Ks_to_plot)),
      Name_ord = factor(Name_ord, levels = name_levels),
      Group    = factor(Group, levels = group_levels)
    )
}

long_all <- map2_dfr(best_rep$K, best_rep$replicate, read_Q_best) %>%
  arrange(K, Group, Name_ord)


# 5. DEFINE COLOR PALETTE AND PLOT STYLE
# --------------------------------------------------------------
# Assign distinct colors to ancestry clusters.
# Set up a minimal clean theme and facet layout.


cluster_levels <- long_all %>% pull(Cluster) %>% unique() %>% sort()
cluster_palette <- setNames(
  RColorBrewer::brewer.pal(max(3, length(cluster_levels)), "Set2")[seq_along(cluster_levels)],
  cluster_levels
)

group_labeller <- labeller(Group = as_labeller(group_pretty))

admix_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.placement    = "outside",
    strip.background   = element_rect(fill = "white", color = NA),
    strip.text.x       = element_text(size = 8, face = "bold"),
    strip.text.y.right = element_text(size = 8, face = "bold"),
    axis.text.x        = element_text(angle = 60, hjust = 1, vjust = 1, size = 6)
  )

# 6. DRAW THE ADMIXTURE BARPLOT
# --------------------------------------------------------------
# Each bar = one individual.
# Height segments = ancestry proportions (sum to 1 per sample).
# Facets by K (rows) and population group (columns).

p_k234 <- ggplot(long_all, aes(x = Name_ord, y = Proportion, fill = Cluster)) +
  geom_col(width = 1, colour = NA) +
  geom_vline(xintercept = seq(1.5, length(levels(long_all$Name_ord)) - 0.5, by = 1),
             linewidth = 0.25, colour = scales::alpha("grey20", 0.12)) +
  facet_grid(rows = vars(K), cols = vars(Group),
             scales = "free_x", space = "free_x", switch = "y",
             labeller = group_labeller) +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  scale_x_discrete(expand = c(0, 0), labels = function(x) gsub("_", " ", x)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02)),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Ancestry proportion") +
  admix_theme

ggsave(out_qsvg, p_k234, width = 8.5, height = 5.0, dpi = 300)
message("Saved admixture K=2–4 plot: ", out_qsvg)
