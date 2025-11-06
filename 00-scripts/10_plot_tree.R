#!/usr/bin/env Rscript
# ==============================================================
# Script: 10_plot_tree.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Read the IQ-TREE output (consensus tree if available), attach
#   sample metadata, optionally re-root on an outgroup tip, and
#   export a publication-ready SVG to 04-results/.
#
# Inputs:
#   - Tree: 03-analysis/08-downstream/03-phylogeny/<RUN>/<RUN>.contree (preferred)
#            or 03-analysis/08-downstream/03-phylogeny/<RUN>/<RUN>.treefile
#   - Metadata: 01-info_files/samples_clean.txt  (expects: IID, Name, Group)
#
# Outputs:
#   - 04-results/<RUN>_tree.svg
#
# Usage:
#   Rscript 00-scripts/10_plot_tree.R r70_noTrentino_noOutgroup
#   # (optional re-root if needed; pass exact tip label as in the tree)
#   Rscript 00-scripts/10_plot_tree.R r70_noTrentino_noOutgroup "49_Trento_B_hirsuta.1"
# ==============================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)    
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(grid)     
})

# USER ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 00-scripts/11_plot_tree.R <RUN> [OUTGROUP_TIP]", call. = FALSE)
}
RUN <- args[[1]]
OUTGROUP_TIP <- if (length(args) >= 2) args[[2]] else NA_character_

# Paths
phy_dir <- file.path("03-analysis", "08-downstream", "03-phylogeny", RUN)
contree_fp  <- file.path(phy_dir, paste0(RUN, ".contree"))
treefile_fp <- file.path(phy_dir, paste0(RUN, ".treefile"))
tree_path <- if (file.exists(contree_fp)) contree_fp else treefile_fp

meta_fp <- file.path("01-info_files", "samples_clean.txt")
out_svg <- file.path("04-results", paste0(RUN, "_tree.svg"))
dir.create("04-results", showWarnings = FALSE)

# Read sample metadata
# Expect columns: IID, Name, Group
meta <- read_delim(meta_fp, delim = "\t", show_col_types = FALSE) %>%
  mutate(across(everything(), ~ trimws(as.character(.))))

# Read tree
tr <- read.tree(tree_path)
# dequote labels if quoted
tr$tip.label <- gsub("^'(.*)'$", "\\1", trimws(tr$tip.label))

# Optional: reroot on outgroup tip
if (!is.na(OUTGROUP_TIP)) {
  if (!requireNamespace("phytools", quietly = TRUE)) {
    install.packages("phytools", quiet = TRUE)
  }
  tr <- phytools::reroot(
    tr,
    node.number = which(tr$tip.label == OUTGROUP_TIP),
    position = 1e-6,
    resolve.root = TRUE
  )
}

# tidy left-to-right
tr <- ladderize(tr, right = FALSE)

# light branch rescale for readability (optional)
if (!is.null(tr$edge.length)) tr$edge.length <- log1p(tr$edge.length)

# add a small basal stem for aesthetics
max_bl <- if (length(tr$edge.length)) max(tr$edge.length, na.rm = TRUE) else 0.1
tr$root.edge <- 0.08 * max_bl

# Relabel tips IID -> Name
tip_df <- tibble(IID = tr$tip.label) %>%
  left_join(meta, by = "IID") %>%
  mutate(Name_lab = ifelse(is.na(Name), IID, Name),
         Name_clean = gsub("_", " ", Name_lab))

# keep original order of tips
tr$tip.label <- tip_df$Name_lab

# Build plot & attach metadata
p <- ggtree(tr, layout = "rectangular") +
  geom_rootedge(linewidth = 0.7)

td <- p$data
max_x <- max(td$x, na.rm = TRUE)

# attach per-tip metadata by current label (Name_lab)
tip_meta <- tibble(label = tr$tip.label) %>%
  left_join(tip_df %>% select(Name_lab, Name_clean, IID, Group),
            by = c("label" = "Name_lab")) %>%
  mutate(Group = ifelse(is.na(Group), "Unknown", Group),
         Group = factor(Group))

p <- p %<+% tip_meta

# color mapping (same palette family you’ve used elsewhere)
group_colors <- c(
  "Scandinavia"        = "purple",
  "Baltic"             = "#a6cee3",
  "Balkan-Carpathian"  = "#66a61e",
  "France"             = "#e41a1c",
  "Central_Europe"     = "#a65628",
  "Western_Europe"     = "#f781bf",
  "Outgroup"           = "orange",
  "Alpine"             = "yellow",
  "Unknown"            = "grey50"
)

# if OUTGROUP_TIP provided, tag that sample as Outgroup in the plot
if (!is.na(OUTGROUP_TIP)) {
  p$data$Group[p$data$IID == OUTGROUP_TIP] <- "Outgroup"
}

# internal node labels (bootstrap) >= 50 if present (contree usually has them)
internal_nodes <- p$data %>%
  filter(!isTip & !is.na(label) & label != "") %>%
  mutate(label_num = suppressWarnings(as.numeric(label))) %>%
  filter(!is.na(label_num), label_num >= 50) %>%
  mutate(x_mid = x - branch.length/2, y_shift = y + 0.5)

p <- p +
  geom_tiplab(aes(label = Name_clean), size = 3, offset = 0.0025) +
  geom_text(data = internal_nodes,
            aes(x = x_mid - 0.001, y = y_shift, label = label),
            size = 2.5) +
  geom_tippoint(aes(color = Group), size = 2.6, show.legend = FALSE) +
  scale_color_manual(values = group_colors, guide = "none") +
  theme_tree2() +
  theme(
    axis.line   = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    plot.margin = unit(c(0.5, 0.6, 0.5, 0.8), "cm")
  ) +
  coord_cartesian(xlim = c(0, max_x + 0.04))

ggsave(out_svg, p, width = 8, height = 10, units = "in")
message("Saved: ", out_svg)
