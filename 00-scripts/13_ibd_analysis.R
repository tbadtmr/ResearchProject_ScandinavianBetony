#!/usr/bin/env Rscript
# ==============================================================
# Script: 13_ibd_analysis.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Analyze isolation by distance (IBD) using PLINK 1–IBS distances
#   and sample coordinates. Performs Mantel tests (9,999 perms)
#   and generates regression plots.
#
# Inputs:
#   - 03-analysis/08-downstream/05-ibd/ibs_<RUN>.mdist(.id)
#   - 01-info_files/samples_clean.txt
#
# Outputs:
#   - 04-results/ibd_mantel_summary.tsv
#   - 04-results/ibd_within_group_mantel.tsv
#   - 04-results/ibd_all_pairs.svg
#   - 04-results/ibd_excluding_western.svg
#   - 04-results/ibd_within_groups_facets.svg
#
# Usage:
#   Rscript 00-scripts/13_ibd_analysis.R r70_noTrentino_noOutgroup
# ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(geosphere)
  library(vegan)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 00-scripts/18_ibd_analysis.R <RUN>")
}
RUN <- args[[1]]

# Paths
ibs_dir  <- file.path("03-analysis","08-downstream","05-ibd")
mdist_fp <- file.path(ibs_dir, paste0("ibs_", RUN, ".mdist"))
ids_fp   <- file.path(ibs_dir, paste0("ibs_", RUN, ".mdist.id"))
meta_fp  <- file.path("01-info_files","samples_clean.txt")

out_dir  <- "04-results"
dir.create(out_dir, showWarnings = FALSE)

# Load data
ids <- read.table(ids_fp, header = FALSE)[[2]]
Dg  <- as.matrix(read.table(mdist_fp))
rownames(Dg) <- colnames(Dg) <- ids

meta <- read.delim(meta_fp, sep = "\t", check.names = FALSE)
meta <- meta %>% filter(IID %in% ids) %>% arrange(match(IID, ids))

coords <- meta %>% select(Longitude, Latitude) %>% as.matrix()
Dgeo_km <- distm(coords, fun = distHaversine) / 1000
rownames(Dgeo_km) <- colnames(Dgeo_km) <- ids

# Mantel tests
set.seed(123)
dg_gen <- as.dist(Dg)
dg_geo <- as.dist(Dgeo_km)
dg_ln  <- as.dist(log(Dgeo_km + 1e-6))

mantel_all  <- mantel(dg_gen, dg_geo,   permutations = 9999)
mantel_lnkm <- mantel(dg_gen, dg_ln,    permutations = 9999)

nonW <- meta$IID[meta$Group != "Western_Europe"]
mantel_nonW <- mantel(as.dist(Dg[nonW,nonW]), as.dist(log(Dgeo_km[nonW,nonW] + 1e-6)), permutations = 9999)

summary_tbl <- tibble(
  Test = c("Global (km)", "Global (ln km)", "Excl Western (ln km)"),
  r = c(mantel_all$statistic, mantel_lnkm$statistic, mantel_nonW$statistic),
  p = c(mantel_all$signif, mantel_lnkm$signif, mantel_nonW$signif)
)
write.table(summary_tbl, file.path(out_dir,"ibd_mantel_summary.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Within-group Mantel
wg <- unique(meta$Group)
wg_tbl <- map_dfr(wg, function(g) {
  ids_g <- meta$IID[meta$Group == g]
  if (length(ids_g) < 5) return(tibble(Group=g, n=length(ids_g), r=NA, p=NA))
  m <- mantel(as.dist(Dg[ids_g,ids_g]), as.dist(log(Dgeo_km[ids_g,ids_g] + 1e-6)), permutations=9999)
  tibble(Group=g, n=length(ids_g), r=m$statistic, p=m$signif)
})
write.table(wg_tbl, file.path(out_dir,"ibd_within_group_mantel.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Plots
pairs_df <- tibble(
  i = rownames(Dgeo_km)[row(Dgeo_km)[upper.tri(Dgeo_km)]],
  j = colnames(Dgeo_km)[col(Dgeo_km)[upper.tri(Dgeo_km)]],
  geo_km = Dgeo_km[upper.tri(Dgeo_km)],
  gen_1ibs = Dg[upper.tri(Dg)],
  group_i = meta$Group[match(rownames(Dgeo_km)[row(Dgeo_km)[upper.tri(Dgeo_km)]], meta$IID)],
  group_j = meta$Group[match(colnames(Dgeo_km)[col(Dgeo_km)[upper.tri(Dgeo_km)]], meta$IID)]
) %>%
  mutate(same_group = group_i == group_j)

# All pairs
p_all <- ggplot(pairs_df, aes(geo_km, gen_1ibs)) +
  geom_point(shape=21, color="black", size=2) +
  geom_smooth(method="lm", color="black", fill="grey80") +
  labs(x="Geographic distance (km)", y="Genetic distance (1–IBS)") +
  theme_classic()
ggsave(file.path(out_dir,"ibd_all_pairs.svg"), p_all, width=6.5, height=5)

# Excluding Western Europe
pairs_nonW <- pairs_df %>% filter(!(group_i == "Western_Europe" | group_j == "Western_Europe"))
p_nonW <- ggplot(pairs_nonW, aes(geo_km, gen_1ibs)) +
  geom_point(shape=21, color="black", size=2) +
  geom_smooth(method="lm", color="black", fill="grey80") +
  labs(x="Geographic distance (km)", y="Genetic distance (1–IBS)") +
  theme_classic()
ggsave(file.path(out_dir,"ibd_excluding_western.svg"), p_nonW, width=6.5, height=5)

# Within-group facets
p_facets <- pairs_df %>% filter(same_group) %>%
  ggplot(aes(log(geo_km + 1e-6), gen_1ibs)) +
  geom_point(shape=21, fill="white", color="black", size=1.8, alpha=0.7) +
  geom_smooth(method="lm", color="black", fill="grey85") +
  facet_wrap(~group_i, scales="free_x") +
  labs(x="ln(geographic distance, km)", y="Genetic distance (1–IBS)") +
  theme_classic()
ggsave(file.path(out_dir,"ibd_within_groups_facets.svg"), p_facets, width=8.5, height=6)

message("IBD analysis complete. Results written to 04-results/")