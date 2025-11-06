#!/usr/bin/env Rscript
# ==============================================================
# Script: 13_ibd_analysis.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Isolation-by-distance analysis: read PLINK 1–IBS distances and
#   sample coordinates, compute geographic distances, run Mantel tests
#   (9999 perms), and plot relationships (global, excl. Western_Europe,
#   and within-group facets). Saves plots and result tables to 04-results/.
#
# Inputs:
#   - 03-analysis/08-downstream/04-diversity/ibs/ibs_<RUN>.mdist
#   - 03-analysis/08-downstream/04-diversity/ibs/ibs_<RUN>.mdist.id
#   - 01-info_files/samples_clean.txt  (must contain: IID, Latitude, Longitude, Group[, Country])
#
# Outputs:
#   - 04-results/ibd_all_pairs.svg
#   - 04-results/ibd_excluding_western.svg
#   - 04-results/ibd_within_groups_facets.svg
#   - 04-results/ibd_mantel_summary.tsv
#   - 04-results/ibd_within_group_mantel.tsv
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
  stop("Usage: Rscript 00-scripts/18_ibd_analysis.R <RUN>\n",
       "Example: Rscript 00-scripts/18_ibd_analysis.R r70_noTrentino_noOutgroup", call. = FALSE)
}
RUN <- args[[1]]

# Paths
base_ibd <- file.path("03-analysis","08-downstream","04-diversity","ibs")
mdist_fp <- file.path(base_ibd, paste0("ibs_", RUN, ".mdist"))
ids_fp   <- file.path(base_ibd, paste0("ibs_", RUN, ".mdist.id"))
meta_fp  <- file.path("01-info_files","samples_clean.txt")

rout_dir <- "04-results"
dir.create(rout_dir, showWarnings = FALSE)

# Read genetic distance (1–IBS)
if (!file.exists(mdist_fp) || !file.exists(ids_fp)) {
  stop("Missing PLINK outputs. Run 17_compute_ibs.sh first.\n",
       "Expected:\n  ", mdist_fp, "\n  ", ids_fp)
}
ids_df <- read.table(ids_fp, header = FALSE, stringsAsFactors = FALSE)
ids <- ids_df[[2]]  # take IID column
Dg  <- as.matrix(read.table(mdist_fp, header = FALSE))
stopifnot(nrow(Dg) == length(ids), ncol(Dg) == length(ids))
rownames(Dg) <- colnames(Dg) <- ids

# Read metadata & geographic distances (km)
meta <- read.delim(meta_fp, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
need_cols <- c("IID","Latitude","Longitude","Group")
if (!all(need_cols %in% names(meta))) {
  stop("samples_clean.txt must contain columns: ", paste(need_cols, collapse=", "))
}

meta_ord <- meta %>% filter(IID %in% ids) %>% arrange(match(IID, ids))
stopifnot(identical(meta_ord$IID, ids))

coords <- meta_ord %>% select(Longitude, Latitude) %>% as.matrix()
Dgeo_km <- distm(coords, fun = distHaversine) / 1000
rownames(Dgeo_km) <- colnames(Dgeo_km) <- ids

# Helpers
eps <- 1e-6
upper_tri <- function(M) M[upper.tri(M)]

dg_gen   <- as.dist(Dg)
dg_geo   <- as.dist(Dgeo_km)
dg_geoLN <- as.dist(log(Dgeo_km + eps))

# Mantel tests (global: raw & ln; excl. Western; within groups)
set.seed(123)
mantel_raw <- mantel(dg_gen, dg_geo,   permutations = 9999, method = "pearson")
mantel_ln  <- mantel(dg_gen, dg_geoLN, permutations = 9999, method = "pearson")

summ_rows <- list(
  tibble(test = "Global (km)",         r = as.numeric(mantel_raw$statistic), p = mantel_raw$signif,
         n_pairs = length(dg_gen)),
  tibble(test = "Global (ln km)",      r = as.numeric(mantel_ln$statistic),  p = mantel_ln$signif,
         n_pairs = length(dg_gen))
)

# Excluding Western_Europe
keep_ids_nonW <- meta_ord %>% filter(Group != "Western_Europe") %>% pull(IID)
if (length(keep_ids_nonW) >= 4) {
  Dg_nonW   <- as.dist(Dg[keep_ids_nonW, keep_ids_nonW])
  Dgeo_nonW <- as.dist(log(Dgeo_km[keep_ids_nonW, keep_ids_nonW] + eps))
  set.seed(123)
  m_nonW <- mantel(Dg_nonW, Dgeo_nonW, permutations = 9999, method = "pearson")
  summ_rows <- append(summ_rows, list(
    tibble(test = "Excl Western (ln km)", r = as.numeric(m_nonW$statistic),
           p = m_nonW$signif, n_pairs = length(Dg_nonW))
  ))
}

mantel_summary <- bind_rows(summ_rows)
write.table(mantel_summary,
            file = file.path(rout_dir, "ibd_mantel_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Within-group Mantels
groups <- sort(unique(meta_ord$Group))
wg_list <- purrr::map_dfr(groups, function(g) {
  ids_g <- meta_ord %>% filter(Group == g) %>% pull(IID)
  if (length(ids_g) < 5) return(tibble(Group = g, n = length(ids_g), r = NA_real_, p = NA_real_))
  Dg_g   <- as.dist(Dg[ids_g, ids_g])
  Dgeo_g <- as.dist(log(Dgeo_km[ids_g, ids_g] + eps))
  set.seed(123)
  res <- mantel(Dg_g, Dgeo_g, permutations = 9999, method = "pearson")
  tibble(Group = g, n = length(ids_g), r = as.numeric(res$statistic), p = res$signif)
})
write.table(wg_list,
            file = file.path(rout_dir, "ibd_within_group_mantel.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Pairwise table for plotting
pairs_df <- tibble::tibble(
  i = rownames(Dgeo_km)[row(Dgeo_km)[upper.tri(Dgeo_km)]],
  j = colnames(Dgeo_km)[col(Dgeo_km)[upper.tri(Dgeo_km)]],
  geo_km   = Dgeo_km[upper.tri(Dgeo_km)],
  ln_geo   = log(Dgeo_km[upper.tri(Dgeo_km)] + eps),
  gen_1ibs = Dg[upper.tri(Dg)]
) %>%
  mutate(
    group_i = meta_ord$Group[match(i, meta_ord$IID)],
    group_j = meta_ord$Group[match(j, meta_ord$IID)],
    same_group = group_i == group_j
  )

# PLOTS
# (P1) All pairs (raw km)
p_all <- ggplot(pairs_df, aes(x = geo_km, y = gen_1ibs)) +
  geom_point(shape = 21, color = "black", fill = NA, size = 2, stroke = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  labs(x = "Geographic distance (km)", y = "Genetic distance (1 – IBS)") +
  theme_classic()
ggsave(file.path(rout_dir, "ibd_all_pairs.svg"), p_all, width = 6.6, height = 5.2, dpi = 300)

# (P2) Excluding Western_Europe
ids_nonW <- meta_ord %>% filter(Group != "Western_Europe") %>% pull(IID)
Dg_nonW_m     <- Dg[ids_nonW, ids_nonW, drop = FALSE]
Dgeo_nonW_m   <- Dgeo_km[ids_nonW, ids_nonW, drop = FALSE]
pairs_nonW <- tibble::tibble(
  i = rownames(Dgeo_nonW_m)[row(Dgeo_nonW_m)[upper.tri(Dgeo_nonW_m)]],
  j = colnames(Dgeo_nonW_m)[col(Dgeo_nonW_m)[upper.tri(Dgeo_nonW_m)]],
  geo_km   = Dgeo_nonW_m[upper.tri(Dgeo_nonW_m)],
  gen_1ibs = Dg_nonW_m[upper.tri(Dg_nonW_m)]
)

p_nonW <- ggplot(pairs_nonW, aes(geo_km, gen_1ibs)) +
  geom_point(shape = 21, color = "black", fill = NA, size = 2, stroke = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  labs(x = "Geographic distance (km)", y = "Genetic distance (1 – IBS)") +
  theme_classic()
ggsave(file.path(rout_dir, "ibd_excluding_western.svg"), p_nonW, width = 6.6, height = 5.2, dpi = 300)

# (P3) Within-group facets (ln km)
p_facets <- pairs_df %>%
  filter(same_group) %>%
  ggplot(aes(ln_geo, gen_1ibs)) +
  geom_point(alpha = .65, shape = 21, fill = "white", color = "black", size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey85") +
  facet_wrap(~ group_i, scales = "free_x") +
  labs(x = "ln(geographic distance, km)", y = "Genetic distance (1 – IBS)") +
  theme_classic()
ggsave(file.path(rout_dir, "ibd_within_groups_facets.svg"), p_facets, width = 9, height = 6.5, dpi = 300)

message("Saved plots and tables in 04-results/")
