#!/usr/bin/env Rscript
# ==============================================================
# Script: 13_ibd_analysis.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
# Isolation-by-distance (IBD) using PLINK 1–IBS and great-circle distances,
# Mantel tests (9,999 perms), diagnostic plots and AMOVA (among broad regions
# and with Sweden split).
#
# Inputs:
#   03-analysis/08-downstream/05-ibd/ibs_<RUN>.mdist(.id)
#   01-info_files/samples_clean.txt
#
# Outputs:
#   04-results/ibd_mantel_summary.tsv
#   04-results/ibd_within_group_mantel.tsv
#   04-results/ibd_all_pairs.svg
#   04-results/ibd_excluding_western.svg
#   04-results/ibd_within_groups_facets.svg
#   04-results/amova_group_summary.tsv
#   04-results/amova_sweden_group2_summary.tsv
#
# Usage:
#   Rscript 00-scripts/13_ibd_analysis.R r70_noTrentino_noOutgroup
# ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(geosphere)  
  library(vegan)      
  library(pegas)      
  library(svglite)
})

# USER Argument
RUN <- commandArgs(trailingOnly = TRUE)[1]
if (is.null(RUN) || RUN == "") stop("Provide RUN, e.g. r70_noTrentino_noOutgroup", call. = FALSE)

# Paths
ibs_dir  <- file.path("03-analysis","08-downstream","05-ibd")
mdist_fp <- file.path(ibs_dir, paste0("ibs_", RUN, ".mdist"))
ids_fp   <- file.path(ibs_dir, paste0("ibs_", RUN, ".mdist.id"))
meta_fp  <- file.path("01-info_files","samples_clean.txt")
out_dir  <- "04-results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read & align
ids <- read.table(ids_fp, header = FALSE, stringsAsFactors = FALSE)[[2]]
Dg  <- as.matrix(read.table(mdist_fp, header = FALSE, check.names = FALSE))
stopifnot(nrow(Dg) == length(ids), ncol(Dg) == length(ids))
rownames(Dg) <- colnames(Dg) <- ids

meta <- read.delim(meta_fp, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
need_cols <- c("IID","Latitude","Longitude","Group","Country","Group2")
stopifnot(all(need_cols %in% names(meta)))
meta <- meta %>% filter(IID %in% ids) %>% arrange(match(IID, ids))
stopifnot(identical(meta$IID, ids))

# Geographic distances (km)
coords     <- meta %>% select(Longitude, Latitude) %>% as.matrix()
Dgeo_km    <- distm(coords, fun = distHaversine) / 1000
rownames(Dgeo_km) <- colnames(Dgeo_km) <- ids

# Mantel tests
set.seed(123)
eps      <- 1e-6
dg_gen   <- as.dist(Dg)
dg_geo   <- as.dist(Dgeo_km)
dg_lngeo <- as.dist(log(Dgeo_km + eps))

m_all_km <- vegan::mantel(dg_gen, dg_geo,   permutations = 9999, method = "pearson")
m_all_ln <- vegan::mantel(dg_gen, dg_lngeo, permutations = 9999, method = "pearson")

# Exclude Western_Europe Group
ids_nonW <- meta$IID[meta$Group != "Western_Europe"]
m_nonW <- if (length(ids_nonW) >= 4) {
  vegan::mantel(as.dist(Dg[ids_nonW, ids_nonW]),
                as.dist(log(Dgeo_km[ids_nonW, ids_nonW] + eps)),
                permutations = 9999, method = "pearson")
} else NULL

mantel_summary <- tibble(
  Test = c("Global (km)", "Global (ln km)", if (!is.null(m_nonW)) "Excl Western (ln km)"),
  r    = c(unname(as.numeric(m_all_km$statistic)),
           unname(as.numeric(m_all_ln$statistic)),
           if (!is.null(m_nonW)) unname(as.numeric(m_nonW$statistic)) else NULL),
  p    = c(m_all_km$signif, m_all_ln$signif, if (!is.null(m_nonW)) m_nonW$signif else NULL)
)
write.table(mantel_summary, file.path(out_dir, "ibd_mantel_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Within-group Mantel on ln(km)
wg <- map_dfr(sort(unique(meta$Group)), function(g) {
  ids_g <- meta$IID[meta$Group == g]
  if (length(ids_g) < 5) return(tibble(Group = g, n = length(ids_g), r = NA_real_, p = NA_real_))
  m <- vegan::mantel(as.dist(Dg[ids_g, ids_g]),
                     as.dist(log(Dgeo_km[ids_g, ids_g] + eps)),
                     permutations = 9999, method = "pearson")
  tibble(Group = g, n = length(ids_g), r = unname(as.numeric(m$statistic)), p = m$signif)
})
write.table(wg, file.path(out_dir, "ibd_within_group_mantel.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Plotting data
UT <- upper.tri(Dgeo_km)
pairs_df <- tibble(
  i        = rownames(Dgeo_km)[row(Dgeo_km)[UT]],
  j        = colnames(Dgeo_km)[col(Dgeo_km)[UT]],
  geo_km   = Dgeo_km[UT],
  ln_geo   = log(Dgeo_km[UT] + eps),
  gen_1ibs = Dg[UT]
) %>%
  mutate(
    group_i    = meta$Group[match(i, meta$IID)],
    group_j    = meta$Group[match(j, meta$IID)],
    same_group = group_i == group_j
  )

# All pairs (km)
p_all <- ggplot(pairs_df, aes(geo_km, gen_1ibs)) +
  geom_point(shape = 21, size = 2, stroke = 0.6, fill = NA, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
  labs(x = "Geographic distance (km)", y = "Genetic distance (1–IBS)") +
  theme_classic()
ggsave(file.path(out_dir, "ibd_all_pairs.svg"), p_all, width = 6.6, height = 5.2, dpi = 300)

# Excluding Western_Europe
pairs_nonW <- pairs_df %>% filter(!(group_i == "Western_Europe" | group_j == "Western_Europe"))
if (nrow(pairs_nonW) > 0) {
  p_nonW <- ggplot(pairs_nonW, aes(geo_km, gen_1ibs)) +
    geom_point(shape = 21, size = 2, stroke = 0.6, fill = NA, color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80") +
    labs(x = "Geographic distance (km)", y = "Genetic distance (1–IBS)") +
    theme_classic()
  ggsave(file.path(out_dir, "ibd_excluding_western.svg"), p_nonW, width = 6.6, height = 5.2, dpi = 300)
}

# Within-group facets (ln km)
p_facets <- pairs_df %>%
  filter(same_group) %>%
  ggplot(aes(ln_geo, gen_1ibs)) +
  geom_point(alpha = .7, shape = 21, size = 1.8, fill = "white", color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey85") +
  facet_wrap(~ group_i, scales = "free_x") +
  labs(x = "ln(geographic distance, km)", y = "Genetic distance (1–IBS)") +
  theme_classic()
ggsave(file.path(out_dir, "ibd_within_groups_facets.svg"), p_facets, width = 9, height = 6.5, dpi = 300)

# AMOVA among groups
D_dist <- as.dist(Dg)
grp    <- factor(setNames(meta$Group, meta$IID)[labels(D_dist)])

set.seed(123)
am1 <- pegas::amova(D_dist ~ grp, nperm = 9999)

tab    <- am1$tab
sigma2 <- tab[, "sigma2"]
perc   <- 100 * sigma2 / sum(sigma2)
phi_st <- as.numeric(sigma2["grp"] / sum(sigma2))

am_summ <- tibble(
  component = rownames(tab),
  df        = unname(tab[, "Df"]),
  SS        = unname(tab[, "SS"]),
  sigma2    = unname(sigma2),
  perc_var  = unname(perc)
) %>%
  bind_rows(tibble(component = "Phi_ST", df = NA, SS = NA, sigma2 = NA, perc_var = phi_st))

write.table(am_summ,
            file.path(out_dir, "amova_group_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# AMOVA within Sweden Split
sw_keep <- meta$IID[meta$Country == "Sweden" & !is.na(meta$Group2)]
if (length(sw_keep) < 6) {
  stop("Required AMOVA (Sweden ~ Group2) cannot run: need >= 6 Swedish samples with non-NA Group2.", call. = FALSE)
}

D_se <- as.dist(Dg[sw_keep, sw_keep])
g2   <- factor(setNames(meta$Group2, meta$IID)[labels(D_se)])

set.seed(123)
am2 <- pegas::amova(D_se ~ g2, nperm = 9999)

tab2     <- am2$tab
sigma2_2 <- tab2[, "sigma2"]
perc2    <- 100 * sigma2_2 / sum(sigma2_2)
phi_st2  <- as.numeric(sigma2_2["g2"] / sum(sigma2_2))

am2_summ <- tibble(
  component = rownames(tab2),
  df        = unname(tab2[, "Df"]),
  SS        = unname(tab2[, "SS"]),
  sigma2    = unname(sigma2_2),
  perc_var  = unname(perc2)
) %>%
  bind_rows(tibble(component = "Phi_ST", df = NA, SS = NA, sigma2 = NA, perc_var = phi_st2))

write.table(am2_summ,
            file.path(out_dir, "amova_sweden_group2_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("IBD + AMOVA complete → outputs in 04-results/")
