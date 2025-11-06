#!/usr/bin/env Rscript
# ==============================================================
# Script: 06_vcftools_report.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Quick visualization and summary of VCFtools QC outputs
#   produced by 07_core_vcf_stats.sh for final runs.
#
#   For each <run_name> given, reads files from:
#     03-analysis/06-finalrun/<run_name>/qc_vcftools/
#   and writes results to:
#     04-results/qc_vcftools/<run_name>/
#
# Outputs (per run):
#   - <run_name>_qc_plots.pdf    → histograms & densities
#   - <run_name>_qc_summary.tsv  → descriptive statistics
#
# Usage:
#   Rscript 00-scripts/08_qc_vcftools_report.R <run_name> [<run_name> ...]
# ==============================================================

library(tidyverse)

# Parse command-line arguments

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 00-scripts/08_qc_vcftools_report.R <run_name> [<run_name> ...]")
}

# input/output paths
input_root  <- "03-analysis/06-finalrun"
output_root <- "04-results/qc_vcftools"

# helper: compute quick descriptive stats for numeric vectors
summarize_vec <- function(x) {
  x <- x[is.finite(x)]
  tibble(
    stat  = c("n","mean","median","sd","min","p25","p75","max"),
    value = c(length(x), mean(x), median(x), sd(x), min(x),
              quantile(x, 0.25), quantile(x, 0.75), max(x))
  )
}

# Loop over each requested run
for (run in args) {

  # define input and output folders
  in_dir  <- file.path(input_root, run, "qc_vcftools")
  out_dir <- file.path(output_root, run)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Read VCFtools outputs for this run
  depth_site <- read_delim(file.path(in_dir, paste0(run, "_depth_site.ldepth.mean")),
                           delim="\t", skip=1, col_names=c("CHROM","POS","MEAN_DEPTH","VAR_DEPTH"))
  depth_indv <- read_delim(file.path(in_dir, paste0(run, "_depth_indv.idepth")),
                           delim="\t", skip=1, col_names=c("IND","NSITES","MEAN_DEPTH"))
  miss_site  <- read_delim(file.path(in_dir, paste0(run, "_missing_site.lmiss")),
                           delim="\t", skip=1, col_names=c("CHROM","POS","NCHR","N_FILTERED","NMISS","F_MISS"))
  miss_indv  <- read_delim(file.path(in_dir, paste0(run, "_missing_indv.imiss")),
                           delim="\t", skip=1, col_names=c("IND","NDATA","NFILTERED","NMISS","F_MISS"))
  freq2      <- read_delim(file.path(in_dir, paste0(run, "_freq2.frq")),
                           delim="\t", skip=1, col_names=c("CHROM","POS","N_ALLELES","N_CHR","A1","A2"))
  het        <- read_delim(file.path(in_dir, paste0(run, "_het.het")),
                           delim="\t", skip=1, col_names=c("IND","O_HET","E_HET","NSITES","F"))
  sitequal   <- tryCatch(
                    read_delim(file.path(in_dir, paste0(run, "_sitequal.lqual")),
                               delim="\t", skip=1, col_names=c("CHROM","POS","QUAL")),
                    error=function(e) NULL)

  # Generate summary plots (multi-page PDF)
  pdf(file.path(out_dir, paste0(run, "_qc_plots.pdf")), width=8.5, height=6)
  theme_set(theme_minimal())

  # mean depth per site
  ggplot(depth_site, aes(MEAN_DEPTH)) +
    geom_histogram(bins=60) +
    labs(title=paste(run,"- Mean Depth per Site"),
         x="Mean Depth per Locus", y="Count") %>%
    print()

  # mean depth per individual
  ggplot(depth_indv, aes(MEAN_DEPTH)) +
    geom_histogram(bins=40) +
    labs(title=paste(run,"- Mean Depth per Individual"),
         x="Mean Depth", y="Count") %>%
    print()

  # missingness per site
  ggplot(miss_site, aes(F_MISS)) +
    geom_density(fill="orange", alpha=0.4) +
    labs(title=paste(run,"- Missingness per Site"),
         x="Proportion Missing per Site", y="Density") %>%
    print()

  # missingness per individual
  ggplot(miss_indv, aes(F_MISS)) +
    geom_histogram(bins=40, fill="tomato", color="black") +
    labs(title=paste(run,"- Missingness per Individual"),
         x="Proportion Missing", y="Count") %>%
    print()

  # minor allele frequency (MAF)
  maf <- freq2 %>% mutate(MAF = pmin(A1, A2))
  ggplot(maf, aes(MAF)) +
    geom_density(fill="purple", alpha=0.4) +
    labs(title=paste(run,"- Minor Allele Frequency"),
         x="MAF", y="Density") %>%
    print()

  # inbreeding coefficient (F)
  ggplot(het, aes(F)) +
    geom_histogram(bins=40, fill="darkolivegreen", color="black") +
    labs(title=paste(run,"- Inbreeding Coefficient (F)"),
         x="F", y="Count") %>%
    print()

  # F per individual
  ggplot(het, aes(x=reorder(IND,F), y=F)) +
    geom_point(color="darkolivegreen", size=2) +
    geom_segment(aes(xend=IND, y=0, yend=F),
                 color="darkolivegreen", alpha=0.5) +
    labs(title=paste(run,"- Inbreeding per Individual"),
         x="Individual", y="F") +
    coord_flip() %>%
    print()


  # Write numeric summaries for quick inspection
  summary_tbl <- bind_rows(
    summarize_vec(depth_site$MEAN_DEPTH) %>% mutate(section="depth_site", variable="MEAN_DEPTH"),
    summarize_vec(depth_indv$MEAN_DEPTH) %>% mutate(section="depth_indv", variable="MEAN_DEPTH"),
    summarize_vec(miss_site$F_MISS) %>% mutate(section="miss_site", variable="F_MISS"),
    summarize_vec(miss_indv$F_MISS) %>% mutate(section="miss_indv", variable="F_MISS"),
    summarize_vec(maf$MAF) %>% mutate(section="maf", variable="MAF"),
    summarize_vec(het$F) %>% mutate(section="het", variable="F")
  )

  if (!is.null(sitequal))
    summary_tbl <- bind_rows(summary_tbl,
                             summarize_vec(sitequal$QUAL) %>% mutate(section="sitequal", variable="QUAL"))

  write_tsv(summary_tbl, file.path(out_dir, paste0(run, "_qc_summary.tsv")))
}

