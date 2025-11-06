#!/usr/bin/env Rscript
# ==============================================================
# Script: 05_parameter_choice.R
# Author: Tabea Dittmar
# Year:   2025
#
# Description:
#   Create combined RADstackshelpR panels for Stacks parameter tuning
#   using outputs produced by the optimization bash scripts:
#
#   m runs: 03-analysis/05-parameters/m/m<m>/populations.snps.vcf
#   M runs: 03-analysis/05-parameters/M/M<M>/populations.snps.vcf
#   n runs: 03-analysis/05-parameters/n/n<n>/populations.snps.vcf
#
#   - mode = "m": plots Depth • SNPs • Loci across all tested m values
#   - mode = "M": plots SNPs • Loci across all tested M values
#   - mode = "n": plots SNPs • Loci for the two tested n values (n=M, n=M+1)
#
# Outputs (PNG only, white bg) written directly to 04-results/:
#   04-results/m_panel.png
#   04-results/gM_panel.png
#   04-results/n_panel.png
#
# Usage:
#   # scan a directory of runs (recommended)
#   Rscript 00-scripts/05_parameter_choice.R m 03-analysis/05-parameters/m
#   Rscript 00-scripts/05_parameter_choice.R M 03-analysis/05-parameters/M
#   Rscript 00-scripts/05_parameter_choice.R n 03-analysis/05-parameters/n
#
#   # OR (n only) give two explicit VCF files:
#   Rscript 00-scripts/05_parameter_choice.R n \
#     03-analysis/05-parameters/n/n4/populations.snps.vcf \
#     03-analysis/05-parameters/n/n5/populations.snps.vcf
# ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(RADstackshelpR)
  library(patchwork)
})

# Parse args & utilities
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("See header for usage.", call. = FALSE)

mode <- args[[1]]   # "m" | "M" | "n"
arg2 <- args[[2]]

out_dir <- "04-results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# small helper to save PNG consistently
save_png <- function(plot, path, width_in, height_in, dpi = 300) {
  ggplot2::ggsave(
    filename = path, plot = plot, device = "png",
    width = width_in, height = height_in, dpi = dpi, bg = "white"
  )
}

# helper: find VCF inside subdirs that match a pattern
find_param_vcfs <- function(root, subdir_regex) {
  # list subdirectories like m3, m4, ... or M2, M3, ... or n4, n5, ...
  subs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  subs <- subs[grepl(subdir_regex, basename(subs))]
  tibble(
    sub  = subs,
    vcf  = file.path(subs, "populations.snps.vcf"),
    tag  = basename(subs)
  ) |>
    filter(file.exists(vcf))
}

# ------------------------------------------------------
# 1) MODE = "m" — expect subdirs: m3, m4, ... each with populations.snps.vcf
# ------------------------------------------------------
if (mode == "m") {
  root <- arg2
  runs <- find_param_vcfs(root, "^m[0-9]+$")
  if (nrow(runs) < 2) stop("Need ≥2 m runs with VCFs under: ", root)

  # parse numeric m value from subdir name
  runs <- runs |>
    mutate(m = as.integer(sub("^m(\\d+)$", "\\1", tag))) |>
    arrange(m)

  # build named argument list for optimize_m: c(m3="path", m4="path", ...)
  args_list <- as.list(runs$vcf)
  names(args_list) <- paste0("m", runs$m)

  # RADstackshelpR: standardizes inputs & builds ggplots with best-fit markers
  res     <- do.call(RADstackshelpR::optimize_m, args_list)
  p_depth <- RADstackshelpR::vis_depth(output = res) + labs(y = "Average depth per sample")
  p_snps  <- RADstackshelpR::vis_snps(output = res, stacks_param = "m")
  p_loci  <- RADstackshelpR::vis_loci(output = res, stacks_param = "m")

  panel <- p_depth | p_snps | p_loci
  out_png <- file.path(out_dir, "m_panel.png")
  save_png(panel, out_png, width_in = 18, height_in = 5)
  cat("Saved:", out_png, "\n")

# ------------------------------------------------------
# 2) MODE = "M" — expect subdirs: M2, M3, ... each with populations.snps.vcf
# ------------------------------------------------------
} else if (mode == "M") {
  root <- arg2
  runs <- find_param_vcfs(root, "^M[0-9]+$")
  if (nrow(runs) < 2) stop("Need ≥2 M runs with VCFs under: ", root)

  runs <- runs |>
    mutate(M = as.integer(sub("^M(\\d+)$", "\\1", tag))) |>
    arrange(M)

  args_list <- as.list(runs$vcf)
  names(args_list) <- paste0("M", runs$M)

  res    <- do.call(RADstackshelpR::optimize_bigM, args_list)
  p_snps <- RADstackshelpR::vis_snps(output = res, stacks_param = "M")
  p_loci <- RADstackshelpR::vis_loci(output = res, stacks_param = "M")

  panel <- p_snps | p_loci
  out_png <- file.path(out_dir, "gM_panel.png")
  save_png(panel, out_png, width_in = 12, height_in = 5)
  cat("Saved:", out_png, "\n")

# ------------------------------------------------------
# 3) MODE = "n" — either:
#    A) arg2 is a root dir with subdirs n<M> and n<M+1>, or
#    B) two explicit VCF paths are provided
# ------------------------------------------------------
} else if (mode == "n") {
  # case A: directory with n-runs
  if (length(args) == 2 && dir.exists(arg2)) {
    root <- arg2
    runs <- find_param_vcfs(root, "^n[0-9]+$")
    if (nrow(runs) != 2) stop("Expected exactly two n runs with VCFs under: ", root)

    runs <- runs |>
      mutate(n = as.integer(sub("^n(\\d+)$", "\\1", tag))) |>
      arrange(n)

    vcf_M  <- runs$vcf[1]  # lower n (assumed n=M)
    vcf_M1 <- runs$vcf[2]  # higher n (assumed n=M+1)

  # case B: two files given explicitly
  } else if (length(args) >= 4) {
    vcf_M  <- arg2
    vcf_M1 <- args[[3]]
    if (!file.exists(vcf_M)  || !grepl("\\.vcf(\\.gz)?$", vcf_M))  stop("First VCF not found or not .vcf(.gz): ", vcf_M)
    if (!file.exists(vcf_M1) || !grepl("\\.vcf(\\.gz)?$", vcf_M1)) stop("Second VCF not found or not .vcf(.gz): ", vcf_M1)
  } else {
    stop('For mode "n", provide either a directory with two runs or two VCF paths.', call. = FALSE)
  }

  res    <- RADstackshelpR::optimize_n(nequalsM = vcf_M, nequalsMplus1 = vcf_M1)
  p_snps <- RADstackshelpR::vis_snps(output = res, stacks_param = "n")
  p_loci <- RADstackshelpR::vis_loci(output = res, stacks_param = "n")

  panel <- p_snps | p_loci
  out_png <- file.path(out_dir, "n_panel.png")
  save_png(panel, out_png, width_in = 12, height_in = 5)
  cat("Saved:", out_png, "\n")

} else {
  stop('mode must be one of: "m", "M", "n"', call. = FALSE)
}

