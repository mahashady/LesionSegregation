#!/usr/bin/env Rscript

source("utils/obs_vs_exp_utils.R")

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(lsa)
  library(pheatmap)
  library(grid)
})

# ----------------------------
# CLI
# ----------------------------
option_list <- list(
  make_option(
    "--treatment",
    type = "character",
    default = NULL,
    help = "Treatment: Alkylating or Platinum"
  ),
  make_option(
    "--results-dir",
    type = "character",
    default = "../results/multi-allelic_sites/spectra",
    help = "Base results directory [default %default]"
  ),
  make_option(
    "--plots-dir",
    type = "character",
    default = "../plots",
    help = "Base plots directory [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$treatment)) {
  stop("--treatment is required")
}
if (!opt$treatment %in% c("Alkylating", "Platinum")) {
  stop("--treatment must be one of: Alkylating, Platinum")
}

# ----------------------------
# Main analysis wrapper
# ----------------------------
run_one_analysis <- function(cfg) {
  print(cfg$expected_labels["enriched"])

  df_expected <- safe_full_join(list(
    read_expected_file(
    cfg$expected_files$enriched,
    unname(cfg$expected_labels["enriched"])
    ),
    read_expected_file(
      cfg$expected_files$nonenriched,
      unname(cfg$expected_labels["nonenriched"])
    )
  ), "context")

  print(head(df_expected))
  df_observed <- safe_full_join(list(
    read_observed_file(
      cfg$observed_files$enriched,
      unname(cfg$observed_labels["enriched"]),
      keep_type = TRUE
    ),
    read_observed_file(
      cfg$observed_files$nonenriched,
      unname(cfg$observed_labels["nonenriched"]),
      keep_type = FALSE
    )
  ), "mutation")

  observed_merge_cols <- c("mutation", unname(cfg$observed_labels))
  if ("Type" %in% colnames(df_observed)) {
    observed_merge_cols <- c("mutation", "Type", unname(cfg$observed_labels))
  }
  print(head(df_observed))
  result <- merge(
    df_expected[, c("context", unname(cfg$expected_labels))],
    df_observed[, observed_merge_cols],
    by.x = "context",
    by.y = "mutation"
  )

  df_signatures <- read_signature_table(
    cfg$signature_file,
    cfg$signature_cols
  )
  print(head(df_signatures))
  result <- merge(result, df_signatures, by = "context")
  print(head(result))

  for (old_name in names(cfg$signature_labels)) {
    colnames(result)[colnames(result) == old_name] <- cfg$signature_labels[[old_name]]
  }

  observed_cols <- unname(cfg$observed_labels)
  comparison_cols <- c(
    unname(cfg$expected_labels),
    unname(cfg$signature_labels)
  )

  mat <- compute_cosine_matrix(result, observed_cols, comparison_cols)

  plot_heatmap(mat, cfg$plot_file)

  print(result)
  print(mat)
  cat("Saved plot to:", cfg$plot_file, "\n")
}

# ----------------------------
# Run
# ----------------------------
configs <- get_config(
  mode = opt$treatment,
  results_dir = opt$`results-dir`,
  plots_dir = opt$`plots-dir`
)

for (cfg in configs) {
  cat("Running:", cfg$name, "\n")
  run_one_analysis(cfg)
}