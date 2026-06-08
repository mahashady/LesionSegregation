library(ggplot2)
library(dplyr)
library(ggridges)
library(patchwork)

# =========================
# Configuration
# =========================
recomb_lambda <- "1"
model <- "no_recombination"  # "no_recombination" or "recombination"

lookup_base_prefix <- "../../results/LAD_inference/likelihood_matrices"

samples_file <- "../../results/LAD_inference/data/multi_stat_by_chr.enriched.platinum.SBS17b_prop_lt_0.1.tsv"

out_dir <- "../../results/LAD_inference/inference_results/examples_with_gt_4_multisites/"

min_M <- 5
n_bins <- 15

# =========================
# Paths
# =========================
make_lookup_path <- function(model, recomb_lambda, lookup_base_prefix) {
  if (model == "no_recombination") {
    file.path(
      lookup_base_prefix,
      "likelihood_matrix_for_LAD_inference.no_recombination.M2-15.LAD1-5.n1000.tsv"
    )
  } else if (model == "recombination") {
    file.path(
      lookup_base_prefix,
      paste0(
        "likelihood_matrix_for_LAD_inference.recombination.M2-15.LAD1-5.lambda",
        recomb_lambda,
        ".0.win10000000.n1000.tsv"
      )
    )
  } else {
    stop("Unknown model: ", model)
  }
}

read_tsv_checked <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }

  read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# =========================
# Plot one sample
# =========================
plot_one_example <- function(sample, df_lookup_sample, observed_C, n_bins = 15) {
  if (nrow(df_lookup_sample) == 0) {
    message("Skipping ", sample, ": no lookup rows for this M.")
    return(NULL)
  }

  ggplot(df_lookup_sample, aes(x = C, y = factor(LAD, levels = c(5, 4, 3, 2, 1)), height = prob)) +
    geom_ridgeline(
      fill = "grey80",
      color = "grey40",
      alpha = 0.7,
      scale = 1.5
    )  +
    geom_vline(
      xintercept = observed_C,
      linewidth = 0.8,
      color = "black"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    labs(
      title = sample,
      x = "N chromosomes with multiallelic sites",
      y = "LAD"
    ) 
  }

# =========================
# Create all sample plots
# =========================
make_sample_plots <- function(df_samples, df_lookup, n_bins = 15) {
  plot_list <- list()

  for (i in seq_len(nrow(df_samples))) {
    sample_id <- df_samples[i, "sample"]
    n_multi <- df_samples[i, "M"]
    observed_C <- df_samples[i, "C"]

    df_lookup_sample <- df_lookup %>%
      filter(M == n_multi)

    p <- plot_one_example(
      sample = sample_id,
      df_lookup_sample = df_lookup_sample,
      observed_C = observed_C,
      n_bins = n_bins
    )

    if (!is.null(p)) {
      plot_list[[sample_id]] <- p
    }
  }

  plot_list
}

# =========================
# Save combined plots
# =========================
save_combined_plots <- function(plot_list, out_dir, model, recomb_lambda) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  if (length(plot_list) == 0) {
    stop("No plots to save.")
  }

  title_text <- ifelse(
    model == "no_recombination",
    "No recombination",
    paste0(recomb_lambda, " recombination")
  )
  combined <- wrap_plots(plot_list, ncol = 3) +
    plot_annotation(
      title = title_text,
      #subtitle = "Observed vs simulated chromosome counts",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )

  model_label <- ifelse(
    model == "no_recombination",
    "no_recombination",
    paste0("recombination_lambda", recomb_lambda)
  )

  out_pdf <- file.path(
    out_dir,
    paste0("LAD_lookup_examples.", model_label, ".pdf")
  )

  out_png <- file.path(
    out_dir,
    paste0("LAD_lookup_examples.", model_label, ".png")
  )

  ggsave(
    out_pdf,
    combined,
    width = 12,
    height = max(6, ceiling(length(plot_list) / 3) * 4),
    units = "in"
  )

  ggsave(
    out_png,
    combined,
    width = 10,
    height = max(4, ceiling(length(plot_list) / 3) * 2),
    units = "in",
    dpi = 300
  )

  message("Saved combined PDF to: ", out_pdf)
  message("Saved combined PNG to: ", out_png)

  combined
}

# =========================
# Main
# =========================
main <- function() {
  lookup_table <- make_lookup_path(
    model = model,
    recomb_lambda = recomb_lambda,
    lookup_base_prefix = lookup_base_prefix
  )

  df_lookup <- read_tsv_checked(lookup_table)
  df_samples <- read_tsv_checked(samples_file) %>%
    filter(M >= min_M)

  message("Loaded lookup rows: ", nrow(df_lookup))
  message("Samples to plot: ", nrow(df_samples))

  required_lookup_cols <- c("M", "C", "LAD")
  required_sample_cols <- c("sample", "M", "C")

  if (!all(required_lookup_cols %in% colnames(df_lookup))) {
    stop("Lookup table missing columns: ",
         paste(setdiff(required_lookup_cols, colnames(df_lookup)), collapse = ", "))
  }

  if (!all(required_sample_cols %in% colnames(df_samples))) {
    stop("Samples table missing columns: ",
         paste(setdiff(required_sample_cols, colnames(df_samples)), collapse = ", "))
  }

  plot_list <- make_sample_plots(
    df_samples = df_samples,
    df_lookup = df_lookup,
    n_bins = n_bins
  )

  save_combined_plots(
    plot_list = plot_list,
    out_dir = out_dir,
    model = model,
    recomb_lambda = recomb_lambda
  )
}

main()
