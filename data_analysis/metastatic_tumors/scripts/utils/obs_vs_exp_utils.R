library(dplyr)
normalize_expected <- function(df, new_col) {
  df[[new_col]] <- df$expected_count / sum(df$expected_count)
  df
}

normalize_observed <- function(df, new_col) {
  df[[new_col]] <- df$N_muts / sum(df$N_muts)
  df
}

read_expected_file <- function(path, new_col) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  normalize_expected(df, new_col) %>%
    select(context, all_of(new_col))
}

read_observed_file <- function(path, new_col, keep_type = FALSE) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- normalize_observed(df, new_col)

  if (keep_type && "Type" %in% colnames(df)) {
    df %>% select(mutation, Type, all_of(new_col))
  } else {
    df %>% select(mutation, all_of(new_col))
  }
}

safe_full_join <- function(df_list, by_col) {
  reduce(df_list, full_join, by = by_col)
}

get_config <- function(mode, results_dir, plots_dir) {

  if (mode == "Alkylating") {
    return(list(
      list(
        name = "Alkylating",
        expected_files = list(
          enriched = file.path(results_dir, "expected",
                               "enriched.Alkylating.expected_triallelic.csv"),
          nonenriched = file.path(results_dir, "expected",
                                  "NONenriched.Alkylating.expected_triallelic.csv")
        ),
        observed_files = list(
          enriched = file.path(results_dir, "multi_spectrum",
                               "multi_spectrum_enriched.Alkylating.tsv"),
          nonenriched = file.path(results_dir, "multi_spectrum",
                                  "multi_spectrum_NONenriched.Alkylating.tsv")
        ),
        signature_file = file.path(results_dir, "expected",
                                   "Signatures.Alkylating.expected_triallelic.csv"),
        signature_cols = c(
          "cyclophosphamide_557117b73fe2",
          "38_SBS2_0.996907_1",
          "20_SBS13_0.948838_1"
        ),
        signature_labels = c(
          "cyclophosphamide_557117b73fe2" = "cyclophosphamide",
          "38_SBS2_0.996907_1" = "SBS2",
          "20_SBS13_0.948838_1" = "SBS13"
        ),
        expected_labels = c(
          enriched = "MAV enriched\nexpected",
          nonenriched = "MAV non-enriched\nexpected"
        ),
        observed_labels = c(
          enriched = "MAV\nenriched samples",
          nonenriched = "MAV\nnon-enriched samples"
        ),
        plot_file = file.path(plots_dir,
                              "obs_vs_exp_multi.Alkylating.jpeg")
      )
    ))
  }

  if (mode == "Platinum") {

    platinum_splits <- c(
      "sigfilt_SBS17b_prop_lt_0.1",
      "sigfilt_SBS17b_prop_gt_0.1"
    )

    cfgs <- lapply(platinum_splits, function(split_name) {
      list(
        name = paste0("Platinum_", split_name),
        expected_files = list(
          enriched = file.path(results_dir, "expected",
                               paste0("enriched.Platinum.", split_name, ".expected_triallelic.csv")),
          nonenriched = file.path(results_dir, "expected",
                                  paste0("NONenriched.Platinum.", split_name, ".expected_triallelic.csv"))
        ),
        observed_files = list(
          enriched = file.path(results_dir, "multi_spectrum",
                               paste0("multi_spectrum_enriched.Platinum.", split_name, ".tsv")),
          nonenriched = file.path(results_dir, "multi_spectrum",
                                  paste0("multi_spectrum_NONenriched.Platinum.", split_name, ".tsv"))
        ),
        signature_file = file.path(results_dir, "expected",
                                   "Signatures.Platinum.expected_triallelic.csv"),
        signature_cols = c(
          "25_1",
          "37_1",
          "21_SBS31_0.953955_1",
          "14_1"
        ),
        signature_labels = c(
          "25_1" = "Carboplatin\n(E-SBS25)",
          "37_1" = "Oxaliplatin\n(E-SBS37)",
          "21_SBS31_0.953955_1" = "Carbo/Cisplatin\n(E-SBS21)",
          "14_1" = "Cis/Oxaliplatin\n(E-SBS14)"
        ),
        expected_labels = c(
          enriched = paste0("MAV enriched exp\n", split_name),
          nonenriched = paste0("MAV non-enriched exp\n", split_name)
        ),
        observed_labels = c(
          enriched = paste0("enriched samples\n", split_name),
          nonenriched = paste0("non-enriched samples\n", split_name)
        ),
        plot_file = file.path(plots_dir,
                              paste0("obs_vs_exp_multi.", split_name, ".jpeg"))
      )
    })

    return(cfgs)
  }

  stop("Unsupported mode")
}

read_signature_table <- function(path, wanted_cols) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  print(head(df))
  if ("context" %in% colnames(df)) {
    context_col <- "context"
  } else if ("X" %in% colnames(df)) {
    context_col <- "X"
  } else {
    context_col <- colnames(df)[1]
  }

  df %>%
    rename(context = all_of(context_col)) %>%
    select(context, all_of(wanted_cols))
}

compute_cosine_matrix <- function(df, row_cols, col_cols) {
  mat <- matrix(
    NA_real_,
    nrow = length(row_cols),
    ncol = length(col_cols),
    dimnames = list(row_cols, col_cols)
  )

  for (i in seq_along(row_cols)) {
    for (j in seq_along(col_cols)) {
      mat[i, j] <- cosine(df[[row_cols[i]]], df[[col_cols[j]]])
    }
  }

  mat
}

plot_heatmap <- function(mat, outfile) {
  if (!dir.exists(dirname(outfile))) {
    dir.create(dirname(outfile), recursive = TRUE)
  }

  jpeg(filename = outfile, width = 20, height = 9, res = 300, units = "cm")

  p <- pheatmap(
    mat,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    clustering_method = "complete",
    color = colorRampPalette(c("grey", "aquamarine4"))(50),
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 10,
    number_color = "black",
    angle_col = 45,
    cellwidth = 50,
    cellheight = 50,
    silent = TRUE
  )

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    3, 3,
    widths = unit.c(unit(1, "cm"), unit(1, "null"), unit(1, "cm")),
    heights = unit.c(unit(1, "cm"), unit(1, "null"), unit(1, "cm"))
  )))

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid.draw(p$gtable)
  popViewport()

  grid.text(
    "cosine",
    x = unit(0.84, "npc"),
    y = unit(0.95, "npc"),
    gp = gpar(fontsize = 12, fontface = "bold")
  )

  dev.off()
}