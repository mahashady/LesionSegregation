library(ggplot2)
library(jsonlite)

# =========================================================
# Constants
# =========================================================

get_muttype_levels <- function() {
  c(
    "A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T",
    "G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
    "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T",
    "G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
    "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T",
    "G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
    "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
    "G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
    "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T",
    "G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
    "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T",
    "G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"
  )
}

get_possible_mutations <- function(muttype_level = get_muttype_levels()) {
  paste0(
    substr(muttype_level, 1, 1),
    substr(muttype_level, 3, 3),
    substr(muttype_level, 7, 7),
    ">",
    substr(muttype_level, 5, 5)
  )
}

get_x_labels <- function() {
  rep(
    c("A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T",
      "G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T"),
    6
  )
}

get_complement_map <- function() {
  c(A = "T", T = "A", G = "C", C = "G")
}

# =========================================================
# Core helpers
# =========================================================

load_genome_composition <- function(json_file) {
  if (!file.exists(json_file)) {
    stop("File not found: ", json_file)
  }

  json_data <- jsonlite::fromJSON(json_file)

  data.frame(
    context = names(json_data),
    counts = as.numeric(unlist(json_data)),
    stringsAsFactors = FALSE
  )
}

reverse_compl <- function(context, complement_map = get_complement_map()) {
  stopifnot(nchar(context) == 3)
  paste0(
    complement_map[substr(context, 3, 3)],
    complement_map[substr(context, 2, 2)],
    complement_map[substr(context, 1, 1)]
  )
}

normalize_mutation_96 <- function(context,
                                  alt_nucleotide,
                                  pyrimidines = c("T", "C"),
                                  complement_map = get_complement_map()) {
  context <- as.character(context)
  alt_nucleotide <- as.character(alt_nucleotide)

  if (is.na(context) || is.na(alt_nucleotide)) {
    return(NA_character_)
  }

  if (nchar(context) != 3 || nchar(alt_nucleotide) < 1) {
    return(NA_character_)
  }

  alt_nucleotide <- substr(alt_nucleotide, 1, 1)
  center_base <- substr(context, 2, 2)

  if (center_base %in% pyrimidines) {
    return(paste0(context, ">", alt_nucleotide))
  }

  rev_context <- reverse_compl(context, complement_map = complement_map)
  paste0(rev_context, ">", complement_map[alt_nucleotide])
}

# =========================================================
# Sample metadata
# =========================================================

load_signature_exposures <- function(sig_cfg) {
  sig_exposures <- read.table(sig_cfg$file, header = TRUE)
  t_sig_exposures <- as.data.frame(t(sig_exposures), stringsAsFactors = FALSE)

  t_sig_exposures$SBS31_prop <- t_sig_exposures[[sig_cfg$sbs31_col]] / rowSums(t_sig_exposures)
  t_sig_exposures$SBS17b_prop <- t_sig_exposures[[sig_cfg$sbs17b_col]] / rowSums(t_sig_exposures)
  t_sig_exposures$sample <- rownames(t_sig_exposures)

  t_sig_exposures
}

get_sample_lists <- function(cfg) {
  df_enriched <- read.table(
    cfg$enriched_file,
    sep = cfg$metadata_sep,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  df_non_enriched <- read.table(
    cfg$non_enriched_file,
    sep = cfg$metadata_sep,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  if (isTRUE(cfg$require_n_multi_non_enriched)) {
    df_non_enriched <- df_non_enriched[
      df_non_enriched[[cfg$n_multi_col]] > cfg$n_multi_min,
      ,
      drop = FALSE
    ]
  }

  if (isTRUE(cfg$signature_filter$enabled)) {
    sig_df <- load_signature_exposures(cfg$signature_filter)

    if (!cfg$sample_col %in% names(df_enriched)) {
      stop("sample_col not found in enriched metadata: ", cfg$sample_col)
    }
    if (!cfg$sample_col %in% names(df_non_enriched)) {
      stop("sample_col not found in non-enriched metadata: ", cfg$sample_col)
    }

    df_enriched <- merge(df_enriched, sig_df, by = cfg$sample_col)
    df_non_enriched <- merge(df_non_enriched, sig_df, by = cfg$sample_col)

    prop_col <- cfg$signature_filter$prop_col
    cutoff <- cfg$signature_filter$cutoff
    direction <- cfg$signature_filter$direction

    keep_fun <- switch(
      direction,
      lt = function(x) x < cutoff,
      le = function(x) x <= cutoff,
      gt = function(x) x > cutoff,
      ge = function(x) x >= cutoff,
      eq = function(x) x == cutoff,
      stop("Unsupported signature_filter$direction: ", direction)
    )

    df_enriched <- df_enriched[
      keep_fun(df_enriched[[prop_col]]),
      ,
      drop = FALSE
    ]
    df_non_enriched <- df_non_enriched[
      keep_fun(df_non_enriched[[prop_col]]),
      ,
      drop = FALSE
    ]
  }
  
  enriched_list <- df_enriched[df_enriched[[cfg$therapy_subset]] == 1, cfg$patient_col]
  non_enriched_list <- df_non_enriched[df_non_enriched[[cfg$therapy_subset]] == 1, cfg$patient_col]

  list(
    df_enriched = df_enriched,
    df_non_enriched = df_non_enriched,
    enriched_list = enriched_list,
    non_enriched_list = non_enriched_list
  )
}

export_sample_lists_json <- function(enriched_list, non_enriched_list, outfile) {
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  json_samples_list <- list(
    enriched = enriched_list,
    NONenriched = non_enriched_list
  )
  write_json(json_samples_list, outfile, pretty = TRUE)
}

# =========================================================
# Mutation extraction
# =========================================================

extract_mutations_from_df <- function(df_sample,
                                      context_col,
                                      alt_cols,
                                      expected_ncol = NULL) {
  if (!is.null(expected_ncol) && ncol(df_sample) != expected_ncol) {
    return(character(0))
  }

  if (nrow(df_sample) == 0) {
    return(character(0))
  }
  context <- as.character(df_sample[[context_col]])
  all_mutations <- character(0)
  #print(context)
  for (alt_col in alt_cols) {
    alt_vals <- substr(as.character(df_sample[[alt_col]]), 1, 1)
    #print(alt_vals)
    muts <- mapply(
      normalize_mutation_96,
      context = context,
      alt_nucleotide = alt_vals,
      USE.NAMES = FALSE
    )
    all_mutations <- c(all_mutations, muts[!is.na(muts)])
  }

  all_mutations
}

# =========================================================
# Spectrum helpers
# =========================================================

mutations_to_spectrum <- function(mut_vec,
                                  possible_mutations = get_possible_mutations(),
                                  muttype_level = get_muttype_levels()) {
  if (length(mut_vec) == 0) {
    mut_df <- data.frame(
      mutation = possible_mutations,
      N_muts = 0,
      stringsAsFactors = FALSE
    )
  } else {
    tab <- table(mut_vec)
    mut_df <- data.frame(
      mutation = names(tab),
      N_muts = as.numeric(tab),
      stringsAsFactors = FALSE
    )

    full_df <- data.frame(
      mutation = possible_mutations,
      stringsAsFactors = FALSE
    )

    mut_df <- merge(full_df, mut_df, by = "mutation", all.x = TRUE, sort = FALSE)
    mut_df$N_muts[is.na(mut_df$N_muts)] <- 0
  }

  mut_df$context <- substr(mut_df$mutation, 1, 3)
  mut_df$nucleotide <- paste0(substr(mut_df$mutation, 2, 2), ">", substr(mut_df$mutation, 5, 5))
  mut_df$Type <- paste0(
    substr(mut_df$context, 1, 1),
    "[",
    mut_df$nucleotide,
    "]",
    substr(mut_df$context, 3, 3)
  )

  mut_df$Type <- factor(mut_df$Type, levels = muttype_level)
  mut_df
}

add_mutation_rates <- function(spectrum_df, genome_composition) {
  out <- merge(spectrum_df, genome_composition, by = "context", all.x = TRUE, sort = FALSE)
  out$Mutrate <- as.numeric(out$N_muts) / as.numeric(out$counts)
  out
}

# =========================================================
# Sample processing
# =========================================================

process_samples <- function(sample_ids,
                            input_dir,
                            input_suffix,
                            per_sample_output_dir = NULL,
                            per_sample_output_suffix = "_spectrum.txt",
                            context_col,
                            alt_cols,
                            expected_ncol,
                            possible_mutations,
                            genome_composition,
                            verbose = TRUE) {
  all_mutations <- character(0)

  if (!is.null(per_sample_output_dir)) {
    dir.create(per_sample_output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  for (sample_id in sample_ids) {
    print(sample_id)
    file_name <- file.path(input_dir, paste0(sample_id, input_suffix))
    print(file_name)
    if (!file.exists(file_name)) {
      next
    }

    if (verbose) {
      message("Processing sample: ", sample_id)
    }

    df_sample <- read.table(file_name, stringsAsFactors = FALSE, sep="\t")
    print(head(df_sample))

    sample_mutations <- extract_mutations_from_df(
      df_sample = df_sample,
      context_col = context_col,
      alt_cols = alt_cols,
      expected_ncol = expected_ncol
    )

    all_mutations <- c(all_mutations, sample_mutations)

    if (!is.null(per_sample_output_dir)) {
      df_sample_mutations <- mutations_to_spectrum(
        sample_mutations,
        possible_mutations = possible_mutations
      )
      df_sample_mutations <- add_mutation_rates(df_sample_mutations, genome_composition)

      write.csv(
        df_sample_mutations,
        file = file.path(per_sample_output_dir, paste0(sample_id, per_sample_output_suffix)),
        row.names = FALSE
      )
    }
  }

  all_mutations
}

build_enriched_vs_nonenriched_result <- function(enriched_mutations,
                                                 non_enriched_mutations,
                                                 possible_mutations = get_possible_mutations()) {
  result_enriched <- mutations_to_spectrum(
    enriched_mutations,
    possible_mutations = possible_mutations
  )
  result_enriched$subset <- "subset_enriched"

  result_non_enriched <- mutations_to_spectrum(
    non_enriched_mutations,
    possible_mutations = possible_mutations
  )
  result_non_enriched$subset <- "subset_non_enriched"

  rbind(result_enriched, result_non_enriched)
}

# =========================================================
# Plotting
# =========================================================

plot_spectrum_comparison <- function(result,
                                     title = NULL,
                                     outfile = NULL,
                                     width = 30,
                                     height = 20,
                                     units = "cm",
                                     dpi = 300,
                                     muttype_level = get_muttype_levels(),
                                     x_labels = get_x_labels()) {
  result$Type <- factor(result$Type, levels = muttype_level)

  p <- ggplot(result) +
    geom_bar(
      aes(
        x = factor(Type, levels = muttype_level),
        y = N_muts,
        fill = nucleotide
      ),
      stat = "identity"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=10),
      legend.position = "top",
      legend.box = "horizontal",
      strip.text = element_text(size = 10),
      legend.key.size = grid::unit(0.4, "cm"),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_discrete(labels = x_labels) +
    scale_fill_manual(values = c(
      "C>A" = "deepskyblue2",
      "C>G" = "black",
      "C>T" = "red",
      "T>A" = "grey",
      "T>C" = "darkolivegreen3",
      "T>G" = "lightpink"
    )) +
    xlab("") +
    facet_wrap(~subset, nrow = 2)

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  if (!is.null(outfile)) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    ggsave(outfile, p, width = width, height = height, dpi = dpi, units = units)
  }

  p
}

make_signature_suffix <- function(cfg) {
  if (!isTRUE(cfg$signature_filter$enabled)) {
    return("")
  }
  
  direction <- cfg$signature_filter$direction
  
  cutoff <- if (!is.null(cfg$signature_filter$cutoff)) {
    cfg$signature_filter$cutoff
  } else if (!is.null(cfg$signature_filter$max_prop)) {
    cfg$signature_filter$max_prop
  } else {
    "NA"
  }
  
  prop_col <- cfg$signature_filter$prop_col
  
  paste0(".sigfilt_", prop_col, "_", direction, "_", cutoff)
}

append_suffix_before_ext <- function(path, suffix) {
  sub("(\\.[^.]+)$", paste0(suffix, "\\1"), path)
}