library(tidyverse)

# ----------------------------
# Config
# ----------------------------
chr_size_file <- "../../data/c3h.chromSize.main"

hmm_multi_dir <- "../results/HMM_multi_result"
hmm_as_dir <- "../results/HMM_as_result"

out_dir <- "../results/HMMmulti_vs_HMM_as_mixtures"
subdirs <- c("plots", "tables")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots", "per_sample"), showWarnings = FALSE, recursive = TRUE)

n_perm <- 1
alpha <- 0.05 #significance pval threshold
set.seed(42)

# ----------------------------
# Read chromosome sizes
# ----------------------------

read_chr_sizes <- function(chr_size_file) {
  chr_size <- read.table(chr_size_file, header = FALSE)
  colnames(chr_size) <- c("chr", "size")

  chr_size %>%
    mutate(chr = as.character(chr))
}

# ----------------------------
# Read WC asymmetry HMM states
# ----------------------------

read_wc_asymmetry <- function(sample, hmm_as_dir, chr_size) {
  file <- file.path(hmm_as_dir, paste0(sample, "_intervals_HMMas.bed"))
  df <- read.table(file, sep = "\t", header = TRUE)

  df %>%
    mutate(
      chr = as.character(chr),
    ) %>%
    left_join(chr_size, by = "chr") %>%
    mutate(
      start_pos_relative = start_pos / size,
      end_pos_relative = end_pos / size
    )
}

# ----------------------------
# Read multiallelic HMM states
# ----------------------------

read_multiallelic_hmm <- function(sample, hmm_multi_dir, chr_size) {
  file <- file.path(hmm_multi_dir, paste0(sample, "_intervals_HMMmulti.bed"))
  df <- read.table(file, header = TRUE)

  df %>%
    mutate(
      chr = as.character(chr),
    ) %>%
    left_join(chr_size, by = "chr") %>%
    mutate(
      start_pos_relative = start_pos / size,
      end_pos_relative = end_pos / size
    )
}

# ----------------------------
# Extract interval boundaries
# ----------------------------

get_boundaries <- function(df) {
  df %>%
    select(chr, start_pos_relative, end_pos_relative) %>%
    pivot_longer(
      cols = c(start_pos_relative, end_pos_relative),
      names_to = "boundary_type",
      values_to = "pos_relative"
    ) %>%
    filter(!is.na(pos_relative)) %>%
    filter(pos_relative > 0, pos_relative < 1) %>%
    distinct(chr, pos_relative)
}

# ----------------------------
# Observed distances
# ----------------------------

calculate_observed_distances <- function(df_as, df_multi) {
  as_boundaries <- get_boundaries(df_as)
  multi_boundaries <- get_boundaries(df_multi)

  chromosomes <- intersect(
    unique(as_boundaries$chr),
    unique(multi_boundaries$chr)
  )

  result <- list()

  for (chr_i in chromosomes) {
    as_pos <- sort(as_boundaries$pos_relative[as_boundaries$chr == chr_i])
    multi_pos <- sort(multi_boundaries$pos_relative[multi_boundaries$chr == chr_i])

    if (length(as_pos) == 0 || length(multi_pos) == 0) {
      next
    }

    result[[chr_i]] <- data.frame(
      chr = chr_i,
      multi_boundary = multi_pos,
      min_distance = sapply(multi_pos, function(x) min(abs(as_pos - x)))
    )
  }

  bind_rows(result)
}

# ----------------------------
# Permutation distances
# ----------------------------

calculate_permutation_distances <- function(
    df_as,
    df_multi,
    n_perm = 1000,
    exclude_same_chr = TRUE
) {
  as_boundaries <- get_boundaries(df_as)
  multi_boundaries <- get_boundaries(df_multi)

  chromosomes <- intersect(
    unique(as_boundaries$chr),
    unique(multi_boundaries$chr)
  )

  result <- list()
  counter <- 1

  for (perm_i in seq_len(n_perm)) {
    for (chr_i in chromosomes) {
      possible_chr <- chromosomes

      if (exclude_same_chr) {
        possible_chr <- setdiff(possible_chr, chr_i)
      }

      if (length(possible_chr) == 0) {
        next
      }

      random_chr <- sample(possible_chr, 1)

      as_pos <- sort(as_boundaries$pos_relative[as_boundaries$chr == random_chr])
      multi_pos <- sort(multi_boundaries$pos_relative[multi_boundaries$chr == chr_i])

      if (length(as_pos) == 0 || length(multi_pos) == 0) {
        next
      }

      result[[counter]] <- data.frame(
        permutation = perm_i,
        chr = chr_i,
        permuted_as_chr = random_chr,
        multi_boundary = multi_pos,
        min_distance = sapply(multi_pos, function(x) min(abs(as_pos - x)))
      )

      counter <- counter + 1
    }
  }

  bind_rows(result)
}
# ----------------------------
# Plot tracks for one sample
# ----------------------------
plot_one_sample_hmm_tracks <- function(
    sample,
    chr_size,
    hmm_as_dir,
    hmm_multi_dir,
    out_dir,
    chr_levels = as.character(1:19)
) {
  df_as <- read_wc_asymmetry(sample, hmm_as_dir, chr_size) %>%
    select(chr, start_pos, end_pos, HMM_as_state) %>%
    rename("value"="HMM_as_state") %>%
    mutate(HMM = "WCas")

  df_multi <- read_multiallelic_hmm(sample, hmm_multi_dir, chr_size) %>%
    select(chr, start_pos, end_pos, HMM_multi_state) %>%
    rename("value"="HMM_multi_state") %>%
    mutate(HMM = "multi")

  result <- bind_rows(df_as, df_multi) %>%
    filter(chr != "X") %>%
    mutate(
      chr = factor(chr, levels = chr_levels),
    )
  print(head(result))
  p <- ggplot(result) +
    geom_segment(
      aes(
        x = start_pos,
        y = value,
        xend = end_pos,
        yend = value
      )
    ) +
    #geom_blank(aes(y = y_max)) +
    facet_grid(
      rows = vars(HMM),
      cols = vars(chr),
      scales = "free"
    ) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    labs(
      title = sample,
      x = "position on chr",
      y = ""
    )

  ggsave(
    file.path(out_dir, paste0(sample, "_WCas_vs_multi_hmm.png")),
    p,
    width = 30,
    height = 10,
    units = "cm",
    dpi = 600
  )

  return(p)
}

# ----------------------------
# Analyze one sample
# ----------------------------

analyze_one_sample <- function(
    sample,
    chr_size,
    hmm_as_dir,
    hmm_multi_dir,
    out_dir,
    n_perm = 1000,
    alpha = 0.05
) {
  message("Processing sample: ", sample)

  df_as <- read_wc_asymmetry(sample, hmm_as_dir, chr_size)
  df_multi <- read_multiallelic_hmm(sample, hmm_multi_dir, chr_size)

  n_states_as <- length(unique(df_as$HMM_as_state))
  n_states_multi <- length(unique(df_multi$HMM_multi_state))
  message("n_states multi=", n_states_multi)
  message("n_states as=", n_states_as)
  if (n_states_as <= 1 || n_states_multi <= 1) {
    warning(
      paste0(
        "Skipping sample ", sample,
        " (insufficient state diversity: ",
        "WC_asymmetry=", n_states_as,
        ", multi=", n_states_multi, ")"
      )
    )

    summary_df <- data.frame(
      sample = sample,
      n_observed = NA,
      n_permutation = NA,
      mean_observed = NA,
      mean_permutation = NA,
      median_observed = NA,
      median_permutation = NA,
      ks_pval = NA,
      wilcox_pval = NA,
      direction = NA,
      significant_ks = NA,
      significant_wilcox = NA,
      n_states_as = n_states_as,
      n_states_multi = n_states_multi,
      skipped = TRUE
    )

    return(list(
      observed = data.frame(),
      permutation = data.frame(),
      summary = summary_df
    ))
  }

  obs <- calculate_observed_distances(df_as, df_multi) %>%
    mutate(sample = sample, type = "observed")

  perm <- calculate_permutation_distances(
    df_as = df_as,
    df_multi = df_multi,
    n_perm = n_perm
  ) %>%
    mutate(sample = sample, type = "permutation")

  if (nrow(obs) == 0 || nrow(perm) == 0) {
    warning("No valid distances for sample: ", sample)

    summary_df <- data.frame(
      sample = sample,
      n_observed = nrow(obs),
      n_permutation = nrow(perm),
      mean_observed = NA,
      mean_permutation = NA,
      median_observed = NA,
      median_permutation = NA,
      ks_pval = NA,
      wilcox_pval = NA,
      direction = NA,
      significant_ks = NA,
      significant_wilcox = NA,
      n_states_as = n_states_as,
      n_states_multi = n_states_multi,
      skipped = TRUE      
    )

    return(list(
      observed = obs,
      permutation = perm,
      summary = summary_df
    ))
  }

  ks_pval <- ks.test(obs$min_distance, perm$min_distance)$p.value

  wilcox_pval <- wilcox.test(
    obs$min_distance,
    perm$min_distance,
    alternative = "less"
  )$p.value

  mean_obs <- mean(obs$min_distance, na.rm = TRUE)
  mean_perm <- mean(perm$min_distance, na.rm = TRUE)

  direction <- ifelse(mean_obs < mean_perm, "closer", "farther_or_equal")

  summary_df <- data.frame(
    sample = sample,
    n_observed = nrow(obs),
    n_permutation = nrow(perm),
    mean_observed = mean_obs,
    mean_permutation = mean_perm,
    median_observed = median(obs$min_distance, na.rm = TRUE),
    median_permutation = median(perm$min_distance, na.rm = TRUE),
    ks_pval = ks_pval,
    wilcox_pval = wilcox_pval,
    direction = direction,
    significant_ks = ks_pval < alpha & direction == "closer",
    significant_wilcox = wilcox_pval < alpha & direction == "closer"
  )
  summary_df$n_states_as <- n_states_as
  summary_df$n_states_multi <- n_states_multi
  summary_df$skipped <- FALSE  

  # Per-sample histogram
  plot_df <- bind_rows(
    obs %>% select(min_distance) %>% mutate(type = "observed"),
    perm %>% select(min_distance) %>% mutate(type = "permutation")
  )

  y_max <- max(table(cut(plot_df$min_distance, breaks = 15)))

  p_hist <- ggplot(plot_df, aes(x = min_distance, fill = type)) +
    geom_histogram(alpha = 0.7, bins = 15, position = "dodge") +
    theme_bw() +
    ggtitle(sample) +
    scale_fill_manual(
      values = c("observed" = "darkred", "permutation" = "darkgrey")
    ) +
    annotate(
      "text",
      x = quantile(plot_df$min_distance, 0.7, na.rm = TRUE),
      y = y_max,
      label = paste("KS p =", signif(ks_pval, 3)),
      hjust = 0
    ) +
    ylab("n observations")

  dir.create(
    file.path(out_dir, "plots", "per_sample"),
    showWarnings = FALSE,
    recursive = TRUE
  )

  ggsave(
    file.path(out_dir, "plots", "per_sample", paste0(sample, "_distance_hist.png")),
    p_hist,
    width = 12,
    height = 10,
    units = "cm",
    dpi = 600
  )

  return(list(
    observed = obs,
    permutation = perm,
    summary = summary_df
  ))
}

# ----------------------------
# Run all samples
# ----------------------------

#df_samples <- read.table("../results/Summary_LAD_with_symmetrical_no_mixtures.txt", sep=",", header=TRUE)
df_samples <- read.table("../results/Summary_mixtures.txt", sep=",", header=TRUE)
samples <- df_samples$Sample
chr_size <- read_chr_sizes(chr_size_file)

plot_one_sample_hmm_tracks(
  sample = "91215_N2",
  chr_size = chr_size,
  hmm_as_dir = hmm_as_dir,
  hmm_multi_dir = hmm_multi_dir,
  out_dir = file.path(out_dir, "plots")
)

plot_one_sample_hmm_tracks(
  sample = "91863_N3",
  chr_size = chr_size,
  hmm_as_dir = hmm_as_dir,
  hmm_multi_dir = hmm_multi_dir,
  out_dir = file.path(out_dir, "plots")
)

plot_one_sample_hmm_tracks(
  sample = "88308_N2",
  chr_size = chr_size,
  hmm_as_dir = hmm_as_dir,
  hmm_multi_dir = hmm_multi_dir,
  out_dir = file.path(out_dir, "plots")
)

results <- lapply(
  samples,
  analyze_one_sample,
  chr_size = chr_size,
  hmm_as_dir = hmm_as_dir,
  hmm_multi_dir = hmm_multi_dir,
  out_dir = out_dir,
  n_perm = n_perm,
  alpha = alpha
)

obs_all <- bind_rows(lapply(results, function(x) x[["observed"]]))
perm_all <- bind_rows(lapply(results, function(x) x[["permutation"]]))
summary_all <- bind_rows(lapply(results, function(x) x[["summary"]]))
#summary_all <- summary_all %>% filter(!skipped)
# ----------------------------
# Meta-summary across samples
# ----------------------------

valid_summary <- summary_all %>%
  filter(!is.na(ks_pval))

n_samples <- nrow(valid_summary)

n_sig_ks <- sum(valid_summary$significant_ks, na.rm = TRUE)
n_sig_wilcox <- sum(valid_summary$significant_wilcox, na.rm = TRUE)

n_closer <- sum(valid_summary$direction == "closer", na.rm = TRUE)

binom_ks_pval <- binom.test(
  n_sig_ks,
  n_samples,
  p = alpha,
  alternative = "greater"
)$p.value

binom_wilcox_pval <- binom.test(
  n_sig_wilcox,
  n_samples,
  p = alpha,
  alternative = "greater"
)$p.value

direction_binom_pval <- binom.test(
  n_closer,
  n_samples,
  p = 0.5,
  alternative = "greater"
)$p.value

meta_summary <- data.frame(
  n_samples = n_samples,

  n_closer = n_closer,
  proportion_closer = n_closer / n_samples,
  direction_binom_pval = direction_binom_pval,

  n_significant_ks = n_sig_ks,
  proportion_significant_ks = n_sig_ks / n_samples,
  binomial_pval_ks = binom_ks_pval,

  n_significant_wilcox = n_sig_wilcox,
  proportion_significant_wilcox = n_sig_wilcox / n_samples,
  binomial_pval_wilcox = binom_wilcox_pval
)

meta_summary

# ----------------------------
# Save outputs
# ----------------------------

write.table(
  summary_all,
  file.path(out_dir, "per_sample_WCas_vs_multi_hmm_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  meta_summary,
  file.path(out_dir, "meta_summary_WCas_vs_multi_hmm.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  obs_all,
  file.path(out_dir, "observed_distances_WCas_vs_multi_hmm.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  perm_all,
  file.path(out_dir, "permutation_distances_WCas_vs_multi_hmm.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# Plot per-sample p-values
# ----------------------------
summary_all <- summary_all %>% filter(!skipped)
p1 <- ggplot(summary_all, aes(x = ks_pval)) +
  geom_histogram(bins = 20, fill = "grey70", color = "black") +
  geom_vline(xintercept = alpha, linetype = "dashed", color = "darkred") +
  scale_x_continuous(trans="log10") +
  theme_bw() +
  labs(
    title = "Per-sample KS p-values",
    x = "KS p-value",
    y = "Number of samples"
  )

ggsave(
  file.path(out_dir, "plots","per_sample_ks_pvalues_WCas_vs_multial_hmm.png"),
  p1,
  width = 12,
  height = 10,
  units = "cm",
  dpi = 600
)

# ----------------------------
# Plot effect sizes per sample
# ----------------------------

p2 <- ggplot(
  summary_all,
  aes(
    x = reorder(sample, mean_observed - mean_permutation),
    y = mean_observed - mean_permutation,
    fill = significant_ks
  )
) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(
    values = c("TRUE" = "darkred", "FALSE" = "grey70"),
    na.value = "grey90"
  ) +
  labs(
    title = "Per-sample co-occurrence effect",
    x = "Sample",
    y = "Difference between\nobserved and permuted min distance",
    fill = "Significant KS"
  )

ggsave(
  file.path(out_dir, "plots","per_sample_effect_WCas_vs_multial_hmm.png"),
  p2,
  width = 25,
  height = 10,
  units = "cm",
  dpi = 300
)

