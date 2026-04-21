source("configs/spectrum_config.R")
source("utils/spectrum_utils.R")

# =========================================================
# User parameters
# =========================================================
therapy_subset <- "Alkylating" 
#therapy_subset <- "Platinum"

mode <- "biallelic"
#mode <- "multiallelic"

cfg <- get_spectrum_config(
  therapy_subset = therapy_subset,
  mode = mode
)

# =========================================================
# Optional overrides
# =========================================================

# If your metadata files are comma-separated instead of tab-separated:
# cfg$metadata_sep <- ","

# If you want the multiallelic script behavior with signature filtering:
#cfg$signature_filter$enabled <- TRUE
#cfg$signature_filter$prop_col <- "SBS17b_prop"
#cfg$signature_filter$cutoff <- 0.1
#fg$signature_filter$direction <- "lt"

# If you do not want JSON export:
# cfg$export_sample_lists_json <- FALSE

# =========================================================
# Construct output names
# =========================================================
sig_suffix <- make_signature_suffix(cfg)
mode_cfg <- cfg$modes[[cfg$mode]]

sample_lists_json <- append_suffix_before_ext(cfg$sample_lists_json, sig_suffix)
combined_output_enriched <- append_suffix_before_ext(mode_cfg$combined_output_enriched, sig_suffix)
combined_output_non_enriched <- append_suffix_before_ext(mode_cfg$combined_output_non_enriched, sig_suffix)
plot_file <- append_suffix_before_ext(mode_cfg$plot_file, sig_suffix)

message("Running spectrum analysis")
message("Mode: ", cfg$mode)
message("Therapy subset: ", cfg$therapy_subset)

# =========================================================
# Load shared inputs
# =========================================================
genome_composition <- load_genome_composition(cfg$genome_counts_file)
possible_mutations <- get_possible_mutations()

sample_lists <- get_sample_lists(cfg)
sig_suffix <- make_signature_suffix(cfg)

message("Number of enriched samples = ", length(sample_lists$enriched_list))
message("Number of NONenriched samples = ", length(sample_lists$non_enriched_list))

if (isTRUE(cfg$export_sample_lists_json)) {
  export_sample_lists_json(
    enriched_list = sample_lists$enriched_list,
    non_enriched_list = sample_lists$non_enriched_list,
    outfile = sample_lists_json <- append_suffix_before_ext(cfg$sample_lists_json, sig_suffix)
  )
}

# =========================================================
# Process enriched samples
# =========================================================
all_mutations_enriched <- process_samples(
  sample_ids = sample_lists$enriched_list,
  input_dir = mode_cfg$input_dir_enriched,
  input_suffix = mode_cfg$sample_input_suffix,
  per_sample_output_dir = mode_cfg$per_sample_output_dir_enriched,
  per_sample_output_suffix = mode_cfg$per_sample_output_suffix,
  context_col = mode_cfg$context_col,
  alt_cols = mode_cfg$alt_cols,
  expected_ncol = mode_cfg$expected_ncol,
  possible_mutations = possible_mutations,
  genome_composition = genome_composition,
  verbose = TRUE
)

message("Enriched done")

# =========================================================
# Process NONenriched samples
# =========================================================
all_mutations_non_enriched <- process_samples(
  sample_ids = sample_lists$non_enriched_list,
  input_dir = mode_cfg$input_dir_non_enriched,
  input_suffix = mode_cfg$sample_input_suffix,
  per_sample_output_dir = mode_cfg$per_sample_output_dir_non_enriched,
  per_sample_output_suffix = mode_cfg$per_sample_output_suffix,
  context_col = mode_cfg$context_col,
  alt_cols = mode_cfg$alt_cols,
  expected_ncol = mode_cfg$expected_ncol,
  possible_mutations = possible_mutations,
  genome_composition = genome_composition,
  verbose = TRUE
)

message("NON enriched done")

# =========================================================
# Build combined result
# =========================================================
result <- build_enriched_vs_nonenriched_result(
  enriched_mutations = all_mutations_enriched,
  non_enriched_mutations = all_mutations_non_enriched,
  possible_mutations = possible_mutations
)

result$N_muts <- as.numeric(result$N_muts)

result_enriched <- subset(result, subset == "subset_enriched")
result_non_enriched <- subset(result, subset == "subset_non_enriched")

# add rates for saved outputs
result_enriched <- add_mutation_rates(result_enriched, genome_composition)
result_non_enriched <- add_mutation_rates(result_non_enriched, genome_composition)

# =========================================================
# Save combined outputs
# =========================================================
combined_output_enriched <- sub(
  "\\.tsv$",
  paste0(sig_suffix, ".tsv"),
  mode_cfg$combined_output_enriched
)

combined_output_non_enriched <- sub(
  "\\.tsv$",
  paste0(sig_suffix, ".tsv"),
  mode_cfg$combined_output_non_enriched
)

plot_file <- sub(
  "\\.jpeg$",
  paste0(sig_suffix, ".jpeg"),
  mode_cfg$plot_file
)

dir.create(dirname(combined_output_enriched), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(combined_output_non_enriched), recursive = TRUE, showWarnings = FALSE)

write.csv(result_enriched, file = combined_output_enriched, row.names = FALSE)
write.csv(result_non_enriched, file = combined_output_non_enriched, row.names = FALSE)

# =========================================================
# Plot
# =========================================================
p <- plot_spectrum_comparison(
  result = result,
  title = paste(cfg$therapy_subset, "-", cfg$mode),
  outfile = plot_file,
  width = cfg$plot_width_cm,
  height = cfg$plot_height_cm,
  dpi = cfg$plot_dpi
)

print(p)
message("Saved plot to: ", mode_cfg$plot_file)