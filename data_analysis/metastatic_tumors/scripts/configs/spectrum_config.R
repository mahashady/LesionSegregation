get_spectrum_config <- function(
  therapy_subset = "Alkylating",
  mode = c("biallelic", "multiallelic"),
  root_folder = "/home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors"
) {
  mode <- match.arg(mode)

  cfg <- list(
    therapy_subset = therapy_subset,
    mode = mode,

    # ---------------------------
    # Input metadata
    # ---------------------------
    enriched_file = paste0(root_folder, "/results/enrichment/all_enriched.chemo_alkyl_immuno.txt"),
    non_enriched_file = paste0(root_folder, "/results/enrichment/all_NONenriched.chemo_alkyl_immuno.txt"),

    metadata_sep = "\t",

    patient_col = "patientIdentifier",
    sample_col = "sample",

    # Keep only NONenriched samples with n_multi > 0
    require_n_multi_non_enriched = TRUE,
    n_multi_col = "n_multi",
    n_multi_min = 0,

    # ---------------------------
    # Genome composition
    # ---------------------------
    genome_counts_file = paste0(root_folder, "/data/genome_counts_tribases.json"),

    # ---------------------------
    # Optional signature exposure filtering
    # ---------------------------
    signature_filter = list(
      enabled = FALSE,
      file = paste0(root_folder, "/data/Pan_full.exposures.tsv"),
      sep = "",
      sbs31_col = "21_SBS31_0.953955_1",
      sbs17b_col = "19_SBS17b_0.932022_1",
      prop_col = "SBS17b_prop",
      cutoff = 0.1,
      direction = "lt"   # "lt" or "gt"    
      ),

    # ---------------------------
    # Optional export of selected sample lists
    # ---------------------------
    export_sample_lists_json = TRUE,
    sample_lists_json = file.path(paste0(root_folder, "/results/samples_lists/samples_lists.", therapy_subset, ".json")),

    # ---------------------------
    # Plot settings
    # ---------------------------
    plot_width_cm = 30,
    plot_height_cm = 20,
    plot_dpi = 300
  )

  cfg$modes <- list(
    biallelic = list(
      input_dir_enriched = paste0(root_folder, "/results/bi-allelic_sites/sites/bi_sites_by_sample_enriched.chemo_alkyl_immuno"),
      input_dir_non_enriched = paste0(root_folder, "/results/bi-allelic_sites/sites/bi_sites_by_sample_NONenriched.chemo_alkyl_immuno"),

      sample_input_suffix = "_bi_sites.txt",

      per_sample_output_dir_enriched = paste0(root_folder, "/results/bi-allelic_sites/spectra/bi_spectrum_by_sample_enriched.chemo_alkyl_immuno"),
      per_sample_output_dir_non_enriched = paste0(root_folder, "/results/bi-allelic_sites/spectra/bi_spectrum_by_sample_NONenriched.chemo_alkyl_immuno"),
      per_sample_output_suffix = "_bi_spectrum.txt",

      combined_output_enriched = paste0(root_folder, "/results/bi-allelic_sites/spectra/bi_spectrum/bi_spectrum_enriched.", therapy_subset, ".tsv"),
      combined_output_non_enriched = paste0(root_folder, "/results/bi-allelic_sites/spectra/bi-spectrum/bi_spectrum_NONenriched.", therapy_subset, ".tsv"),

      plot_file = paste0(root_folder, "/plots/biallelic_spectrum_enriched_vs_NONenriched_", therapy_subset, ".jpeg"),

      # file layout assumptions
      expected_ncol = 6,
      context_col = 5,
      alt_cols = c(6)
    ),

    multiallelic = list(
      input_dir_enriched = paste0(root_folder, "/results/multi-allelic_sites/sites/multi_sites_by_sample_enriched.chemo_alkyl_immuno"),
      input_dir_non_enriched = paste0(root_folder, "/results/multi-allelic_sites/sites/multi_sites_by_sample_NONenriched.chemo_alkyl_immuno"),

      sample_input_suffix = "_multi_sites.txt",

      per_sample_output_dir_enriched = paste0(root_folder, "/results/multi-allelic_sites/spectra/multi_spectrum_by_sample_enriched.chemo_alkyl_immuno"),
      per_sample_output_dir_non_enriched = paste0(root_folder, "/results/multi-allelic_sites/spectra/multi_spectrum_by_sample_NONenriched.chemo_alkyl_immuno"),
      per_sample_output_suffix = "_multi_spectrum.txt",

      combined_output_enriched = paste0(root_folder, "/results/multi-allelic_sites/spectra/multi_spectrum/multi_spectrum_enriched.", therapy_subset, ".tsv"),
      combined_output_non_enriched = paste0(root_folder, "/results/multi-allelic_sites/spectra/multi_spectrum/multi_spectrum_NONenriched.", therapy_subset, ".tsv"),

      plot_file = paste0(root_folder, "/plots/multiallelic_spectrum_enriched_vs_NONenriched_", therapy_subset, ".jpeg"),

      # file layout assumptions
      expected_ncol = 8,
      context_col = 5,
      alt_cols = c(6, 7)
    )
  )

  cfg
}