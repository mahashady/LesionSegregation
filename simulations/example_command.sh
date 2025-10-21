# Initialize cell
python simulations/initialize_cell.py --working_dir ~/LesionSegregation/simulations/  --chr_file chr_length_files/grcm38.p6_rounded.txt  --driver_freq 7  --outprefix test_run

# Draw lesions: (assumes provided cell is present in directory, can change --cell_dir argument passed to use a different initialized cell 
python simulations/draw_lesions.py  --working_dir ~/LesionSegregation/simulations/ --cell_dir cells/grcm38.p6_rounded_driver_loci_freq_5/2024-06-02_13-23/ --lesion_counts 2,0 --drivers_same_chrom False --drivers_same_strand False --driver_s_first 0.4 --driver_s_max 0.5 --outsuffix s_0.5_0.4

# Run simulations
python simulations/simulations.py --working_dir ~/LesionSegregation/simulations/ --cell_dir cells/grcm38.p6_rounded_driver_loci_freq_5/2024-06-02_13-23/ --lesions_dir cells/grcm38.p6_rounded_driver_loci_freq_5/2024-06-02_13-23/lesions_2_drivers_0_other_s_0.5_0.4/ --epistasis_type add --driver_strand tx --two_strands_required True --error_free_replication 0.25 --error_free_transcription 0 --repair_rate 0.25 --shared_mutations_through_repair 0.05 --repair_process lesion_mut --n_track_segregations 5 --min_detailed_generations 10 --max_generations 300 --carrying_capacity 1e8  --min_tumor_size 1e6 --outprefix test_sim  --seed_min 0 --seed_max 100000 --sim_count 100000 --verbose False

