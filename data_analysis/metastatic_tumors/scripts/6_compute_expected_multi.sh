#This script contains examples of configurations to run 6_compute_expected_triallelic.py
#used in the analysis

# To calculate expected tri-allelic spectra based on the spectra of bi-allelic sites
python 6_compute_expected_multi.py samples --run-all

#To calculate expected tri-allelic spectra based on the spectra of known treatment signatures
python 6_compute_expected_multi.py signatures \
    --treatment Platinum 