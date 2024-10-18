#!/bin/bash

num_sim=100
subset=1

Rscript src/R_scripts/Data_generation.R $PWD $num_sim
python3 src/python_scripts/2d_non-Gaussian-cov.py $num_sim
python3 src/python_scripts/2d_nonstationary.py $num_sim

Rscript src/R_scripts/Gaussian_kriging_with_covariates.R $PWD $num_sim $subset
Rscript src/R_scripts/Gaussian_kriging_nonstat.R $PWD $num_sim $subset
Rscript src/R_scripts/calculate_validation_results_GKrigning.R $PWD $num_sim
Rscript src/R_scripts/plot_results_reproducible.R $PWD
