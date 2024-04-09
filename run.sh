#!/bin/bash

num_sim=100

Rscript src/R_scripts/Data_generation.R "/home/nagp/Desktop/Bivariate_DeepKriging/" $num_sim
echo ""
echo "######################   Biv.DeepKriging both simulations   ##################################"
echo ""
python3 src/python_scripts/2d_non-Gaussian-cov.py $num_sim
python3 src/python_scripts/2d_nonstationary.py $num_sim
python3 src/python_scripts/pred_interval_nonGaussian-cov.py $num_sim
python3 src/python_scripts/pred_interval_nonstationary-cov.py $num_sim

echo ""
echo "######################   Cokriging non-Gaussian   ##################################"
echo ""

Rscript src/R_scripts/Gaussian_kriging_with_covariates.R "/home/nagp/Desktop/Bivariate_DeepKriging/" $num_sim

echo ""
echo "######################   Cokriging nonstationary   ##################################"
echo ""

Rscript src/R_scripts/Gaussian_kriging_nonstat.R "/home/nagp/Desktop/Bivariate_DeepKriging/" $num_sim

Rscript src/R_scripts/calculate_validation_results_GKrigning.R "/home/nagp/Desktop/Bivariate_DeepKriging/" $num_sim
