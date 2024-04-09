#!/bin/bash

Rscript R_scripts/Data_generation.R "/home/nagp/Desktop/Bivariate_DeepKriging/" 3
python3 python_scripts/2d_non-Gaussian-cov.py 3
python3 python_scripts/2d_nonstationary.py 3
python3 python_scripts/pred_interval_nonGaussian-cov.py 3
python3 python_scripts/pred_interval_nonstationary-cov.py 3
