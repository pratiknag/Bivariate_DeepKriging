This is the Readme file for reproducing the results when of the paper [Bivariate DeepKriging](https://arxiv.org/abs/2307.08038). The commands we have provided below. Running all of them sequentially will generate two plots in .pdf format inside the plot_results/ directory. These two plots are represented in the paper Fig 2,3. 

## Pre-requisites

Please ensure that you have R (>= 4.4.1) installed on your system, along with the following libraries:
```
geoR
MASS
fields
parallel
doSNOW
```
Additionally, please verify if Python 3+ is installed. If not please download and install python from [here](https://www.python.org/downloads/)

## Install python virtual env to run the code

Clone the github repo 
`git clone https://github.com/pratiknag/Bivariate_DeepKriging.git`

cd to the directory `reproducible/`

Check if `pip` is installed

`$ pip --version`

If `pip` is not installed, follow steps below:

```
$ cd ~
$ curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
$ python3 get-pip.py
```

Install virtual environment first & then activate:

```
$ python3 -m pip install --user virtualenv #Install virtualenv if not installed in your system
$ python3 -m virtualenv env #Create virtualenv for your project
$ source env/bin/activate #Activate virtualenv for linux/MacOS
```

Install all dependencies for your project from `requirements.txt` file:

```
$ pip install -r requirements.txt
```

## Reproducing results

To generate synthetic non-Gaussian and non-stationary datasets run

```
num_sim=10
Rscript src/R_scripts/Data_generation.R $PWD $num_sim
```
Note here num_sim is number of simulations. 

The results related to Bivariate deepkriging can be produced by running the following python scripts 

```
num_sim=10
python3 src/python_scripts/2d_non-Gaussian-cov.py $num_sim
python3 src/python_scripts/2d_nonstationary.py $num_sim
```
The above scripts will generate two .csv files inside the plot_results/ directory. These .csv files contain the mean square error for each simulation modeled with bivariate DeepKriging.

The results related to Gaussian Kriging can be produced by running the following R scripts. Note that these scripts require running the code in multiple cores in parallel and hence can be computationally intensive. Please note that in the following code we have introduced a binary parameter called subset which takes values 0 or 1. If set to 1 the training of the Gaussian Kriging will happen on a smaller subset which is randomly selected from the whole training dataset. To run the 100 simulations on a linux based system with 56 cores took us 4+ hours for the training subset. If set to 0 the training will happen on the whole training dataset.   

```
num_sim=10
subset=1
Rscript src/R_scripts/Gaussian_kriging_with_covariates.R $PWD $num_sim $subset
Rscript src/R_scripts/Gaussian_kriging_nonstat.R $PWD $num_sim $subset
Rscript src/R_scripts/calculate_validation_results_GKrigning.R $PWD $num_sim
``` 

In the final step the plots can be generated by running the following R script 
```
Rscript src/R_scripts/plot_results_reproducible.R $PWD
```

