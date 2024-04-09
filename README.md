# Bivariate DeepKriging
Deep learning based spatial interpolation for bivariate spatial processes.

## Pre-requisites

Please ensure that you have R installed on your system, along with the following libraries:
```
geoR
MASS
fields
parallel
doParallel
```
Additionally, please verify if Python 3+ is installed. If not please download and install python from [here](https://www.python.org/downloads/)

## Install python virtual env to run the code

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
$ git clone git@github.com:pratik187/Bivariate_DeepKriging.git
$ cd Bivariate_DeepKriging
$ python3 -m pip install --user virtualenv #Install virtualenv if not installed in your system
$ python3 -m virtualenv env #Create virtualenv for your project
$ source env/bin/activate #Activate virtualenv for linux/MacOS
```

Install all dependencies for your project from `requirements.txt` file:

```
$ pip install -r requirements.txt
```

## Reproducing results

Results can be reproduced by running the following command:
```
bash run.sh
```
Running the above will create files related to validation results for the non-Gaussian and non-stationary simulation scenarios in the `plot_results` directory. 
To create the boxplots as shown in the paper run the file `R_scripts/plot_results_reproducible.R` in R Studio.
