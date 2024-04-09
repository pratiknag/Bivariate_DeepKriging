rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
file_path = args[1]
setwd(file_path)
library(parallel)
library(doParallel)
library(geoR)
library(MASS)
library(fields)
source("R_scripts/GKriging_functions.R")
num_sim = args[2]
file_path = paste0("synthetic_data_simulations_nonGaussian_cov/training_data/2D_nonGaussian_1200_",1,"-train.csv")
  data_train = read.csv(file_path, header = T)
  file_path = paste0("synthetic_data_simulations_nonGaussian_cov/testing_data/2D_nonGaussian_1200_",1,"-test.csv")
  data_test = read.csv(file_path, header = T)
  init.ind = c(1,1,1,1,1,1,0,0)
  indmat.estim = optim_indmat_loglik(init.ind)
cl <- makeCluster(parallel::detectCores())
registerDoParallel(cl)   #create a cluster
parallel.function.lmc = function(sim){
	source("R_scripts/GKriging_functions.R")
  file_path = paste0("synthetic_data_simulations_nonGaussian_cov/training_data/2D_nonGaussian_1200_",sim,"-train.csv")
  data_train = read.csv(file_path, header = T)
  file_path = paste0("synthetic_data_simulations_nonGaussian_cov/testing_data/2D_nonGaussian_1200_",sim,"-test.csv")
  data_test = read.csv(file_path, header = T)

  init.lmc = c(indmat.estim$par[c(1:6)],0,0,0,0,indmat.estim$par[7:8])
  lmc.loglikelihood=function(par)
{

  return(lmc.loglikelihood_allcomponents(p=par,data_train)$mlv)
}

optim_lmc.loglikelihood = function(par){
  optim(par=par,
        fn = lmc.loglikelihood,
        hessian=FALSE,
        control=list(trace=4,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=300))
}
  fit.Model.lmc = optim_lmc.loglikelihood(par=init.lmc)
  pred = lmc.pred.summary(fit.Model.lmc$par, data_train, data_test)
  # mse1 = mean((pred$pred[1:120] - data_test$var1)^2)
  # mse2 = mean((pred$pred[121:240] - data_test$var2)^2)
  return(pred)
}
parallel.function.par_Matern = function(sim){
  source("R_scripts/GKriging_functions.R")
  file_path = paste0("synthetic_data_simulations_nonGaussian_cov/training_data/2D_nonGaussian_1200_",sim,"-train.csv")
  data_train = read.csv(file_path, header = T)
  file_path = paste0("synthetic_data_simulations_nonGaussian_cov/testing_data/2D_nonGaussian_1200_",sim,"-test.csv")
  data_test = read.csv(file_path, header = T)
  
  init.par_Matern = c(indmat.estim$par[c(1:8)],0)
  par_Matern.loglikelihood=function(par)
  {
    
    return(par_Matern.loglikelihood_allcomponents(p=par,data_train)$mlv)
  }
  
  optim_par_Matern.loglikelihood = function(par){
    optim(par=par,
          fn = par_Matern.loglikelihood,
          hessian=FALSE,
          control=list(trace=4,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=300))
  }
  fit.Model.par_Matern = optim_par_Matern.loglikelihood(par=init.par_Matern)
  pred = par_Matern.pred.summary(fit.Model.par_Matern$par, data_train, data_test)
  # mse1 = mean((pred$pred[1:120] - data_test$var1)^2)
  # mse2 = mean((pred$pred[121:240] - data_test$var2)^2)
  return(pred)
}


print(getwd())
mainDir <- "."
subDir <- "GKriging_cov_results/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
fx1 = foreach(sim=1:num_sim) %dopar% parallel.function.lmc(sim)
saveRDS(fx1, "GKriging_cov_results/pred_list_LMC_cov.rds")
fx2 = foreach(sim=1:num_sim) %dopar% parallel.function.par_Matern(sim)
saveRDS(fx2, "GKriging_cov_results/pred_list_par_Matern_cov.rds")


