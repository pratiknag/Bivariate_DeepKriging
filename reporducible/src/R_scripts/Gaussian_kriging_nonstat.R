rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
file_path = args[1]
setwd(file_path)
library(parallel)
library(doSNOW)
library(geoR)
library(MASS)
library(fields)
print("[*][*][*][*][*] Fitting co-Kriging for Non-Stationary Data [*][*][*][*][*]")
print("")
mainDir <- "."
subDir <- "GKriging_nonstat_results/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
num_sim = args[2]
subset = args[3]
source("src/R_scripts/GKriging_functions.R")

file_path = paste0("synthetic_data_simulations_non-Stationary/training_data/2D_nonstationary_1200_",1,"-train.csv")
data_train = read.csv(file_path, header = T)
if(subset == 1){
  sample1 = sample(1:1080, 200)
  data_train = data_train[sample1,]
}
file_path = paste0("synthetic_data_simulations_non-Stationary/testing_data/2D_nonstationary_1200_",1,"-test.csv")
data_test = read.csv(file_path, header = T)
init.ind = c(1,1,1,1,1,1,0,0)
 
print("[*][*][*][*][*] Fitting Independent likelihood to get initial parameter estimates [*][*][*][*][*]")
indmat.estim = optim_indmat_loglik(init.ind)
print("[*][*][*][*][*] Fitting Independent likelihood ---- Process finised [*][*][*][*][*]")
  
  
parallel.function.lmc = function(sim){
	source("src/R_scripts/GKriging_functions.R")
  file_path = paste0("synthetic_data_simulations_non-Stationary/training_data/2D_nonstationary_1200_",sim,"-train.csv")
  data_train = read.csv(file_path, header = T)
  if(subset == 1){
    sample1 = sample(1:1080, 200)
    data_train = data_train[sample1,]
  }
  file_path = paste0("synthetic_data_simulations_non-Stationary/testing_data/2D_nonstationary_1200_",sim,"-test.csv")
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
  
  return(pred)
}
parallel.function.par_Matern = function(sim){
  source("src/R_scripts/GKriging_functions.R")
  file_path = paste0("synthetic_data_simulations_non-Stationary/training_data/2D_nonstationary_1200_",sim,"-train.csv")
  data_train = read.csv(file_path, header = T)
  if(subset == 1){
    sample1 = sample(1:1080, 200)
    data_train = data_train[sample1,]
  }
  file_path = paste0("synthetic_data_simulations_non-Stationary/testing_data/2D_nonstationary_1200_",sim,"-test.csv")
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
  # print(paste0("mse for variable 1 is ",mse1, " and variable 2 is ", mse2))
  return(pred)
}


print("[*][*][*][*][*] Fitting cokriging LMC [*][*][*][*][*]")
print("")
a = c(1:num_sim)
pb = txtProgressBar(max=length(a), style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
cl = makeCluster(parallel::detectCores() - 2)
registerDoSNOW(cl) #create a cluster

fx1 = foreach(sim=1:num_sim, .options.snow=opts) %dopar% 
  {
    s = parallel.function.lmc(sim)
    # setTxtProgressBar(pb, sim) 
    return(s)
  }
close(pb)
stopCluster(cl) 
saveRDS(fx1, "GKriging_nonstat_results/pred_list_nonstat_lmc.rds")

print("[*][*][*][*][*] Fitting cokriging parsimonious Matern [*][*][*][*][*]")
print("")
a = c(1:num_sim)
pb = txtProgressBar(max=length(a), style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
cl = makeCluster(parallel::detectCores() - 2)
registerDoSNOW(cl) #create a cluster
# print("########## Inside foreach function ############")
fx2 = foreach(sim=1:num_sim, .options.snow=opts) %dopar% 
  { 
    s = parallel.function.par_Matern(sim)
    # setTxtProgressBar(pb, sim) 
    return(s)
  }
close(pb)
stopCluster(cl) 
saveRDS(fx2, "GKriging_nonstat_results/pred_list_nonstat_par_Matern.rds")

