rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
file_path = args[1]
setwd(file_path)

num_sim = 100

mainDir <- "."
subDir <- "plot_results/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

# MSE for covariates LMC
a = readRDS("GKriging_cov_results/pred_list_LMC_cov.rds")
mse1 = rep(NA,num_sim)
mse2 = rep(NA,num_sim)
mpiw1 = rep(NA,num_sim)
mpiw2 = rep(NA,num_sim)
picp1 = rep(NA, num_sim)
picp2 = rep(NA,num_sim)
for(i in 1:num_sim){
  file_path = paste0("synthetic_data_simulations_nonGaussian_cov/testing_data/2D_nonGaussian_1200_",i,"-test.csv")
  test_data = read.csv(file_path,  header = T)
  mse1[i] = mean((a[[i]]$pred[1:120] - test_data$var1)^2)
  mse2[i] = mean((a[[i]]$pred[121:240] - test_data$var2)^2)
  ub = (a[[i]]$pred[1:120]+1.96*sqrt(diag(a[[i]]$conditional_var)[1:120]))
  lb = (a[[i]]$pred[1:120]-1.96*sqrt(diag(a[[i]]$conditional_var)[1:120]))
  width = ub - lb
  count_var = 0
  for(j in 1:dim(test_data)[1]){
    if((test_data$var1[j] > lb[j]) & (test_data$var1[j] < ub[j])) count_var = count_var + 1
  }
  picp1[i] = count_var/120
  mpiw1[i] = mean(width)
  ub = (a[[i]]$pred[121:240]+1.96*sqrt(diag(a[[i]]$conditional_var)[121:240]))
  lb = (a[[i]]$pred[121:240]-1.96*sqrt(diag(a[[i]]$conditional_var)[121:240]))
  count_var = 0
  for(j in 1:dim(test_data)[1]){
    if((test_data$var1[j] > lb[j]) & (test_data$var1[j] < ub[j])) count_var = count_var + 1
  }
  width = ub - lb
  mpiw2[i] = mean(width)
  picp2[i] = count_var/120
}
print(paste0("mean picp Variable 1 ",mean(picp1),", mean mpiw Variable 1 ",mean(mpiw1)))
print(paste0("mean picp Variable 2 ",mean(picp2),", mean mpiw Variable 2 ",mean(mpiw2)))
mse_df = data.frame(mse_var1 = mse1, picp1 = picp1, mpiw1 = mpiw1, 
                    mse_var2 = mse2, picp2 = picp2, mpiw2 = mpiw2)
write.csv(mse_df, "plot_results/validation_LMC_cov.csv",row.names = F)

# MSE for covariates parsimonious Gaussian
a = readRDS("GKriging_cov_results/pred_list_par_Matern_cov.rds")
mse1 = rep(NA,num_sim)
mse2 = rep(NA,num_sim)
mpiw1 = rep(NA,num_sim)
mpiw2 = rep(NA,num_sim)
picp1 = rep(NA, num_sim)
picp2 = rep(NA,num_sim)
for(i in 1:num_sim){
  file_path = paste0("synthetic_data_simulations_nonGaussian_cov/testing_data/2D_nonGaussian_1200_",i,"-test.csv")
  test_data = read.csv(file_path,  header = T)
  mse1[i] = mean((a[[i]]$pred[1:120] - test_data$var1)^2)
  mse2[i] = mean((a[[i]]$pred[121:240] - test_data$var2)^2)
  ub = (a[[i]]$pred[1:120]+1.96*sqrt(diag(a[[i]]$conditional_var)[1:120]))
  lb = (a[[i]]$pred[1:120]-1.96*sqrt(diag(a[[i]]$conditional_var)[1:120]))
  width = ub - lb
  count_var = 0
  for(j in 1:dim(test_data)[1]){
    if((test_data$var1[j] > lb[j]) & (test_data$var1[j] < ub[j])) count_var = count_var + 1
  }
  picp1[i] = count_var/120
  mpiw1[i] = mean(width)
  ub = (a[[i]]$pred[121:240]+1.96*sqrt(diag(a[[i]]$conditional_var)[121:240]))
  lb = (a[[i]]$pred[121:240]-1.96*sqrt(diag(a[[i]]$conditional_var)[121:240]))
  count_var = 0
  for(j in 1:dim(test_data)[1]){
    if((test_data$var1[j] > lb[j]) & (test_data$var1[j] < ub[j])) count_var = count_var + 1
  }
  width = ub - lb
  mpiw2[i] = mean(width)
  picp2[i] = count_var/120
}
print(paste0("mean picp Variable 1 ",mean(picp1),", mean mpiw Variable 1 ",mean(mpiw1)))
print(paste0("mean picp Variable 2 ",mean(picp2),", mean mpiw Variable 2 ",mean(mpiw2)))

mse_df = data.frame(mse_var1 = mse1, picp1 = picp1, mpiw1 = mpiw1, 
                    mse_var2 = mse2, picp2 = picp2, mpiw2 = mpiw2)
write.csv(mse_df, "plot_results/validation_par_Matern_cov.csv",row.names = F)








# MSE for nonstationary LMC
a = readRDS("GKriging_nonstat_results/pred_list_nonstat_lmc.rds")
mse1 = rep(NA,num_sim)
mse2 = rep(NA,num_sim)
mpiw1 = rep(NA,num_sim)
mpiw2 = rep(NA,num_sim)
picp1 = rep(NA, num_sim)
picp2 = rep(NA,num_sim)
for(i in 1:num_sim){
  file_path = paste0("synthetic_data_simulations_non-Stationary/testing_data/2D_nonstationary_1200_",i,"-test.csv")
  test_data = read.csv(file_path,  header = T)
  mse1[i] = mean((a[[i]]$pred[1:120] - test_data$var1)^2)
  mse2[i] = mean((a[[i]]$pred[121:240] - test_data$var2)^2)
  ub = (a[[i]]$pred[1:120]+1.96*sqrt(diag(a[[i]]$conditional_var)[1:120]))
  lb = (a[[i]]$pred[1:120]-1.96*sqrt(diag(a[[i]]$conditional_var)[1:120]))
  width = ub - lb
  count_var = 0
  for(j in 1:dim(test_data)[1]){
    if((test_data$var1[j] > lb[j]) & (test_data$var1[j] < ub[j])) count_var = count_var + 1
  }
  picp1[i] = count_var/120
  mpiw1[i] = mean(width)
  ub = (a[[i]]$pred[121:240]+1.96*sqrt(diag(a[[i]]$conditional_var)[121:240]))
  lb = (a[[i]]$pred[121:240]-1.96*sqrt(diag(a[[i]]$conditional_var)[121:240]))
  count_var = 0
  for(j in 1:dim(test_data)[1]){
    if((test_data$var1[j] > lb[j]) & (test_data$var1[j] < ub[j])) count_var = count_var + 1
  }
  width = ub - lb
  mpiw2[i] = mean(width)
  picp2[i] = count_var/120
}
print(paste0("mean picp Variable 1 ",mean(picp1),", mean mpiw Variable 1 ",mean(mpiw1)))
print(paste0("mean picp Variable 2 ",mean(picp2),", mean mpiw Variable 2 ",mean(mpiw2)))
mse_df = data.frame(mse_var1 = mse1, picp1 = picp1, mpiw1 = mpiw1, 
                    mse_var2 = mse2, picp2 = picp2, mpiw2 = mpiw2)
write.csv(mse_df, "plot_results/validation_LMC_nonstat.csv",row.names = F)

# MSE for nonstat parsimonious Gaussian
a = readRDS("GKriging_nonstat_results/pred_list_nonstat_par_Matern.rds")
mse1 = rep(NA,num_sim)
mse2 = rep(NA,num_sim)
mpiw1 = rep(NA,num_sim)
mpiw2 = rep(NA,num_sim)
picp1 = rep(NA, num_sim)
picp2 = rep(NA,num_sim)
for(i in 1:num_sim){
  file_path = paste0("synthetic_data_simulations_non-Stationary/testing_data/2D_nonstationary_1200_",i,"-test.csv")
  test_data = read.csv(file_path,  header = T)
  mse1[i] = mean((a[[i]]$pred[1:120] - test_data$var1)^2)
  mse2[i] = mean((a[[i]]$pred[121:240] - test_data$var2)^2)
  ub = (a[[i]]$pred[1:120]+1.96*sqrt(diag(a[[i]]$conditional_var)[1:120]))
  lb = (a[[i]]$pred[1:120]-1.96*sqrt(diag(a[[i]]$conditional_var)[1:120]))
  width = ub - lb
  count_var = 0
  for(j in 1:dim(test_data)[1]){
    if((test_data$var1[j] > lb[j]) & (test_data$var1[j] < ub[j])) count_var = count_var + 1
  }
  picp1[i] = count_var/120
  mpiw1[i] = mean(width)
  ub = (a[[i]]$pred[121:240]+1.96*sqrt(diag(a[[i]]$conditional_var)[121:240]))
  lb = (a[[i]]$pred[121:240]-1.96*sqrt(diag(a[[i]]$conditional_var)[121:240]))
  count_var = 0
  for(j in 1:dim(test_data)[1]){
    if((test_data$var1[j] > lb[j]) & (test_data$var1[j] < ub[j])) count_var = count_var + 1
  }
  width = ub - lb
  mpiw2[i] = mean(width)
  picp2[i] = count_var/120
}
print(paste0("mean picp Variable 1 ",mean(picp1),", mean mpiw Variable 1 ",mean(mpiw1)))
print(paste0("mean picp Variable 2 ",mean(picp2),", mean mpiw Variable 2 ",mean(mpiw2)))

mse_df = data.frame(mse_var1 = mse1, picp1 = picp1, mpiw1 = mpiw1, 
                    mse_var2 = mse2, picp2 = picp2, mpiw2 = mpiw2)
write.csv(mse_df, "plot_results/validation_par_Matern_nonstat.csv",row.names = F)




