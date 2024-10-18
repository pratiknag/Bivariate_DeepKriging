rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
file_path = args[1] #"/Users/nagp/Desktop/Biv.DeepKriging/Bivariate_DeepKriging/"
setwd(file_path)

library(geoR)
library(MASS)
library(fields)


# mainDir <- "."
# subDir <- "synthetic_data_simulations_nonGaussian/"
# dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
mainDir <- "."
subDir <- "Model_Example/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "."
subDir <- "plot_results/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "."
subDir <- "synthetic_data_simulations_nonGaussian_cov/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "."
subDir <- "synthetic_data_simulations_non-Stationary/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "synthetic_data_simulations_non-Stationary/"
subDir <- "training_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "synthetic_data_simulations_non-Stationary/"
subDir <- "testing_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "synthetic_data_simulations_nonGaussian_cov/"
subDir <- "training_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "synthetic_data_simulations_nonGaussian_cov/"
subDir <- "testing_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

# mainDir <- "synthetic_data_simulations_nonGaussian/"
# subDir <- "training_data/"
# dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
# 
# mainDir <- "synthetic_data_simulations_nonGaussian/"
# subDir <- "testing_data/"
# dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

num_sim = args[2]


print("####################### Simulating non-stationary datasets ########################")
##############################################
############# Non-Stationary #################
##############################################

set.seed(18)

x = seq(0,1, length.out = 80)
y = seq(0,1, length.out = 80)

d1 <- expand.grid(x = x, y = y)

sample1 = sample(1:6400,1200)
d1 = d1[sample1,]
s = data.matrix(d1)
N = dim(s)[1]
K = 0
num_basis = c(2^2,3^2,5^2)
knots_1d = vector("list",length(num_basis))
for(i in 1:length(num_basis)) knots_1d[[i]] = seq(0,1,length.out = sqrt(num_basis[i]))

phi = matrix(0, nrow = N, ncol = sum(num_basis))

for(res in 1:length(num_basis)){
  theta = 1
  d = data.matrix(expand.grid(x = knots_1d[[res]],y = knots_1d[[res]]))
  for(i in 1:num_basis[res]){
    d_norm = sqrt(rowSums((s-d[i,])^2))
    for(j in 1:length(d_norm)){
      if(d_norm[j] >=0 & d_norm[j] <= 1) phi[j,i + K] = (1-d_norm[j])^6 * (35 * d_norm[j]^2 + 18 * d_norm[j] + 3)/3
      else phi[j,i + K] = 0
    }
  }
  K = K + num_basis[res]
}

write.csv(phi,"src/python_scripts/phi.csv",  row.names = F)

for(i in 1:num_sim){
  even = seq(2,sum(num_basis),2)
  even_cols = phi[,even]
  odd_cols = phi[,-even]
  a = sample(seq(-2.5,2.5,length.out = 100),1)
  b = sample(seq(-2.5,2.5,length.out = 100),1)
  c = sample(seq(-2.5,2.5,length.out = 100),1)
  d = sample(seq(-2.5,2.5,length.out = 100),1)
  e = sample(seq(-2.5,2.5,length.out = 100),1)
  f = sample(seq(-2.5,2.5,length.out = 100),1)
  
  var1 = rowSums((a*even_cols^(3/2) + c*odd_cols) - b*sqrt(even_cols*odd_cols))
  var2 = rowSums((d*even_cols - f*odd_cols^(3/2)))
  # print(cor(var1,var2))
  # par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))
  # 
  # fields::quilt.plot(s[,1],s[,2], var1, main =
  #              "variable 1", nx = 30, ny = 30)
  # fields::quilt.plot(s[,1],s[,2], var2,main = "variable 2", nx = 30, ny = 30)
  
  df = data.frame(x=s[,1],y=s[,2],var1 = var1, var2 = var2)
  write.csv(x = df,file = paste0("synthetic_data_simulations_non-Stationary/2d_nonstationary_1200_",i,".csv"), 
            row.names = F)
  sample1 = sample(1:1200, 1080)
  df_train = df[sample1,]
  df_test = df[-sample1,]
  write.csv(x = df_train,file = paste0("synthetic_data_simulations_non-Stationary/training_data/2D_nonstationary_1200_",i,"-train.csv"), 
            row.names = F)
  write.csv(x = df_test,file = paste0("synthetic_data_simulations_non-Stationary/testing_data/2D_nonstationary_1200_",i,"-test.csv"), 
            row.names = F)

}

print("####################### Simulating non-Gaussian datasets ########################")
##############################################
######## NonGaussian with covariates #########
##############################################
set.seed(12345567)

x = seq(0,1, length.out = 80)
y = seq(0,1, length.out = 80)

d1 <- expand.grid(x = x, y = y)

sample1 = sample(1:6400,1200)
d1 = d1[sample1,]
X = d1$x              # X, Y co-ordinates getting generated here
Y = d1$y
m <- as.matrix(dist(data.frame(X=X,Y=Y)))

## the nonGaussian field
R = 0.5
del1 = 0.7
del2 = 0.9
R_corr1 = 0.3
R_corr2 = 0.88

s11 = 0.7
s22 = 0.8

nu11 = 0.3 #0.2 ## for bivariate stationary the neu value was 0.4 0.6
nu22 = 0.6 #0.7
nu12 = (nu11 + nu22) /2 #+ (del1*(1-R_corr1))

alpha11 = 0.05
alpha22 = 0.1
alpha12 = (alpha11 + alpha22)/2 #+ (del2*(1-R_corr2))


s12 = sqrt(s11*s22)*(alpha11^(nu11/2))*(alpha22^(nu22/2))*gamma(nu12)*R / 
  ((alpha12^nu12)*sqrt(gamma(nu11)*gamma(nu22)))

constant = (sqrt(s11*s22)*R*(gamma(nu12)/sqrt(gamma(nu11)*gamma(nu22)))) /
  ((alpha12^nu12)/sqrt(alpha11^nu11 * alpha22^nu22))

matern_cov1 = s11*matern(m,sqrt(alpha11),nu11)
matern_cov2 <- constant*matern(m,sqrt(alpha12),nu12)
matern_cov4 = s22*matern(m,sqrt(alpha22),nu22)
full_matern_cov = rbind(cbind(matern_cov1,matern_cov2),cbind(t(matern_cov2),matern_cov4))

## mean function with covariates 
a = 0.1
s = 0.9
nu = 0.5
cov = s*matern(m,a,nu)  

for(i in 1:num_sim){
  x1 = mvrnorm(1,rep(sample(seq(-2.5,2.5,length.out = 100),1),1200),cov)
  x2 = mvrnorm(1,rep(sample(seq(-2.5,2.5,length.out = 100),1),1200),cov)
  x3 = mvrnorm(1,rep(sample(seq(-2.5,2.5,length.out = 100),1),1200),cov)
  x4 = mvrnorm(1,rep(sample(seq(-2.5,2.5,length.out = 100),1),1200),cov)
  x5 = mvrnorm(1,rep(sample(seq(-2.5,2.5,length.out = 100),1),1200),cov)
  
  mean = x1^2 -x2^2 + x3^2 - x4^2 - x5^2 + 2*x1*x2 + 3*x2*x3 -2*x3*x5 + 10*x1*x4 + sin(x1)*x2*x3 + 
    cos(x2)*x3*x5 + (x1*x2*x4*x5)
  simulation = mvrnorm(1,rep(0,2*1200),full_matern_cov)
  
  var1 = simulation[1:1200]
  var2 = simulation[1201:(2*1200)]
  g = 0.8
  h = 0.5
  
  tukey_var1 = ((exp(g*var1)-1)/g) * exp(h*(var1^2)/2) + mean
  
  g = -0.8
  h = 0.5
  tukey_var2 = ((exp(g*var2)-1)/g) * exp(h*(var2^2)/2) + mean 
  
  # par(mfrow=c(1,3), mar=c(5.1,4.1,4.1,2.1))
  # 
  # quilt.plot(X,Y, mean, main =
  #              "variable 1", nx = 35, ny = 35)
  # quilt.plot(X,Y, tukey_var1, main =
  #              "variable 1", nx = 35, ny = 35)
  # quilt.plot(X,Y, tukey_var2,main = "variable 2", nx = 35, ny = 35)
  
  df = data.frame(cov1 = x1, cov2 = x2, 
                  cov3 = x3, cov4 = x4,
                  cov5 = x5, x=X,y=Y,
                  var1 = tukey_var1, var2 = tukey_var2)
  write.csv(x = df,file = paste0("synthetic_data_simulations_nonGaussian_cov/2d_nongaussian_1200_",i,".csv"), 
            row.names = F)
  sample1 = sample(1:1200, 1080)
  df_train = df[sample1,]
  df_test = df[-sample1,]
  write.csv(x = df_train,file = paste0("synthetic_data_simulations_nonGaussian_cov/training_data/2D_nonGaussian_1200_",i,"-train.csv"), 
            row.names = F)
  write.csv(x = df_test,file = paste0("synthetic_data_simulations_nonGaussian_cov/testing_data/2D_nonGaussian_1200_",i,"-test.csv"), 
            row.names = F)
}


