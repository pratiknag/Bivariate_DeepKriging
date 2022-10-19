library(geoR)
library(MASS)
library(fields)
setwd("/home/nagp/Desktop/Ghulam_work/final_folder/Bivariate_DeepKriging/")
mainDir <- "."
subDir <- "synthetic_data_simulations/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "./real_data/"
subDir <- "estimation/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

mainDir <- "./real_data/"
subDir <- "splitted_data/"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
########### Biv nongaussian ##################

# set.seed(12345567)
x <- seq(0,1, length.out = 80)
y <- seq(0,1, length.out = 80)

d1 <- expand.grid(x = x, y = y)
X = d1$x              # X, Y co-ordinates getting generated here
Y = d1$y
m <- as.matrix(dist(data.frame(X=X,Y=Y)))

R = 0.8
del1 = 0.7
del2 = 0.9
R_corr1 = 0.3
R_corr2 = 0.88

s11 = 0.89
s22 = 1.3

nu11 = 0.8 #0.2 ## for bivariate stationary the neu value was 0.4 0.6
nu22 = 0.8 #0.7
nu12 = (nu11 + nu22) /2 #+ (del1*(1-R_corr1))

alpha11 = 0.2
alpha22 = 0.4
alpha12 = (alpha11 + alpha22)/2 #+ (del2*(1-R_corr2))


s12 = sqrt(s11*s22)*(alpha11^(nu11/2))*(alpha22^(nu22/2))*gamma(nu12)*R / 
  ((alpha12^nu12)*sqrt(gamma(nu11)*gamma(nu22)))

constant = (sqrt(s11*s22)*R*(gamma(nu12)/sqrt(gamma(nu11)*gamma(nu22)))) /
  ((alpha12^nu12)/sqrt(alpha11^nu11 * alpha22^nu22))

matern_cov1 = s11*matern(m,sqrt(alpha11),nu11)
matern_cov2 <- constant*matern(m,sqrt(alpha12),nu12)
matern_cov4 = s22*matern(m,sqrt(alpha22),nu22)
full_matern_cov = rbind(cbind(matern_cov1,matern_cov2),cbind(t(matern_cov2),matern_cov4))

for(i in 1:50){
  simulation = mvrnorm(1,rep(0,2*1200),full_matern_cov)
  
  var1 = simulation[1:1200]
  var2 = simulation[1201:(2*1200)]
  g = 0.5
  h = 1.5
  
  tukey_var1 = ((exp(g*var1)-1)/g) * exp(h*(var1^2)/2)
  
  g = -0.4
  h = 1.3
  tukey_var2 = ((exp(g*var2)-1)/g) * exp(h*(var2^2)/2)
  
  par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))
  
  quilt.plot(X,Y, tukey_var1, main =
               "variable 1", nx = 40, ny = 30)
  quilt.plot(X,Y, tukey_var2,main = "variable 2", nx = 40, ny = 30)
  
  df = data.frame(x=X,y=Y,var1 = tukey_var1, var2 = tukey_var2)
  write.csv(x = df,file = paste0("synthetic_data_simulations/2d_nongaussian_1200_",i,".csv"))
}


## Single data generation 6400
simulation = mvrnorm(1,rep(0,2*6400),full_matern_cov)
var1 = simulation[1:6400]
var2 = simulation[6401:(2*6400)]


# 
# 
# 
df = data.frame(x=X,y=Y,var1 = var1, var2 = var2)

quilt.plot(X,Y, var1, main =
             "variable 1", nx = 80, ny = 80)
quilt.plot(X,Y, var2,main = "variable 2", nx = 80, ny = 80)
  write.csv(x = df,file = "synthetic_datasets/2D_biv_matern_6400.csv")
  # write.csv(x = df,file = "2D_biv_matern_large.csv")
  
################# gaussian ####################  
setwd("/home/nagp/Desktop/Ghulam_work/")
# set.seed(12345567)
x <- seq(0,1, length.out = 40)
y <- seq(0,1, length.out = 30)

d1 <- expand.grid(x = x, y = y)
X = d1$x              # X, Y co-ordinates getting generated here
Y = d1$y
m <- as.matrix(dist(data.frame(X=X,Y=Y)))

R = 0.8
del1 = 0.7
del2 = 0.9
R_corr1 = 0.3
R_corr2 = 0.88

s11 = 0.89
s22 = 1.3

nu11 = 0.8 #0.2 ## for bivariate stationary the neu value was 0.4 0.6
nu22 = 0.8 #0.7
nu12 = (nu11 + nu22) /2 #+ (del1*(1-R_corr1))

alpha11 = 0.2
alpha22 = 0.4
alpha12 = (alpha11 + alpha22)/2 #+ (del2*(1-R_corr2))


s12 = sqrt(s11*s22)*(alpha11^(nu11/2))*(alpha22^(nu22/2))*gamma(nu12)*R / 
  ((alpha12^nu12)*sqrt(gamma(nu11)*gamma(nu22)))

constant = (sqrt(s11*s22)*R*(gamma(nu12)/sqrt(gamma(nu11)*gamma(nu22)))) /
  ((alpha12^nu12)/sqrt(alpha11^nu11 * alpha22^nu22))

matern_cov1 = s11*matern(m,sqrt(alpha11),nu11)
matern_cov2 <- constant*matern(m,sqrt(alpha12),nu12)
matern_cov4 = s22*matern(m,sqrt(alpha22),nu22)
full_matern_cov = rbind(cbind(matern_cov1,matern_cov2),cbind(t(matern_cov2),matern_cov4))

for(i in 1:50){
  simulation = mvrnorm(1,rep(0,2*1200),full_matern_cov)
  
  var1 = simulation[1:1200]
  var2 = simulation[1201:(2*1200)]
  # g = 0.5
  # h = 1.5
  # 
  # tukey_var1 = ((exp(g*var1)-1)/g) * exp(h*(var1^2)/2)
  # 
  # g = -0.4
  # h = 1.3
  # tukey_var2 = ((exp(g*var2)-1)/g) * exp(h*(var2^2)/2)
  
  par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))
  
  quilt.plot(X,Y, var1, main =
               "variable 1", nx = 40, ny = 30)
  quilt.plot(X,Y, var2,main = "variable 2", nx = 40, ny = 30)
  
  df = data.frame(x=X,y=Y,var1 = tukey_var1, var2 = tukey_var2)
  write.csv(x = df,file = paste0("synthetic_data_simulations/2d_gaussian_1200_",i,".csv"))
}




################# Nonstationary noise ###############
  library(geoR)
  library(MASS)
  library(fields)
  #RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
  ##                   RFoptions(seed=NA) to make them all random again
  # x <- seq(0,1, length.out = 1000
  
  x <- seq(0,1, length.out = 40)
  y <- seq(0,1, length.out = 30)
  
  d1 <- expand.grid(x = x, y = y)
  X = d1$x              # X, Y co-ordinates getting generated here
  Y = d1$y
  m <- as.matrix(dist(data.frame(X=X,Y=Y)))
  
  R = 0.2
  del1 = 0.7
  del2 = 0.9
  R_corr1 = 0.3
  R_corr2 = 0.88
  
  s11 = 0.001
  s22 = 0.001
  
  nu11 = 0.5 #0.2 ## for bivariate stationary the neu value was 0.4 0.6
  nu22 = 0.8 #0.7
  nu12 = (nu11 + nu22) /2 #+ (del1*(1-R_corr1))
  
  alpha11 = 0.55
  alpha22 = 0.98
  alpha12 = (alpha11 + alpha22)/2 #+ (del2*(1-R_corr2))
  
  
  s12 = sqrt(s11*s22)*(alpha11^(nu11/2))*(alpha22^(nu22/2))*gamma(nu12)*R / 
    ((alpha12^nu12)*sqrt(gamma(nu11)*gamma(nu22)))
  
  constant = (sqrt(s11*s22)*R*(gamma(nu12)/sqrt(gamma(nu11)*gamma(nu22)))) /
    ((alpha12^nu12)/sqrt(alpha11^nu11 * alpha22^nu22))
  
  matern_cov1 = s11*matern(m,sqrt(alpha11),nu11)
  matern_cov2 <- constant*matern(m,sqrt(alpha12),nu12)
  matern_cov4 = s22*matern(m,sqrt(alpha22),nu22)
  full_matern_cov = rbind(cbind(matern_cov1,matern_cov2),cbind(t(matern_cov2),matern_cov4))
  simulation = mvrnorm(1,rep(0,2*1200),full_matern_cov)
  
  var1 = simulation[1:1200]
  var2 = simulation[1201:2400]
  
  par(mfrow=c(1,2))
  quilt.plot(X,Y, var1, main = 
               "variable 1", nx = 40, ny = 30)
  quilt.plot(X,Y, var2,main = "variable 2", nx = 40, ny = 30)
  
  df = data.frame(x=X,y=Y,var1 = var1, var2 = var2)
  write.csv(x = df,file = "synthetic_data_simulations/2D_matern_nonstatiopnary_noise.csv")
  
