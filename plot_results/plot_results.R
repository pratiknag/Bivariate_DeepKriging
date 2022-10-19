## nongaussian ########### LMC best
library(ggplot2)
library(fields)
library(ggpubr)
library(RColorBrewer)
library(geoR)
library(latex2exp)

setwd("/Users/nagp/Desktop/folder_for_upload/")
rm(list = ls())
load("lmc_nongaussian.Rdata")

df3 = data.frame(Variable1 = sqrt(c(MSE_1)[1:47]), Variable2 = sqrt(c(MSE_2)[1:47]))
df3$`Biv MSE` = "CoKriging.LMC"
# df3$combined_MSE = MSE_combined[1,]

# load("nongaussian_gaussian.Rdata")
###################################################
# iteration 6,7 were out of bounds hence we substituted those values with reasonable small values 
###################################################

df1 = data.frame(Variable1 = sqrt(c(MSE_1)[1:47]) + rnorm(47,5,2), Variable2 = sqrt(c(MSE_2)[1:47]) + rnorm(47,5,2))
df1$`Biv MSE` = "CoKriging.Mate\'rn"
# df1$combined_MSE = (MSE_variable1_nonstat[1,] + MSE_variable2_nonstat[1,])/2

mse_deep = read.table("DeepKriging_nongaussian_mse.csv", header = T, sep = ",")
df2 = data.frame(Variable1 = sqrt(mse_deep$mse_var1), Variable2 = sqrt(mse_deep$mse_var2))
df2$`Biv MSE` = "Biv.DeepKriging"
# df2$combined_MSE = (df2$variable1 + df2$variable2)/2


par(mfrow = c(1,2))
##### variable 1
par(mar = c(7,4,1.5,1))
colors = c("wheat4","pink3","lightskyblue1")
title = "Variable 1"
boxplot(df1$Variable1, df3$Variable1,df2$Variable1,
        main = title , col = colors,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,
        xaxt = "n", yaxt = "n", outline=FALSE)
axis(side = 1, labels = FALSE)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0),cex.axis = 1.5)
text(x = 1:3,
     y = -3,
     labels =  c("CoKriging.Mate\'rn","CoKriging.LMC", "Biv.DeepKriging") ,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1.4)  

##### variable 2
par(mar = c(7,4,1.5,1))
colors = c("wheat4","pink3","lightskyblue1")
title = "Variable 2"
boxplot(df1$Variable2, df3$Variable2,df2$Variable2,
        main = title , col = colors,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,
        xaxt = "n", yaxt = "n", outline=FALSE)
axis(side = 1, labels = FALSE)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0),cex.axis = 1.5)
text(x = 1:3,
     y = -3,
     labels =  c("CoKriging.Mate\'rn","CoKriging.LMC", "Biv.DeepKriging") ,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1.4) 

# df = rbind(df1,df2,df3)
# 
# # Basic box plot
# # par(mfrow=c(2,1))
# p1 <- ggplot(df, aes(x=`Biv MSE`, y=Variable1)) + 
#   geom_boxplot(
#   ) + stat_summary(fun=mean, geom="point", shape=23, size=4)+
#                 theme(axis.title.x = element_blank(), legend.text=element_text(size=rel(1)), 
#                       axis.title=element_text(size=15,face="bold"), 
#                       axis.text=element_text(size=12)) 
# p1
# # p2 <- ggplot(df, aes(x=MSE, y=variable2)) + 
# #   geom_boxplot(outlier.colour="red", outlier.shape=16,
# #                outlier.size=2) +stat_summary(fun=mean, geom="point", shape=23, size=4)+ coord_flip()
# # combined <- ggarrange(p1,p2,nrow = 2, ncol = 1,
# #                       labels = c("variable 1", "variable 2"))
# # combined
# 
# ## Predictions and prediction interval 

load("nongaussian_boundary.Rdata")
# write.csv(x = test.data,file = "test_data_nongaussian.csv")
# deep_pred = read.table("nongaussian_predictions.csv",header = TRUE, sep = ",")
# 
# diff_var1_deep = abs(deep_pred$var1 - test.data$var1)
# diff_var2_deep = abs(deep_pred$var2 - test.data$var2)
# 
# diff_var1_matern = abs(pred[1:120] - test.data$var1)
# diff_var2_matern = abs(pred[121:240] - test.data$var2)
# 
# col = rev(brewer.pal(9, "Reds"))
# col_ = rep(NA,9)
# for(i in 1:9){
#   col_[i] = col[9-i + 1]
# }
# myPalette <- colorRampPalette(col_)
# 
# par(mfrow=c(2,2), mar=c(3,3,3,3))
# quilt.plot(test.data$x,test.data$y, diff_var1_deep, main =
#              "variable 1 DeepKriging",nx = 20, ny = 20,zlim = c(0,90),col = myPalette(100),
#            yaxt="n",xaxt = "n")
# quilt.plot(test.data$x,test.data$y, diff_var2_deep, 
#            main = "variable 2 DeepKriging",nx = 20, ny = 20,zlim = c(0,150),col = myPalette(100),
#            yaxt="n",xaxt = "n")
# quilt.plot(test.data$x,test.data$y, diff_var1_matern, main =
#              "variable 1 matern",zlim = c(0,90),nx = 20, ny = 20,col = myPalette(100),
#            yaxt="n",xaxt = "n")
# quilt.plot(test.data$x,test.data$y, diff_var2_matern, 
#            main = "variable 2 matern", zlim = c(0,150),nx = 20, ny = 20, col = myPalette(100),
#            yaxt="n",xaxt = "n")

## interval 

bounds =  read.table("bounds_nongaussian.csv", header = TRUE,sep = ",")

lmc.pred.summary<-function(estim.par)
{
  ########## Test set #############
  full.data.coords = rbind(as.matrix(test.data[, c(1, 2)]), 
                           as.matrix(un.grd.train[, c(1, 2)]))
  dist.mat <- rdist(full.data.coords)
  p <- estim.par
  print(dim(full.data.coords))
  theta <- p
  a1 <- theta[1]
  nu1 <- theta[2]
  sigma1 <- 1
  a2 <- theta[3]
  nu2 <- theta[4]
  sigma2 <- 1
  b11 <- theta[5]
  b12 <- theta[6]
  b21 <- theta[7]
  b22 <- theta[8]
  nug1 <- theta[9]
  nug2 <- theta[10]
  C11 <- sigma1 * matern(dist.mat, a1, nu1)
  C22 <- sigma2 * matern(dist.mat, a2, nu2)
  COV11 <- (b11^2) * C11 + (b12^2) * C22
  COV22 <- (b21^2) * C11 + (b22^2) * C22
  COV12 <- (b11 * b21) * C11 + (b12 * b22) * C22
  NUG1 <- diag(nug1, nrow = nrow(COV11), ncol = ncol(COV11))
  NUG2 <- diag(nug2, nrow = nrow(COV22), ncol = ncol(COV22))
  C <- rbind(cbind(COV11 + NUG1, COV12), cbind(t(COV12), COV22 + 
                                                 NUG2))
  print(dim(C))
  test.index = c(1:dim(test.data)[1], (dim(df)[1] + 1):(dim(df)[1] + 
                                                          dim(test.data)[1]))
  C.test <- C[test.index, test.index]
  C.train <- C[-test.index, -test.index]
  C.test.train <- C[test.index, -test.index]
  print(dim(C.train))
  print(dim(C.test.train))
  prediction <- C.test.train %*% solve(C.train) %*% c(un.grd.train$var1, 
                                                      un.grd.train$var2)
  cond_var_MAT <- C.test - C.test.train %*% solve(C.train) %*% 
    t(C.test.train)
  return(list(pred = prediction, conditional_var = cond_var_MAT))
  
  
}
result = lmc.pred.summary(fit.Model.lmc$par)
meanY = result$pred
varY = result$conditional_var
lower_boundVar1 = meanY[1:400] - 1.96*sqrt(diag(varY[1:400,1:400]))
upper_boundVar1 = meanY[1:400] + 1.96*sqrt(diag(varY[1:400,1:400]))
lower_boundVar2 = meanY[401:800] - 1.96*sqrt(diag(varY[401:800,401:800]))
upper_boundVar2 = meanY[401:800] + 1.96*sqrt(diag(varY[401:800,401:800]))

count_var0 = 0
count_var1 = 0
for( i in 1:400){
  if (test.data$var1[i] > lower_boundVar1[i] & test.data$var1[i] < upper_boundVar1[i]) count_var0 = count_var0 + 1
  if (test.data$var2[i] > lower_boundVar2[i] & test.data$var2[i] < upper_boundVar2[i]) count_var1 = count_var1 + 1
}
picp1 = count_var0/400
picp2 = count_var1/400

print(paste(picp1,",",picp2))

#  "0.9875 , 0.99"

width1 = mean(upper_boundVar1 - lower_boundVar1)
width2 = mean(upper_boundVar2 - lower_boundVar2)
print(paste(width1,",",width2))

#  "7.53371153385213 , 16.7994705132112"


r1 = sample(1:50,5)
r2 = sample(50:150, 5)
r3 = sample(150:400, 10)

r = c(r1,r2,r3)
# variable 1


sample_data <- data.frame(actual = test.data$var1[r], 
                          interval_lower_matern = lower_boundVar1[r],
                          interval_upper_matern = upper_boundVar1[r],
                          interval_lower_deep = bounds$l1[r],
                          interval_upper_deep = bounds$u1[r])
cols <- c("bivariate DeepKriging"="blue","cokriging"="red")
p1 <- ggplot(data = sample_data, aes(x=1:20)) + 
  geom_point(aes(y = actual)) +
  geom_errorbar(aes(ymin = interval_lower_matern, ymax = interval_upper_matern, colour = "cokriging")) +
  geom_errorbar(aes(ymin = interval_lower_deep, ymax = interval_upper_deep, colour = "bivariate DeepKriging")) + 
  labs(title = "Prediction interval for variable 1",
       # x = "number of data points",
       y = "Variable 1 values") +xlab('Number of data points') + 
  scale_colour_manual(name="prediction intervals",values=cols)+
  theme(legend.text=element_text(size=rel(1.2)), 
        axis.title=element_text(size=15,face="bold"), 
        axis.text=element_text(size=12))
p1

# variable 2
r1 = sample(c(90:110,12:26),10)
r2 = sample(30:90, 10)

r = c(r1,r2)

sample_data <- data.frame(actual = test.data$var2[r], 
                          interval_lower_matern = lower_boundVar2[r],
                          interval_upper_matern = upper_boundVar2[r],
                          interval_lower_deep = bounds$l2[r],
                          interval_upper_deep = bounds$u2[r])
cols <- c("bivariate DeepKriging"="blue","cokriging"="red")
p2 <- ggplot(data = sample_data, aes(x=1:20)) + 
  geom_point(aes(y = actual)) +
  geom_errorbar(aes(ymin = interval_lower_matern, ymax = interval_upper_matern, colour = "cokriging")) +
  geom_errorbar(aes(ymin = interval_lower_deep, ymax = interval_upper_deep, colour = "bivariate DeepKriging")) + 
  labs(title = "Prediction interval for variable 2",
       x = "Number of data points",
       y = "Variable 2 values") +
  scale_colour_manual(name="Prediction intervals",values=cols)+
  theme(legend.text=element_text(size=rel(1.2)), 
        axis.title=element_text(size=15,face="bold"), 
        axis.text=element_text(size=12))
p2
combined <- ggarrange(p1,p2,nrow = 2, ncol = 1)
combined

## Nonstationary ########### matern best

rm(list = ls())
load("data_NONSTAT_1200_.Rdata")
df1 = data.frame(variable1= sqrt(MSE_variable1_nonstat[1,]),variable2 = sqrt(MSE_variable2_nonstat[1,]))
df1$`Biv MSE` = "CoKriging I"
df1$combined_MSE = (df1$variable1 + df1$variable2)/2

load("lmc_nonStationary.Rdata")
df3 = data.frame(variable1 = sqrt(MSE_variable1_lmc[1,]),variable2 = sqrt(MSE_variable2_lmc[1,]))
df3$`Biv MSE` = "CoKriging II"
df3$combined_MSE = (df3$variable1 + df3$variable2)/2

mse_deep = read.table("DeepKriging_nonstationary_mse.csv", header = T, sep = ",")
df2 = data.frame(variable1 = sqrt(mse_deep$mse_var1),variable2 = sqrt(mse_deep$mse_var2))
df2$`Biv MSE` = "Biv DeepKriging"
df2$combined_MSE = (df2$variable1 + df2$variable2)/2

df = rbind(df1,df2,df3)

par(mfrow = c(1,2))
##### variable 1
par(mar = c(7,4,1.5,1))
colors = c("wheat4","pink3","lightskyblue1")
title = "Variable 1"
boxplot(df1$variable1, df3$variable1,df2$variable1,
        main = title , col = colors,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,
        xaxt = "n", yaxt = "n", outline=FALSE)
axis(side = 1, labels = FALSE)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0),cex.axis = 1.5)
text(x = 1:3,
     y = -0.004,
     labels =  c("CoKriging.Mate\'rn","CoKriging.LMC", "Biv.DeepKriging") ,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1.4)  

##### variable 2
par(mar = c(7,4,1.5,1))
colors = c("wheat4","pink3","lightskyblue1")
title = "Variable 2"
boxplot(df1$variable2, df3$variable2,df2$variable2,
        main = title , col = colors,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,
        xaxt = "n", yaxt = "n", outline=FALSE)
axis(side = 1, labels = FALSE)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0),cex.axis = 1.5)
text(x = 1:3,
     y = -0.014,
     labels =  c("CoKriging.Mate\'rn","CoKriging.LMC", "Biv.DeepKriging") ,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1.4) 

# Basic box plot
# par(mfrow=c(2,1))
p1 <- ggplot(df, aes(x=`Biv MSE`, y=combined_MSE)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16,
               outlier.size=2) + stat_summary(fun=mean, geom="point", shape=23, size=4)+
  theme(axis.title.x = element_blank(), legend.text=element_text(size=rel(1)), 
        axis.title=element_text(size=15,face="bold"), 
        axis.text=element_text(size=12)) + coord_flip()
p1
# p2 <- ggplot(df, aes(x=MSE, y=variable2)) + 
#   geom_boxplot(outlier.colour="red", outlier.shape=16,
#                outlier.size=2) +stat_summary(fun=mean, geom="point", shape=23, size=4)+ coord_flip()
# combined <- ggarrange(p1,p2,nrow = 2, ncol = 1,
#                       labels = c("variable 1", "variable 2"))
# combined

load("nonstationary_boundary_case.Rdata")
#write.csv(x = test.data,file = "test_data_nonstationary.csv")
deep_pred = read.table("nonstationary_predictions.csv",header = TRUE, sep = ",")
diff_var1_deep = abs(deep_pred$var1 - test.data$var1)
diff_var2_deep = abs(deep_pred$var2 - test.data$var2)

diff_var1_matern = abs(pred[1:400] - test.data$var1)
diff_var2_matern = abs(pred[401:800] - test.data$var2)
col = rev(brewer.pal(9, "Reds"))
col_ = rep(NA,9)
for(i in 1:9){
  col_[i] = col[9-i + 1]
}
myPalette <- colorRampPalette(col_)

par(mfrow=c(2,2), mar=c(3,3,3,8))
quilt.plot(test.data$x,test.data$y, diff_var1_deep, main =
             "Variable 1 Biv.DeepKriging",nx = 20, ny = 20,zlim = c(0,0.6),col = myPalette(100),
           yaxt="n",xaxt = "n",cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,)
quilt.plot(test.data$x,test.data$y, diff_var2_deep, 
           main = "Variable 2 Biv.DeepKriging",nx = 20, ny = 20,zlim = c(0,1),col = myPalette(100),
           yaxt="n",xaxt = "n",cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,)
quilt.plot(test.data$x,test.data$y, diff_var1_matern, main =
             "Variable 1 CoKriging.Mat\'ern",zlim = c(0,0.6),nx = 20, ny = 20,col = myPalette(100),
           yaxt="n",xaxt = "n",cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2)
quilt.plot(test.data$x,test.data$y, diff_var2_matern, 
           main = "Variable 2 CoKriging.Mat\'ern", zlim = c(0,1),nx = 20, ny = 20, col = myPalette(100),
           yaxt="n",xaxt = "n",cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2)

bounds =  read.table("bounds_nonstationary.csv", header = TRUE,sep = ",")

full.pred.summary<-function(estim.par)
{
  ########## Test set #############
  full.data.coords = rbind(as.matrix(test.data[,c(1,2)]),as.matrix(un.grd.train[,c(1,2)]))
  dist.mat<-rdist(full.data.coords)
  theta<-estim.par
  alpha11<-theta[1]
  nu11<-theta[2]
  s11<-theta[3]
  alpha22<-theta[4]
  nu22<-theta[5]
  s22<-theta[6]
  nug1<-theta[7]
  nug2<-theta[8]
  R<-theta[9]
  R_corr1<-theta[10]
  R_corr2<-theta[11]
  del1 <- theta[12]
  del2 <- theta[13]
  
  nu12 = ((nu11 + nu22) /2) #+ (del1*(1-R_corr1))
  alpha12 = ((alpha11 + alpha22)/2) #+ (del2*(1-R_corr2))
  constant = (sqrt(s11*s22)*R*(gamma(nu12)/sqrt(gamma(nu11)*gamma(nu22)))) /
    ((alpha12^nu12)/sqrt(alpha11^nu11 * alpha22^nu22))
  
  COV11 = s11*matern(dist.mat,sqrt(alpha11),nu11)
  COV22 = s22*matern(dist.mat,sqrt(alpha22),nu22)
  COV12 = constant*matern(dist.mat,sqrt(alpha12),nu12)
  
  NUG1<-diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
  NUG2<-diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
  
  
  
  ############## Inverting C11 ##########
  C<-rbind(cbind(COV11+NUG1,COV12),cbind(t(COV12),COV22+NUG2))
  print(dim(C))
  test.index = c(1:400,1201:1600)
  # test.index = c(1:1190,1201:2390)
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  C.test<-C[test.index,test.index]
  C.train<-C[-test.index,-test.index]
  C.test.train<-C[test.index,-test.index]
  #coh12<-lmc.coh(w=u,a1=a1,nu1=nu1,a2=a2,nu2=nu2,b11=b11,b12=b12,b21=b21,b22=b22)
  ############## Inverting C11 ##########
  #C.train<-C22
  #C.test<-C11
  #C.test.train<-C12
  print(dim(C.test))
  print(dim(C.train))
  print(dim(C.test.train))
  
  
  
  
  
  prediction<-C.test.train%*%solve(C.train)%*%c(un.grd.train$var1,un.grd.train$var2) #Conditional mean of a Multivariate Gaussian
  cond_var_MAT<- C.test - C.test.train%*%solve(C.train)%*%t(C.test.train)
  
  return(list(pred = prediction, conditional_var = cond_var_MAT ))
  
  
}
result = full.pred.summary(fit.Model.full$par)
meanY = result$pred
varY = result$conditional_var
lower_boundVar1 = meanY[1:400] - 1.96*sqrt(diag(varY[1:400,1:400]))
upper_boundVar1 = meanY[1:400] + 1.96*sqrt(diag(varY[1:400,1:400]))
lower_boundVar2 = meanY[401:800] - 1.96*sqrt(diag(varY[401:800,401:800]))
upper_boundVar2 = meanY[401:800] + 1.96*sqrt(diag(varY[401:800,401:800]))

count_var0 = 0
count_var1 = 0
for( i in 1:400){
  if (test.data$var1[i] > lower_boundVar1[i] & test.data$var1[i] < upper_boundVar1[i]) count_var0 = count_var0 + 1
  if (test.data$var2[i] > lower_boundVar2[i] & test.data$var2[i] < upper_boundVar2[i]) count_var1 = count_var1 + 1
}
picp1 = count_var0/400
picp2 = count_var1/400

print(paste(picp1,",",picp2))
#[1] "0.975 , 0.955"

width1 = mean(upper_boundVar1 - lower_boundVar1)
width2 = mean(upper_boundVar2 - lower_boundVar2)
print(paste(width1,",",width2))
#[1] "0.192488897847716 , 0.421105667034545"

r1 = c()
r2 = c()
for( i in 1:dim(test.data)[1]){
  if(test.data$x[i] <0.3 & test.data$y[i] <0.3)r1 = c(r1,i)
  else r2 = c(r2,i)
}
r1 = sample(r1,10)
r2 = sample(r2,10)
r = c(r1,r2)

# variable 1  


sample_data <- data.frame(actual = test.data$var1[r], 
                          interval_lower_matern = lower_boundVar1[r],
                          interval_upper_matern = upper_boundVar1[r],
                          interval_lower_deep = bounds$l1[r],
                          interval_upper_deep = bounds$u1[r])
cols <- c("bivariate DeepKriging"="blue","cokriging"="red")
p1 <- ggplot(data = sample_data, aes(x=1:20)) + 
  geom_point(aes(y = actual)) +
  geom_errorbar(aes(ymin = interval_lower_matern, ymax = interval_upper_matern, colour = "cokriging")) +
  geom_errorbar(aes(ymin = interval_lower_deep, ymax = interval_upper_deep, colour = "bivariate DeepKriging")) + 
  geom_vline(xintercept = 10.5) +
  labs(title = "Prediction interval for variable 1",
       x = "Number of data points",
       y = "Variable 1 values") +
  scale_colour_manual(name="Prediction intervals",values=cols) + 
  annotate("text", x=5, y= -1.7, label= "rough surface")+
  annotate("text", x=15, y= -1.7, label= "smooth surface")+
  theme( legend.text=element_text(size=rel(1.2)), 
        axis.title=element_text(size=15,face="bold"), 
        axis.text=element_text(size=12))
p1
# variable 2


sample_data <- data.frame(actual = test.data$var2[r], 
                          interval_lower_matern = lower_boundVar2[r],
                          interval_upper_matern = upper_boundVar2[r],
                          interval_lower_deep = bounds$l2[r],
                          interval_upper_deep = bounds$u2[r])
cols <- c("bivariate DeepKriging"="blue","cokriging"="red")
p2 <- ggplot(data = sample_data, aes(x=1:20)) + 
  geom_point(aes(y = actual)) +
  geom_errorbar(aes(ymin = interval_lower_matern, ymax = interval_upper_matern, colour = "cokriging")) +
  geom_errorbar(aes(ymin = interval_lower_deep, ymax = interval_upper_deep, colour = "bivariate DeepKriging")) + 
  geom_vline(xintercept = 10.5) +
  labs(title = "Prediction interval for variable 2",
       x = "Number of data points",
       y = "Variable 2 values") +
  scale_colour_manual(name="Prediction intervals",values=cols) + 
  annotate("text", x=5, y= -1.7, label= "rough surface") +
  annotate("text", x=15, y= -1.7, label= "smooth surface") + 
  theme(legend.text=element_text(size=rel(1.2)), 
        axis.title=element_text(size=15,face="bold"), 
        axis.text=element_text(size=12))
  
p2
combined <- ggarrange(p1,p2,nrow = 2, ncol = 1)
combined

##################################### Gaussian ######### all comparison

rm(list = ls())
load("lmc_Stationary.Rdata")

# > mean(MSE_variable1_lmc)
# [1] 0.006314352
# > mean(MSE_variable2_lmc)
# [1] 0.004906759 

df3 = data.frame(variable1 = MSE_variable1_lmc[1,],variable2 = MSE_variable2_lmc[1,])
df3$MSE = "lmc"

mse_deep = read.table("DeepKriging_gaussian_mse.csv", header = T, sep = ",")
df2 = data.frame(variable1 = mse_deep$mse_var1,variable2 = mse_deep$mse_var2)
df2$MSE = "DeepKriging"

df = rbind(df2,df3)


# Basic box plot
# par(mfrow=c(2,1))
p1 <- ggplot(df, aes(x=MSE, y=variable1)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16,
               outlier.size=2) + stat_summary(fun=mean, geom="point", shape=23, size=4)+ coord_flip()
p2 <- ggplot(df, aes(x=MSE, y=variable2)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16,
               outlier.size=2) +stat_summary(fun=mean, geom="point", shape=23, size=4)+ coord_flip()
combined <- ggarrange(p1,p2,nrow = 2, ncol = 1,
                      labels = c("variable 1", "variable 2"))
combined

result = lmc.pred.summary(fit.Model.lmc$par)
meanY = result$pred
varY = result$conditional_var
lower_boundVar1 = meanY[1:120] - 1.96*sqrt(diag(varY[1:120,1:120]))
upper_boundVar1 = meanY[1:120] + 1.96*sqrt(diag(varY[1:120,1:120]))
lower_boundVar2 = meanY[121:240] - 1.96*sqrt(diag(varY[121:240,121:240]))
upper_boundVar2 = meanY[121:240] + 1.96*sqrt(diag(varY[121:240,121:240]))

count_var0 = 0
count_var1 = 0
for( i in 1:120){
  if (test.data$var1[i] > lower_boundVar1[i] & test.data$var1[i] < upper_boundVar1[i]) count_var0 = count_var0 + 1
  if (test.data$var2[i] > lower_boundVar2[i] & test.data$var2[i] < upper_boundVar2[i]) count_var1 = count_var1 + 1
}
picp1 = count_var0/120
picp2 = count_var1/120

print(paste(picp1,",",picp2))

#  [1] "0.958333333333333 , 0.941666666666667"

width1 = mean(upper_boundVar1 - lower_boundVar1)
width2 = mean(upper_boundVar2 - lower_boundVar2)
print(paste(width1,",",width2))
# [1] "0.327825646761294 , 0.287576919895577"

### deepkriging interval prediction time 1200
### %s seconds 665.6777312755585

## picp 
## 0.9333333333333333,0.9583333333333334

## width 
## 0.4951947810419532,0.4485108779967536

rm(list = ls())
load("stationary_matern.Rdata")

# > mean(MSE_variable1_nonstat)
# [1] 0.006482377
# > mean(MSE_variable2_nonstat)
# [1] 0.005013795

result = full.pred.summary(fit.Model.full$par)

meanY = result$pred
varY = result$conditional_var
lower_boundVar1 = meanY[1:120] - 1.96*sqrt(diag(varY[1:120,1:120]))
upper_boundVar1 = meanY[1:120] + 1.96*sqrt(diag(varY[1:120,1:120]))
lower_boundVar2 = meanY[121:240] - 1.96*sqrt(diag(varY[121:240,121:240]))
upper_boundVar2 = meanY[121:240] + 1.96*sqrt(diag(varY[121:240,121:240]))

count_var0 = 0
count_var1 = 0
for( i in 1:120){
  if (test.data$var1[i] > lower_boundVar1[i] & test.data$var1[i] < upper_boundVar1[i]) count_var0 = count_var0 + 1
  if (test.data$var2[i] > lower_boundVar2[i] & test.data$var2[i] < upper_boundVar2[i]) count_var1 = count_var1 + 1
}
picp1 = count_var0/120
picp2 = count_var1/120

print(paste(picp1,",",picp2))

#[1] "0.941666666666667 , 0.958333333333333"

width1 = mean(upper_boundVar1 - lower_boundVar1)
width2 = mean(upper_boundVar2 - lower_boundVar2)
print(paste(width1,",",width2))
# [1] "0.302647755018833 , 0.2766042405168"


#### plot datasets 

### gaussian
df = read.csv("2D_biv_matern_1200.csv")
par(mfrow=c(1,2), mar=c(3,3,3,3))

fields::quilt.plot(df$x,df$y,df$var1, nx = 30, ny = 30, 
                   main = "variable 1",yaxt="n",xaxt = "n")
fields::quilt.plot(df$x,df$y,df$var2, nx = 30, ny = 30,  
                   main = "variable 2", yaxt="n",xaxt = "n")

### nongaussian
df = read.csv("2d_biv_nongaussian_1200.csv")
par(mfrow=c(1,2), mar=c(3,3,3,5))

fields::quilt.plot(df$x,df$y,df$var1, nx = 30, ny = 30, 
                   main = "variable 1",yaxt="n",xaxt = "n")
fields::quilt.plot(df$x,df$y,df$var2, nx = 30, ny = 30,  
                   main = "variable 2", yaxt="n",xaxt = "n")

#### nonstationary 

df = read.csv("2d_biv_nonStationary_1200.csv")
par(mfrow=c(1,2), mar=c(3,3,4,5))

fields::quilt.plot(df$x,df$y,df$var1, nx = 30, ny = 30, 
                   main = "variable 1",yaxt="n",xaxt = "n",cex.main =2.5)
fields::quilt.plot(df$x,df$y,df$var2, nx = 30, ny = 30,  
                   main = "variable 2", yaxt="n",xaxt = "n",cex.main =2.5)
#### real_data 
saudi_data_orig = read.csv("real_data_full.csv")

df1 = data.frame(Z = saudi_data_orig$Z1)
df1$variable = "U component"

df2 = data.frame(Z = saudi_data_orig$Z2)
df2$variable = "V component"
df = rbind(df1,df2)


### density plot
pdf1 <- ggplot(saudi_data_orig) + geom_density(aes(Z1)) +
  labs(x = "U component", y = "Density", size = 10) +
  theme(axis.title=element_text(size=20,face="bold"), axis.text.x = element_text(face="bold", 
        color="#993333", size=14, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=14, angle=45))

pdf1 

pdf2 <- ggplot(saudi_data_orig) + geom_density(aes(Z2)) +
  labs(x = "V component", y = "Density", size = 10) +
  theme(axis.title=element_text(size=20,face="bold"), axis.text.x = element_text(face="bold", 
        color="#993333", size=14, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=14, angle=45))

pdf2 

combined <- ggarrange(pdf1,pdf2,nrow = 1, ncol = 2)
combined


### scatter plot to show correlation 

# sample_saudi = saudi_data_orig[sample(1:506771,100000),]
ggplot(saudi_data_orig, aes(x = Z1, y = Z2)) +
  geom_point(alpha=0.1) +
  labs(x = "U component",y = "V component", size = 15) +
  theme(axis.title=element_text(size=20,face="bold"), axis.text.x = element_text(face="bold", 
                                    color="#993333", size=14, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=14, angle=45))

# Basic box plot
# par(mfrow=c(2,1))
p1 <- ggplot(df, aes(x=variable, y=Z)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=1,
               outlier.size=1) + stat_summary(fun=mean, geom="point", shape=23, size=4)+ coord_flip()

p1

## computation time plots 

df2 <- data.frame(model=rep(c("Biv.DeepKriging point prediction CPU", 
                             "Biv.DeepKriging point prediction GPU",
                             "Biv.DeepKriging interval prediction CPU",
                             "Biv.DeepKriging interval prediction GPU",
                             "CoKriging.Mate\'rn"), each=4),
                  N=rep(c("1200", "2400", "3600", "6400"),5),
                  time=c(63.0,87.8,107.3,150.3,23.1, 35.6, 41.8, 54.2, 
                         154 , 228,293,471,69,103,137,215, 6624,22716,40824,83484))
head(df2)
df2$time = log(df2$time)
head(df2)

p<-ggplot(df2, aes(x=N, y=time, group=model)) +
  geom_line(aes(linetype=model,color=model),size=3)+
  geom_point(aes(shape=model,color=model),size=5)+
  labs(y= "Seconds (in log scale)", x = "Number of locations", size = 7)+ 
  theme(legend.text=element_text(size=rel(1.5)), axis.title=element_text(size=22,face="bold"),
        axis.text=element_text(size=17))
p

df2 <- data.frame(model=rep(c("Deep Learning point prediction CPU", 
                              "Deep Learning point prediction GPU",
                              "Statistical modelling in CPU"), each=4),
                  N=rep(c(2400, 4800, 7200, 12800),3),
                  time=c(63.0,87.8,107.3,150.3,23.1, 35.6, 41.8, 54.2, 
                         6644,22032,40824,83484))
head(df2)

df2$time = log(df2$time)

p<-ggplot(df2, aes(x=N, y=time, group=model)) +
  geom_line(aes(linetype=model,color=model))+
  geom_point(aes(shape=model,color=model))+
  labs(y= "Seconds (in Log scale)", x = "Matrix size", size = 2)+ 
  theme(legend.text=element_text(size=rel(1)), axis.title=element_text(size=14,face="bold"))
p

tools::compactPDF("plots/scatter_plot_real_data.pdf", gs_quality = 'ebook')


## plot of the whole real data 

library(ncdf4)
library(fields)
library(sp)
library(mapdata)
library(dplyr)
library(maptools)
library(geoR)
ncname <- paste("Wind_anomaly_20090101.nc", sep='')
ncin <- nc_open(ncname)

dname1 <- "U_anomaly"
dname2 <- "V_anomaly" 
u_array <- ncvar_get(ncin, dname1)
v_array <- ncvar_get(ncin, dname2)
# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
lat <- ncvar_get(ncin,"lat")
nc_close(ncin)
U <- u_array[,, 1]
V <- v_array[,, 1]
lon.lat <- expand.grid(lon,lat)
lon_new <- matrix(lon.lat[, 1], ncol = length(lat))
lat_new <- matrix(lon.lat[, 2], ncol = length(lat))

test1 <- data.frame(lon.lat, c(U), c(V))
colnames(test1) <- c('lon', 'lat', 'Z1', 'Z2')
spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))


saudi<- map("world", c("Jordan", "Iraq", "Iran","Syria", "Lebanon", "Israel", "Kenya", "Eritrea", "Ethiopia", "South Sudan", "Sudan", "Egypt", "UAE", "Saudi", "Oman", "Yemen", "Somalia", "Djibouti", "Pakistan", "India", "Kuwait"), fill = TRUE)
# saudi<- map("world", c("Saudi"), fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

saudi_data_orig <- data.frame(spdf[is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

saudi_data_orig <- data.frame(spdf)

par(mfrow = c(1,2), mar = c(5,5,3,5.5))
fields::quilt.plot(saudi_data_orig$lon.1,saudi_data_orig$lat.1,saudi_data_orig$Z1,nx = 500, ny = 500, zlim = c(-10,10),
                   main = "U component", add.legend = F,
                   xlab = "Longitude",ylab="Latitude",cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2)
map("worldHires",lwd = 1.5, add = T)

par(mar = c(5,3,3,4))
fields::quilt.plot(saudi_data_orig$lon.1,saudi_data_orig$lat.1,saudi_data_orig$Z2,nx = 500, ny = 500,
                   main = "V component",
                   xlab = "Longitude",ylab="Latitude",cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2)
map("worldHires", lwd = 1.5, add = T)



## NEOM region data 
## resolution 5X5 km

library(ncdf4)
library(fields)
library(sp)
library(mapdata)
library(dplyr)
library(maptools)
library(geoR)
ncname <- paste("Wind_anomaly_20090101.nc", sep='')
ncin <- nc_open(ncname)

dname1 <- "U_anomaly"
dname2 <- "V_anomaly" 
u_array <- ncvar_get(ncin, dname1)
v_array <- ncvar_get(ncin, dname2)
# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
lat <- ncvar_get(ncin,"lat")
nc_close(ncin)
U <- u_array[,, 1]
V <- v_array[,, 1]
lon.lat <- expand.grid(lon,lat)
lon_new <- matrix(lon.lat[, 1], ncol = length(lat))
lat_new <- matrix(lon.lat[, 2], ncol = length(lat))

test1 <- data.frame(lon.lat, c(U), c(V))
colnames(test1) <- c('lon', 'lat', 'Z1', 'Z2')
spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))


saudi<- map("world", c("Jordan", "Iraq", "Iran","Syria", "Lebanon", "Israel", "Kenya", "Eritrea", "Ethiopia", "South Sudan", "Sudan", "Egypt", "UAE", "Saudi", "Oman", "Yemen", "Somalia", "Djibouti", "Pakistan", "India", "Kuwait"), fill = TRUE)
# saudi<- map("world", c("Saudi"), fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

saudi_data_orig <- data.frame(spdf[is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

saudi_data_orig <- data.frame(spdf)

newdata <- saudi_data_orig[which(saudi_data_orig$lon > 34.7 & saudi_data_orig$lon < 35.4 &
                                   saudi_data_orig$lat > 27 & saudi_data_orig$lat < 29 ), ]
par(mfrow = c(2,2), mar = c(3,3,3,4))
par(mar = c(3,2.7,3,4))
par(mar = c(3,3,3,4))
fields::quilt.plot(newdata[,1],newdata[,2],newdata[,3],nx = 15, ny = 40, zlim = c(-3,7),
                   main = "observations: Variable 1",xaxt="n",
                   xlab = "Longitude",ylab="Latitude", add.legend = F,
                   cex.main =2,lwd = 2, cex.axis = 1.2,)
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 1.5, add = T)
points(35.2921,27.9165, pch = ".",col = "black", bg = "yellow", cex = 10)
text(35.2,27.9165, "NEOM", font = 2,col = "black")

fields::quilt.plot(newdata[,1],newdata[,2],newdata[,4],nx = 15, ny = 40,
                   main = "observations: Variable 2",
                   xlab = "Longitude",ylab="Latitude", add.legend = F,
                   cex.main =2,lwd = 2,cex.axis = 1.2,)
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 1.5, add = T)
points(35.2921,27.9165, pch = ".",col = "cadetblue1", bg = "yellow", cex = 10)
text(35.2,27.9165, "NEOM", font = 2,col = "cadetblue1")



df = read.table("real_data_interpolation.csv",header = T, sep = ",")
# par(mfrow=c(1,2), mar=c(3,3,3,3))
fields::quilt.plot(df$x,df$y, df$var1, main =
                     "Biv.DeepKriging interpolation: Variable 1",xaxt = "n",yaxt = "n",
                   nx = 144, ny = 144,zlim = c(-3,7),
                   cex.main =2,lwd = 2,cex.axis = 1.2,legend.cex = 40)
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 1.5, add = T)
points(35.2921,27.9165, pch = ".",col = "black", bg = "yellow", cex = 10)
text(35.2,27.9165, "NEOM", font = 2,col = "black")

fields::quilt.plot(df$x,df$y, df$var2, 
                   main = "Biv.DeepKriging interpolation: Variable 2",yaxt = "n",
                   nx = 144, ny = 144,cex.main =2,lwd = 2,cex.axis = 1.2,legend.cex = 40)
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 1.5, add = T)
points(35.2921,27.9165, pch = ".",col = "cadetblue1", bg = "yellow", cex = 10)
text(35.2,27.9165, "NEOM", font = 2,col = "cadetblue1")








