library(mapdata)
library(maptools)
library(dplyr)
library(sp)
library(ncdf4)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(geoR)
library(fields)

setwd("/home/nagp/Desktop/Ghulam_work/final_folder/Bivariate_DeepKriging/")
df = read.table("real_data/testing_real_dataset_450000.csv",header = TRUE, sep = ",")

fields::quilt.plot(df$lon.1, df$lat.1, df$Z2)
rand_index = sample(1:dim(df)[1],1000)
test.data = df[rand_index,]
files_data <- list.files(path="splitted_data/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
files_result <- list.files(path="estimation/", pattern="*.Rdata", full.names=TRUE, recursive=FALSE)

  MSE_var1 = c()
  MSE_var2 = c()
  lb1 = c()
  ub1 = c()
  lb2 = c()
  ub2 = c()
  
  predictions1 = c()
  predictions2 = c()
  
par(mfrow=c(1,2))
fields::quilt.plot(test.data$lon.1, test.data$lat.1, test.data$Z2)
fields::quilt.plot(test.data$lon.1, test.data$lat.1, predictions2)
  
  for(i in 1:1000){
    index = test.data[i,]$name
    test_data = test.data[i,]
    data_path =  paste0("splitted_data/dataset_",index-1,".csv")
    un.grd.train = read.table(data_path, header = TRUE, sep = ",")
    un.grd.train = un.grd.train[sample(1:dim(un.grd.train)[1],800),]
    
    est_path = paste0("estimation/",index,".Rdata")
    fit.Model.full = readRDS(est_path)
    
    full.pred.summary<-function(estim.par)
    {
      ########## Test set #############
      full.data.coords = rbind(as.matrix(test_data[1,][,c(3,4)]),as.matrix(un.grd.train[,c(3,4)]))
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
      test.index = c(1,801)
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
      # print(dim(C.test))
      # print(dim(C.train))
      # print(dim(C.test.train))
      
      
      
      
      
      prediction<-C.test.train%*%solve(C.train)%*%c(un.grd.train$Z1,un.grd.train$Z2) #Conditional mean of a Multivariate Gaussian
      cond_var_MAT<- C.test - C.test.train%*%solve(C.train)%*%t(C.test.train)
      
      return(list(pred = prediction, conditional_var = cond_var_MAT ))
      
      
    }
    result = full.pred.summary(fit.Model.full$par)
    pred = result$pred
    predictions1 = cbind(predictions1,pred[1])
    predictions2 = cbind(predictions2,pred[2])
    varY = result$conditional_var
    lower_boundVar1 = pred[1] - 1.96*sqrt(varY[1,1])
    upper_boundVar1 = pred[1] + 1.96*sqrt(varY[1,1])
    lower_boundVar2 = pred[2] - 1.96*sqrt(varY[2,2])
    upper_boundVar2 = pred[2] + 1.96*sqrt(varY[2,2])
    lb1 = cbind(lb1,lower_boundVar1)
    ub1 = cbind(ub1,upper_boundVar1)
    lb2 = cbind(lb2,lower_boundVar2)
    ub2 = cbind(ub2,upper_boundVar2)

    mse_var1 = (pred[1]-test_data$Z1)^2
    mse_var2 = (pred[2]-test_data$Z2)^2
    
    MSE_var1 = cbind(MSE_var1,mse_var1)
    MSE_var2 = cbind(MSE_var2,mse_var2)
    
    print(paste(i,"th step done !!!!!!!!!!!!!!!!!!!!!!!!!!"))
  }
  count_var0 = 0
  count_var1 = 0
  for( i in 1:100){
    if (test.data$Z1[i] > lb1[i] & test.data$Z1[i] < ub1[i]) count_var0 = count_var0 + 1
    if (test.data$Z2[i] > lb2[i] & test.data$Z2[i] < ub2[i]) count_var1 = count_var1 + 1
  }
  picp1 = count_var0/100
  picp2 = count_var1/100
  
  print(paste(picp1,",",picp2))
  
  # [1] "0.6 , 0.12"
  
  width1 = mean(upper_boundVar1 - lower_boundVar1)
  width2 = mean(upper_boundVar2 - lower_boundVar2)
  print(paste(width1,",",width2))
  # "1.67128870091983 , 1.34285631205146"

## Prameter estimation time difference 

#Time difference of 2.18597 days
  
## deepkriging 450000 data points time taken for point prediction 
  ## %s seconds 3373.135583639145
  
## deepkriging 147000 data points interval prediction 
  ## %s seconds 33098.696251153946
  
## picp and width 
## 0.970000,0.950000
##  1.222264529004371,1.3361754733446265

  
## deepkriging 450000 data points time taken for interval
  
saveRDS(MSE_var1,"real_MSE_var1.rds")
saveRDS(MSE_var2,"real_MSE_var2.rds")

mse_var1 = readRDS("real_MSE_var1.rds")
mse_var2 = readRDS("real_MSE_var2.rds")

df2 = read.table("real_data/real_data_147000_mse_deep.csv", header = TRUE, sep = ",")

df1 = data.frame(variable1 = sqrt(mse_var1[1,]),variable2 = sqrt(mse_var2[1,]))
df1$standard_error = "cokriging"
df2 = data.frame(variable1 = sqrt(df2$mse1),variable2 = sqrt(df2$mse2))
df2$standard_error = "DeepKriging_147000"


## 450000 data points estimation result

df3 = read.table("real_data/real_data_450000_mse_deep.csv", header = TRUE, sep = ",")
df3 = df3[sample(1:56771,1000),]
df3 = data.frame(variable1 = sqrt(df3$mse1),variable2 = sqrt(df3$mse2))
df3$standard_error = "DeepKriging_450000"
df = rbind(df1,df2,df3)

# Basic box plot
p1 <- ggplot(df, aes(x=standard_error, y=variable1)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16,
               outlier.size=2) + stat_summary(fun=mean, geom="point", shape=23, size=4) + coord_flip()

p2 <- ggplot(df, aes(x=standard_error, y=variable2)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16,
               outlier.size=2) +stat_summary(fun=mean, geom="point", shape=23, size=4) + coord_flip()


combined <- ggarrange(p1,p2,nrow = 2, ncol = 1,
                      labels = c("variable 1", "variable 2"))
combined

## prediction results 
df1 = read.table("training_real_dataset_450000.csv",header = T, sep = ",")
df1 = do.call(rbind, Map(data.frame,x = df1$lon, y = df1$lat, var1=df1$Z1,var2 = df1$Z2))
df2 = read.table("real_data_450000_predictions.csv", header = T, sep = ",")
df2 = df2[,2:5]
df = rbind(df1,df2)

par(mfrow=c(1,2), mar=c(3,3,3,3))
fields::quilt.plot(df$x,df$y, df$var1, main =
             "variable 1 predictions biv DeepKriging",nx = 600, ny = 600,zlim = c(-11,11))
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)
fields::quilt.plot(df$x,df$y, df$var2, 
           main = "variable 2 predictions biv DeepKriging",nx = 600, ny = 600,zlim = c(-11,15))
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)

## NEOM region data 
## resolution 5X5 km
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

saudi_data_orig <- data.frame(spdf[is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

saudi_data_orig <- data.frame(spdf)

newdata <- saudi_data_orig[which(saudi_data_orig$lon > 34.7 & saudi_data_orig$lon < 35.4 &
                                   saudi_data_orig$lat > 27 & saudi_data_orig$lat < 29 ), ]
par(mfrow = c(2,2), mar = c(3,3,3,3))
fields::quilt.plot(newdata[,1],newdata[,2],newdata[,3],nx = 15, ny = 40, zlim = c(-3,7),
                   main = "observations: Variable 1",
                   xlab = "Longitude",ylab="Latitude")
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 0.75, add = T)
points(35.2921,27.9165, pch = ".",col = "black", bg = "yellow", cex = 10)
text(35.2,27.9165, "NEOM", font = 2,col = "black")
fields::quilt.plot(newdata[,1],newdata[,2],newdata[,4],nx = 15, ny = 40,
                   main = "observations: Variable 2",
                   xlab = "Longitude",ylab="Latitude")
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 0.75, add = T)
points(35.2921,27.9165, pch = ".",col = "cadetblue1", bg = "yellow", cex = 10)
  text(35.2,27.9165, "NEOM", font = 2,col = "cadetblue1")



df = read.table("real_data_interpolation.csv",header = T, sep = ",")
# par(mfrow=c(1,2), mar=c(3,3,3,3))
fields::quilt.plot(df$x,df$y, df$var1, main =
                     "variable 1 interpolation biv DeepKriging",nx = 144, ny = 144,zlim = c(-3,7),)
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 0.75, add = T)
points(35.2921,27.9165, pch = ".",col = "black", bg = "yellow", cex = 10)
text(35.2,27.9165, "NEOM", font = 2,col = "black")
fields::quilt.plot(df$x,df$y, df$var2, 
                   main = "variable 2 interpolation biv DeepKriging",nx = 144, ny = 144)
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 0.75, add = T)
points(35.2921,27.9165, pch = ".",col = "cadetblue1", bg = "yellow", cex = 10)
text(35.2,27.9165, "NEOM", font = 2,col = "cadetblue1")

