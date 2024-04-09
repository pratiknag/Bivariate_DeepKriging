library(mapdata)
library(maptools)
library(dplyr)
library(sp)
library(ncdf4)
library(ggplot2)

setwd("/home/nagp/Desktop/Bivariate_DeepKriging/")


## resolution 5X5 km
ncname <- paste("real_data/Wind_anomaly_20090101.nc", sep='')   

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
#Temp <- temp_array[,, 1]
#Z <- z_array[,, 4]

lon.lat <- expand.grid(lon,lat)

lon_new <- matrix(lon.lat[, 1], ncol = length(lat))
lat_new <- matrix(lon.lat[, 2], ncol = length(lat))


###########################   GET ONLY ARABIAN SEA IN STENCHIKOV DATA
par(mfrow=c(1,1))
saudi<- map("world", c("Jordan", "Iraq", "Iran","Syria", "Lebanon", "Israel", "Kenya", "Eritrea", "Ethiopia", "South Sudan", "Sudan", "Egypt", "UAE", "Saudi", "Oman", "Yemen", "Somalia", "Djibouti", "Pakistan", "India", "Kuwait"), fill = TRUE)
# saudi<- map("world", c("Saudi"), fill = TRUE)
IDs <- sapply(strsplit(saudi$names, ":"), function(x) x[1])
saudi <- map2SpatialPolygons(saudi, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

test1 <- data.frame(lon.lat, c(U), c(V))
colnames(test1) <- c('lon', 'lat', 'Z1', 'Z2')
spdf <- SpatialPointsDataFrame(coords = test1[, c("lon", "lat")], data = test1,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))

saudi_data_orig <- data.frame(spdf[is.na(over(spdf, as(saudi, "SpatialPolygons"))), ])

saudi_data_orig <- data.frame(spdf)

dim(saudi_data_orig)

write.csv(saudi_data_orig,"real_data/real_data_full.csv")

rand_index = sample(1:506771,450000)

saudi_data_orig_training <- saudi_data_orig[rand_index,]

write.csv(saudi_data_orig_training,file = "real_data/training_real_dataset_450000.csv")
saudi_data_orig_testing <- saudi_data_orig[-rand_index,]
write.csv(saudi_data_orig_testing,file = "real_data/testing_real_dataset_450000.csv")


df1 = data.frame(Z = saudi_data_orig$Z1)
df1$variable = "variable1"

df2 = data.frame(Z = saudi_data_orig$Z2)
df2$variable = "variable2"
df = rbind(df1,df2)

# Basic box plot
# par(mfrow=c(2,1))
p1 <- ggplot(df, aes(x=variable, y=Z)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=1,
               outlier.size=1) + stat_summary(fun=mean, geom="point", shape=23, size=4)+ coord_flip()
# p2 <- ggplot(df, aes(x=variable, y=Z2)) + 
#   geom_boxplot(outlier.colour="red", outlier.shape=16,
#                outlier.size=2) +stat_summary(fun=mean, geom="point", shape=23, size=4)+ coord_flip()
# combined <- ggarrange(p1,p2,nrow = 2, ncol = 1,
#                       labels = c("variable 1", "variable 2"))
p1

qqnorm(df1$Z, pch = 1, frame = T)
qqline(df1$Z, col = "steelblue", lwd = 2)

qqnorm(df2$Z, pch = 1, frame = T)
qqline(df2$Z, col = "steelblue", lwd = 2)

p <- ggplot(df1, aes(sample = Z))
p + stat_qq() + stat_qq_line()















newdata <- saudi_data_orig[which(saudi_data_orig$lon > 34.7 & saudi_data_orig$lon < 35.4 &
                                   saudi_data_orig$lat > 27 & saudi_data_orig$lat < 29 ), ]
par(mfrow = c(1,1))
fields::quilt.plot(newdata[,1],newdata[,2],newdata[,4],nx = 15, ny = 40,
                   xlab = "Longitude",ylab="Latitude")
map("worldHires", xlim = c(34.7, 35.4), ylim = c(27, 29), lwd = 0.75, add = T)

lon_1 = seq(min(newdata$lon),max(newdata$lon),length.out = 5*29)
lat_1 = seq(min(newdata$lat),max(newdata$lat),length.out = 5*29)
lon.lat <- expand.grid(lon_1,lat_1)

write.csv(lon.lat,"real_data/interpolation_location.csv")


saudi_data_orig_testing <- saudi_data_orig[-rand_index,]

training_loc = data.frame(x = saudi_data_orig_training$lon, y = saudi_data_orig_training$lat)
training_var1 = data.frame(z1 = saudi_data_orig_training$Z1)
training_var2 = data.frame(z2 = saudi_data_orig_training$Z2)

testing_loc = data.frame(x = saudi_data_orig_testing$lon, y = saudi_data_orig_testing$lat)
testing_var1 = data.frame(z1 = saudi_data_orig_testing$Z1)
testing_var2 = data.frame(z2 = saudi_data_orig_testing$Z2)


setwd("/Users/nagp/Desktop/Directed research/research works going on/final_scripts/real_data/")
write.csv(training_loc, file = "training_locations.csv")
write.csv(training_var1, file = "training_variable1.csv")
write.csv(training_var2, file = "training_variable2.csv")

write.csv(testing_loc, file = "testing_locations.csv")
write.csv(testing_var1, file = "testing_variable1.csv")
write.csv(testing_var2, file = "testing_variable2.csv")
write.csv(saudi_data_orig_testing, file = "testing.csv")



#plot(saudi_data_orig[,1:2])
exa_data = read.table("real_data_full_pred.csv", header = TRUE, sep = ",")
par(mfrow=c(1,2))
fields::quilt.plot(saudi_data_orig_training[,1],saudi_data_orig_training[,2],saudi_data_orig_training[,4],zlim = c(-11,11),nx = 600, ny =600,
                   xlab = "Longitude",ylab="Latitude")
map("worldHires", xlim = c(26.719, 85.078), ylim = c(5.625, 42.188), lwd = 0.75, add = T)


  
z = residuals(lm(saudi_data_orig[,3]~saudi_data_orig[,1]+saudi_data_orig[,2]))
fields::quilt.plot(saudi_data_orig[,1],saudi_data_orig[,2],z,nx = 200, ny =200,
                   cex.lab = 2, cex.axis = 2,cex.main =1,xlab = "Longitude",ylab="Latitude")
map("worldHires",lwd = 0.75, add = T)

write.csv(x = saudi_data_orig,file = "real_data_sample.csv")
write.csv(x = saudi_data_orig,file = "real_data_full.csv")



