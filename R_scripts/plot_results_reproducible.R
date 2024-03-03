rm(list = ls())
## nongaussian ########### LMC best
library(ggplot2)
library(fields)
library(ggpubr)
library(RColorBrewer)
library(geoR)
library(latex2exp)
file_path = "/Users/nagp/Desktop/Biv.DeepKriging/Bivariate_DeepKriging/plot_results"
setwd(file_path)

### boxplot for nonGaussian with covariates
mse_deep = read.csv("DeepKriging_nongaussian-cov_mse.csv", header = T)
mse_lmc = read.csv("validation_LMC_cov.csv", header = T)
mse_par_matern = read.csv("validation_par_Matern_cov.csv", header = T)
par(mfrow = c(1,2))
##### variable 1
par(mar = c(7,4,1.5,1))
colors = c("wheat4","pink3","lightskyblue1")
title = "Variable 1"
boxplot(mse_par_matern$mse_var1, mse_lmc$mse_var1, mse_deep$mse_var1,
        main = title , col = colors,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,
        xaxt = "n", yaxt = "n", outline=FALSE)
axis(side = 1, labels = FALSE)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0),cex.axis = 1.5)
text(x = 1:3,
     y = 0,
     labels =  c("CoKriging.Mate\'rn","CoKriging.LMC", "Biv.DeepKriging") ,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1.4)  
## Variable 2
par(mar = c(7,4,1.5,1))
colors = c("wheat4","pink3","lightskyblue1")
title = "Variable 2"
boxplot(mse_par_matern$mse_var2, mse_lmc$mse_var2, mse_deep$mse_var2,
        main = title , col = colors,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,
        xaxt = "n", yaxt = "n", outline=FALSE)
axis(side = 1, labels = FALSE)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0),cex.axis = 1.5)
text(x = 1:3,
     y = 0,
     labels =  c("CoKriging.Mate\'rn","CoKriging.LMC", "Biv.DeepKriging") ,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1.4)  


### boxplot for nonstationary
mse_deep = read.csv("DeepKriging_nonstationary_mse.csv", header = T)
mse_lmc = read.csv("validation_LMC_nonstat.csv", header = T)
mse_par_matern = read.csv("validation_par_Matern_nonstat.csv", header = T)
par(mfrow = c(1,2))
##### variable 1
par(mar = c(7,4,1.5,1))
colors = c("wheat4","pink3","lightskyblue1")
title = "Variable 1"
boxplot(mse_par_matern$mse_var1, mse_lmc$mse_var1, mse_deep$mse_var1,
        main = title , col = colors,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,
        xaxt = "n", yaxt = "n", outline=FALSE)
axis(side = 1, labels = FALSE)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0),cex.axis = 1.5)
text(x = 1:3,
     y = 0,
     labels =  c("CoKriging.Mate\'rn","CoKriging.LMC", "Biv.DeepKriging") ,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1.4)  
## Variable 2
par(mar = c(7,4,1.5,1))
colors = c("wheat4","pink3","lightskyblue1")
title = "Variable 2"
boxplot(mse_par_matern$mse_var2, mse_lmc$mse_var2, mse_deep$mse_var2,
        main = title , col = colors,
        cex.lab = 2, cex.axis = 2,cex.main =2,lwd = 2,
        xaxt = "n", yaxt = "n", outline=FALSE)
axis(side = 1, labels = FALSE)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0),cex.axis = 1.5)
text(x = 1:3,
     y = 0,
     labels =  c("CoKriging.Mate\'rn","CoKriging.LMC", "Biv.DeepKriging") ,
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1.4)  

