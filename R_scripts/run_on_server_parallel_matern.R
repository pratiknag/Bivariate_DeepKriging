###### Bivariate stationary ##############
library(fields)
library(geoR)
library(parallel)
setwd("/home/nagp/Desktop/Ghulam_work/final_folder/Bivariate_DeepKriging/")
#load("nongaussian.Rdata")
#exa_data = read.table("2d_biv_nongaussian_1200_test2.csv",header = TRUE, sep = ",")

#MSE_combined = matrix(data=NA,nrow = 1,ncol = 10)
MSE_1 = matrix(data=NA,nrow = 1,ncol = 50)
MSE_2 = matrix(data=NA,nrow= 1, ncol= 50)

rand_index = c()
for(i in 1:50){
  rand_index = cbind(rand_index,sample(1:1200,800))
  #rand_index = cbind(rand_index,sample(1:1200,10))
}
print("inside mclapply")
start.time = Sys.time()
numCores = detectCores()
fx = function(i){
  exa_data = read.table(paste0("synthetic_data_simulations/2d_nongaussian_1200_",i,".csv"),header = TRUE, sep = ",")
  #print(paste0("synthetic_data_simulations/2d_gaussian_1200_",i,".csv")
  df = do.call(rbind, Map(data.frame,x = exa_data$x, y = exa_data$y, var1=exa_data$var1,var2 = exa_data$var2))
  un.grd.train = df[rand_index[,i],]
  dist.mat.train<-rdist(un.grd.train[,-c(3,4)])
  
  mle_ind_mat<-function(p,z,dmat.ml)
  {
    a1<-p[1]
    nu1<-p[2]
    sigma1<-p[3]
    a2<-p[4]
    nu2<-p[5]
    sigma2<-p[6]
    nug1<-p[7]
    nug2<-p[8]
    if(sum(p[1:6]<0)!=0||nug1<0||nug2<0)
    {
      nloglikelihood<-10000000
      return(list(mlv=nloglikelihood,params=NULL))
    }
    else
    {
      
      dist.mat<-dmat.ml
      C11<-sigma1*matern(dist.mat,a1,nu1)
      C22<-sigma2*matern(dist.mat,a2,nu2)
      
      COV12<-matrix(0,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      
      NUG1<-diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      NUG2<-diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      
      ############## Inverting C11 ##########
      C<-rbind(cbind(C11+NUG1,COV12),cbind(t(COV12),C22+NUG2))
      Cchol<-chol(C)
      Cinv<-chol2inv(Cchol)
      logD<-determinant(C)$modulus
      
      nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
      if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
      return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,sigma1=sigma1,sigma2=sigma2,full.cov=C))
    }
    
  }
  
  
  init.ind<-c(1,1,1,1,1,1,0,0)
  
  
  mle_ind_mlv<-function(pars)
  {
    return(mle_ind_mat(p=pars,z=c(un.grd.train$var1,un.grd.train$var2),dmat.ml=dist.mat.train)$mlv)
  }
  
  optim_indmat_loglik <- function(par){
    optim(par=par,
          fn = mle_ind_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=700))
  }
  
  indmat.estim<-optim_indmat_loglik(init.ind)
  
  #####################################################################
  ############ LMC calculation ########################################
  #####################################################################
  fullMat.loglikelihood_allcomponents<-function(p,z,dmat.ml){
    
    
    theta<-p
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
    
    
    
    ######## Putting hard constraints on the parameters #############
    if( alpha11<=0 | nu11<=0 | s11<=0 | alpha22<=0 | nu22<=0 | s22<=0 |
        nug1<0|nug2<0 | del1 <=0 |del2<=0 | R < -1 | R >1 | R_corr1<0 | R_corr1>1 |
        R_corr2<0 | R_corr2>1)
    {
      return(list(mlv=Inf))
    }
    else
    {
      nu12 = ((nu11 + nu22) /2) + (del1*(1-R_corr1))
      alpha12 = ((alpha11 + alpha22)/2) + (del2*(1-R_corr2))
      constant = (sqrt(s11*s22)*R*(gamma(nu12)/sqrt(gamma(nu11)*gamma(nu22)))) /
        ((alpha12^nu12)/sqrt(alpha11^nu11 * alpha22^nu22))
      dist.mat<-dmat.ml
      n<-nrow(dist.mat)
      
      COV11 = s11*matern(dist.mat,sqrt(alpha11),nu11)
      COV22 = s22*matern(dist.mat,sqrt(alpha22),nu22)
      COV12 = constant*matern(dist.mat,sqrt(alpha12),nu12)
      
      NUG1<-diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      NUG2<-diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      # ##################################################################################
      # ############## Creating full covariance matrix ###################################
      # ##################################################################################
      # COV11op<-COV22op<-COV12op<-matrix(NA,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      # 
      # COV11op[index.mat.ml]<-COV11[index.val.ml]
      # COV22op[index.mat.ml]<-COV22[index.val.ml]
      # COV12op[index.mat.ml]<-COV12[index.val.ml]
      
      
      
      ############## Inverting C11 ##########
      C<-rbind(cbind(COV11+NUG1,COV12),cbind(t(COV12),COV22+NUG2))
      
      
      
      H = eigen(C)
      e_vec = H$vectors
      D = H$values
      D = 1/D
      for(i in 1:length(D)){
        if(D[i] > 1e+05) D[i] = 0
      }
      D.inv = diag(D)
      Cinv<-e_vec%*%D.inv%*%t(e_vec)
      logD<-determinant(C)$modulus
      
      nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
      if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
      return(list(mlv=nloglikelihood,full.cov=C))
    }
  }
  
  init.full<-c(indmat.estim$par[1:8],0,0,0,0.01,0.01)
  ##########################################################################################
  ##### Now we write the code for log_likelihood of Linear model of coregionalization ######
  ##########################################################################################
  
  
  Full.loglikelihood<-function(par)
  {
    
    return(fullMat.loglikelihood_allcomponents(p=par,z=c(un.grd.train$var1,un.grd.train$var2),
                                               dmat.ml=dist.mat.train)$mlv)
  }
  
  
  
  
  ### Finding mle parametrs for lmc model ####
  optim_full.loglikelihood <- function(par){
    optim(par=par,
          fn = Full.loglikelihood,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=1500))
  }
  
  
  
  fit.Model.full <- optim_full.loglikelihood(par=init.full)
  test.data = df[-rand_index[,i],]
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
    
    nu12 = ((nu11 + nu22) /2) + (del1*(1-R_corr1))
    alpha12 = ((alpha11 + alpha22)/2) + (del2*(1-R_corr2))
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
    test.index = c(1:dim(test.data)[1],(dim(df)[1]+1):(dim(df)[1]+dim(test.data)[1]))
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
    
    H = eigen(C.train)
    e_vec = H$vectors
    D = H$values
    D = 1/D
    for(i in 1:length(D)){
      if(D[i] > 1e+05) D[i] = 0
    }
    D.inv = diag(D)
    Cinv<-e_vec%*%D.inv%*%t(e_vec)
    
    
    
    
    
    prediction<-C.test.train%*%Cinv%*%c(un.grd.train$var1,un.grd.train$var2) #Conditional mean of a Multivariate Gaussian
    cond_var_MAT<- C.test - C.test.train%*%Cinv%*%t(C.test.train)
    
    return(list(pred = prediction, conditional_var = cond_var_MAT ))
    
    
  }
  
  pred = full.pred.summary(fit.Model.full$par)$pred
  
  MSE_var1_LMC = mean((pred[1:120]-test.data$var1)^2)
  MSE_var2_LMC = mean((pred[121:240]-test.data$var2)^2)
  # MSE_var1_LMC = mean((pred[1:1190]-test.data$var1)^2)
  # MSE_var2_LMC = mean((pred[1191:2380]-test.data$var2)^2)
  
  MSE_variable1_nonstat = cbind(MSE_variable1_nonstat,MSE_var1_LMC)
  MSE_variable2_nonstat = cbind(MSE_variable2_nonstat,MSE_var2_LMC)
  print(paste(MSE_var1_LMC,",",MSE_var2_LMC))
  print(paste(i,"th step done !!!!!!!"))
  
}

index = seq(1,50)

mclapply(index,fx, mc.cores = numCores - 2)

# par(mfrow=c(1,2))
# quilt.plot(df$x,df$y, df$var1, main =
#              "variable 1")
# quilt.plot(df$x,df$y, df$var2,main = "variable 2")


end.time = Sys.time()

time.taken_nonstat = end.time - start.time
print(time.taken_nonstat)
save.image("plot_results/nongaussian_gaussian.Rdata")
