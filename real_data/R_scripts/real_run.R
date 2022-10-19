###### Bivariate stationary ###############
library(fields)
library(geoR)
setwd("/home/nagp/Desktop/Ghulam_work/final_folder/Bivariate_DeepKriging/")
start.time = Sys.time()
for(i in 1:100){
  file_name = paste0("real_data/splitted_data/dataset_",i-1,".csv")
  exa_data = read.table(file_name,header = TRUE, sep = ",")
  print(exa_data[1,])
	#break
  df = do.call(rbind, Map(data.frame,x = exa_data$lon, y = exa_data$lat, var1=exa_data$Z1,var2 = exa_data$Z2))
  un.grd.train = df[sample(1:dim(df)[1],900),]
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
                       maxit=400))
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
   R_corr1<-0.1
   R_corr2<- 0.1
   del1 <-0.1
   del2 <-0.1
    
    
    
    ######## Putting hard constraints on the parameters #############
    if( alpha11<=0 | nu11<=0 | s11<=0 | alpha22<=0 | nu22<=0 | s22<=0 |
        nug1<0|nug2<0 | del1 <=0 |del2<=0 | R < -1 | R >1 | R_corr1<0 | R_corr1>1 |
        R_corr2<0 | R_corr2>1)
    {
      return(list(mlv=Inf))
    }
    else
    {
      nu12 = ((nu11 + nu22) /2) # + (del1*(1-R_corr1))
      alpha12 = ((alpha11 + alpha22)/2) # + (del2*(1-R_corr2))
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
  
  init.full<-c(indmat.estim$par[1:8],0)
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
                       maxit=800))
  }
  
  
  
  fit.Model.full <- optim_full.loglikelihood(par=init.full)
  file_name = paste0("real_data/estimation/",i,".Rdata")
  saveRDS(fit.Model.full,file_name)


}
# par(mfrow=c(1,2))
# quilt.plot(df$x,df$y, df$var1, main =
#              "variable 1")
# quilt.plot(df$x,df$y, df$var2,main = "variable 2")


end.time = Sys.time()

time.taken_nonstat = end.time - start.time
print(time.taken_nonstat)
#save.image("stationary_matern.Rdata")
