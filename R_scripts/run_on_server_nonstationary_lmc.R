############################################################################################
library(fields)
library(geoR)

setwd("/home/nagp/Desktop/Ghulam_work/final_folder/Bivariate_DeepKriging/")
exa_data = read.table("synthetic_data_simulations/2d_biv_nonStationary_1200.csv",header = TRUE, sep = ",")

df = do.call(rbind, Map(data.frame,x = exa_data$x, y = exa_data$y, var1=exa_data$var1,var2 = exa_data$var2))
rand_index = c()
for(i in 1:10){
  rand_index = cbind(rand_index,sample(1:1200,800))
  # rand_index = cbind(rand_index,sample(1:1200,10))
}

MSE_variable1_lmc = c()
MSE_variable2_lmc = c()


start.time = Sys.time()
for(i in 1:10){
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
                       maxit=6000))
  }
  
  indmat.estim<-optim_indmat_loglik(init.ind)
  
  #####################################################################
  ############ LMC calculation ########################################
  #####################################################################
  lmc.loglikelihood_allcomponents<-function(p,z,dmat.ml){
    
    
    theta<-p
    a1<-theta[1]
    nu1<-theta[2]
    sigma1<-1
    a2<-theta[3]
    nu2<-theta[4]
    sigma2<-1
    b11<-theta[5]
    b12<-theta[6]
    b21<-theta[7]
    b22<-theta[8]
    nug1<-theta[9]
    nug2<-theta[10]
    ######## Putting hard constraints on the parameters #############
    if( a1<=0 | nu1<=0 | sigma1<=0 | a2<=0 | nu2<=0 | sigma2<=0 |nug1<0|nug2<0)
    {
      return(list(mlv=Inf))
    }
    else
    {
      
      dist.mat<-dmat.ml
      n<-nrow(dist.mat)
      C11<-sigma1*matern(dist.mat,a1,nu1)
      C22<-sigma2*matern(dist.mat,a2,nu2)
      
      COV11<-(b11^2)*C11+(b12^2)*C22
      COV22<-(b21^2)*C11+(b22^2)*C22
      COV12<-(b11*b21)*C11+(b12*b22)*C22
      
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
      
      
      
      Cchol<-chol(C)
      Cinv<-chol2inv(Cchol)
      logD<-determinant(C)$modulus
      
      nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
      if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
      return(list(mlv=nloglikelihood,full.cov=C, Cchol=Cchol))
    }
  }
  
  init.lmc<-c(indmat.estim$par[c(1,2,4,5)],indmat.estim$par[3],0,0,indmat.estim$par[6],indmat.estim$par[7:8])
  ##########################################################################################
  ##### Now we write the code for log_likelihood of Linear model of coregionalization ######
  ##########################################################################################
  
  
  lmc.loglikelihood<-function(par)
  {
    
    return(lmc.loglikelihood_allcomponents(p=par,z=c(un.grd.train$var1,un.grd.train$var2),
                                           dmat.ml=dist.mat.train)$mlv)
  }
  
  
  
  
  ### Finding mle parametrs for lmc model ####
  optim_lmc.loglikelihood <- function(par){
    optim(par=par,
          fn = lmc.loglikelihood,
          hessian=FALSE, 
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=3000))
  }
  
  
  
  fit.Model.lmc <- optim_lmc.loglikelihood(par=init.lmc)
  
  test.data = df[-rand_index[,i],]
  
  lmc.pred.summary<-function(estim.par)
  {
    ########## Test set #############
    full.data.coords = rbind(as.matrix(test.data[,c(1,2)]),as.matrix(un.grd.train[,c(1,2)]))
    dist.mat<-rdist(full.data.coords)
    p<-estim.par
    print(dim(full.data.coords))
    
    theta<-p
    a1<-theta[1]
    nu1<-theta[2]
    sigma1<-1
    a2<-theta[3]
    nu2<-theta[4]
    sigma2<-1
    b11<-theta[5]
    b12<-theta[6]
    b21<-theta[7]
    b22<-theta[8]
    nug1<-theta[9]
    nug2<-theta[10]
    ######## Putting hard constraints on the parameters #############
    #n<-nrow(dist.mat)
    C11<-sigma1*matern(dist.mat,a1,nu1)
    C22<-sigma2*matern(dist.mat,a2,nu2)
    
    COV11<-(b11^2)*C11+(b12^2)*C22
    COV22<-(b21^2)*C11+(b22^2)*C22
    COV12<-(b11*b21)*C11+(b12*b22)*C22
    
    NUG1<-diag(nug1,nrow = nrow(COV11),ncol=ncol(COV11))
    NUG2<-diag(nug2,nrow = nrow(COV22),ncol=ncol(COV22))
    ##################################################################################
    ############## Creating full covariance matrix ###################################
    ##################################################################################
    
    
    
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
    #C.test.train<-C12    print(dim(C.test))
    print(dim(C.train))
    print(dim(C.test.train))
    
    
    
    
    
    prediction<-C.test.train%*%solve(C.train)%*%c(un.grd.train$var1,un.grd.train$var2) #Conditional mean of a Multivariate Gaussian
    cond_var_MAT<- C.test - C.test.train%*%solve(C.train)%*%t(C.test.train)
    
    return(list(pred = prediction, conditional_var = cond_var_MAT ))
    
    
  }
  
  pred = lmc.pred.summary(fit.Model.lmc$par)$pred
  
 MSE_var1_LMC = mean(abs(pred[1:400]-test.data$var1))
 MSE_var2_LMC = mean(abs(pred[401:800]-test.data$var2))
  #MSE_var1_LMC = mean((pred[1:1190]-test.data$var1)^2)
  #MSE_var2_LMC = mean((pred[1191:2380]-test.data$var2)^2)
  
 # MSE_variable1_lmc = cbind(MSE_variable1_lmc,MSE_var1_LMC)
 # MSE_variable2_lmc = cbind(MSE_variable2_lmc,MSE_var2_LMC)
  
 # print(paste(i,"th step done !!!!!!!"))
  
 print(paste(MSE_var1_LMC,",",MSE_var2_LMC))
}



end.time = Sys.time()

time.taken = end.time - start.time
print(time.taken)
save.image("plot_results/lmc_nonStationary.Rdata")



