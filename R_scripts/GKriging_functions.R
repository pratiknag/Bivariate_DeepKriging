#file_path = "/home/nagp/Bivariate_DeepKriging/"
#setwd(file_path)
library(geoR)
library(MASS)
library(fields)
mle_ind_mat=function(p,data_train)
{
  if(dim(data_train)[2] > 4){ un.grd.train = data_train[,c(6:9)]}
  else{un.grd.train = data_train[,c(1:4)]}
  dist.mat = rdist(un.grd.train[,-c(3,4)])
  z=c(un.grd.train$var1, un.grd.train$var2)
  a1=p[1]
  nu1=p[2]
  sigma1=p[3]
  a2=p[4]
  nu2=p[5]
  sigma2=p[6]
  nug1=p[7]
  nug2=p[8]
  if(sum(p[1:6]<0)!=0||nug1<0||nug2<0)
  {
    nloglikelihood=10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  {
    C11=sigma1*matern(dist.mat,a1,nu1)
    C22=sigma2*matern(dist.mat,a2,nu2)
    
    COV12=matrix(0,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    NUG1=diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    NUG2=diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    ############## Inverting C11 ##########
    C=rbind(cbind(C11+NUG1,COV12),cbind(t(COV12),C22+NUG2))
    Cchol=chol(C)
    Cinv=chol2inv(Cchol)
    logD=determinant(C)$modulus
    
    nloglikelihood =(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood = 1e+08}
    return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,
                sigma1=sigma1,sigma2=sigma2,full.cov=C, df_train=un.grd.train))
  }
  
}


mle_ind_mlv=function(pars)
{
  return(mle_ind_mat(p=pars,data_train)$mlv)
}

optim_indmat_loglik = function(par){
  optim(par=par,
        fn = mle_ind_mlv,
        hessian=FALSE, 
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=50))
}

########################################################################
###################### FUNCTIONS FOR LMC ###############################
########################################################################

lmc.loglikelihood_allcomponents=function(p,data_train){
  ############# Covariate X ##############
  if(dim(data_train)[2] > 4){
    un.grd.train = data_train[,c(6:9)]
    a = data.matrix(data_train[,c(1:5)])
    b = matrix(0, dim(a)[1], dim(a)[2])
    X = rbind(cbind(a,b),cbind(b,a)) 
  }
  else{
    un.grd.train = data_train[,c(1,2,3,4)]
    X = matrix(c(rep(1,dim(un.grd.train)[1]),rep(0,dim(un.grd.train)[1]),
                 rep(0,dim(un.grd.train)[1]),rep(1,dim(un.grd.train)[1])), 
               nrow = 2*(dim(un.grd.train)[1]),
               ncol = 2, byrow = F)
  }
  
  dist.mat =fields::rdist(un.grd.train[,-c(3,4)])
  z=c(un.grd.train$var1 , un.grd.train$var2)
  theta=p
  a1=theta[1]
  nu1=theta[2]
  sigma1=theta[3]
  a2=theta[4]
  nu2=theta[5]
  sigma2=theta[6]
  b11=theta[7]
  b12=theta[8]
  b21=theta[9]
  b22=theta[10]
  nug1=theta[11]
  nug2=theta[12]
  ######## Putting hard constraints on the parameters #############
  if( a1<=0 | nu1<=0 | sigma1<=0 | a2<=0 | nu2<=0 | sigma2<=0 |nug1<0|nug2<0)
  {
    return(list(mlv=Inf))
  }
  else
  {
    n=nrow(dist.mat)
    C11=sigma1*matern(dist.mat,a1,nu1)
    C22=sigma2*matern(dist.mat,a2,nu2)
    
    COV11=(b11^2)*C11+(b12^2)*C22
    COV22=(b21^2)*C11+(b22^2)*C22
    COV12=(b11*b21)*C11+(b12*b22)*C22
    
    NUG1=diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    NUG2=diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    
    
    ############## Inverting C11 ##########
    C=rbind(cbind(COV11+NUG1,COV12),cbind(t(COV12),COV22+NUG2))
    
    # print("I am here 1")
    
    ## Inverse calculation of C
    H = eigen(C)
    e_vec = H$vectors
    D = H$values
    D = 1/D
    for(i in 1:length(D)){
      if(D[i] > 1e+05) D[i] = 0
    }
    D.inv = diag(D)
    Cinv = e_vec%*%D.inv%*%t(e_vec)
    #################
    
    logD=determinant(C)$modulus
    # print("I am here 2")
    
    ## Inverse calculation of V
    H = eigen(t(X)%*%Cinv%*%X)
    e_vec = H$vectors
    D = H$values
    D = 1/D
    for(i in 1:length(D)){
      if(D[i] > 1e+05) D[i] = 0
    }
    D.inv = diag(D)
    V  = e_vec%*%D.inv%*%t(e_vec)
    #################
    
    beta = V%*%t(X)%*%Cinv%*%z
    # print("I am here 3")
    
    nloglikelihood =(0.5 * logD + 0.5 * t(z - X%*%beta) %*% Cinv %*% (z - X%*%beta)+0.5*length(z)*log(2*pi))
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood = 1e+08}
    return(list(mlv=nloglikelihood,full.cov=C))
  }
}

lmc.loglikelihood=function(par)
{
  
  return(lmc.loglikelihood_allcomponents(p=par,data_train)$mlv)
}

optim_lmc.loglikelihood = function(par){
  optim(par=par,
        fn = lmc.loglikelihood,
        hessian=FALSE, 
        control=list(trace=4,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=300))
}

lmc.pred.summary=function(estim.par, data_train, test_data)
{
  if(dim(data_train)[2] > 4){
    un.grd.train = data_train[,c(6:9)]
    a = data.matrix(data_train[,c(1:5)])
    b = matrix(0, dim(a)[1], dim(a)[2])
    X = rbind(cbind(a,b),cbind(b,a))    
    test.data = test_data[,c(6:9)]
    a = data.matrix(test_data[,c(1:5)])
    b = matrix(0, dim(a)[1], dim(a)[2])
    X0 = rbind(cbind(a,b),cbind(b,a))
  }
  else{
    un.grd.train = data_train[,c(1,2,3,4)]
    X = matrix(c(rep(1,dim(un.grd.train)[1]),rep(0,dim(un.grd.train)[1]),
                 rep(0,dim(un.grd.train)[1]),rep(1,dim(un.grd.train)[1])), 
               nrow = 2*(dim(un.grd.train)[1]),
               ncol = 2, byrow = F)
    test.data = test_data[,c(1,2,3,4)]
    X0 = matrix(c(rep(1,dim(test.data)[1]),rep(0,dim(test.data)[1]),
                 rep(0,dim(test.data)[1]),rep(1,dim(test.data)[1])), 
               nrow = 2*(dim(test.data)[1]),
               ncol = 2, byrow = F)
  }
  full.data.coords = rbind(as.matrix(test.data[,c(1,2)]),as.matrix(un.grd.train[,c(1,2)]))
  dist.mat=rdist(full.data.coords)
  p=estim.par
  print(dim(full.data.coords))
  
  theta=p
  a1=theta[1]
  nu1=theta[2]
  sigma1=theta[3]
  a2=theta[4]
  nu2=theta[5]
  sigma2=theta[6]
  b11=theta[7]
  b12=theta[8]
  b21=theta[9]
  b22=theta[10]
  nug1=theta[11]
  nug2=theta[12]
  
  C11=sigma1*matern(dist.mat,a1,nu1)
  C22=sigma2*matern(dist.mat,a2,nu2)
  
  COV11=(b11^2)*C11+(b12^2)*C22
  COV22=(b21^2)*C11+(b22^2)*C22
  COV12=(b11*b21)*C11+(b12*b22)*C22
  
  NUG1=diag(nug1,nrow = nrow(COV11),ncol=ncol(COV11))
  NUG2=diag(nug2,nrow = nrow(COV22),ncol=ncol(COV22))
  
  ############## Inverting C11 ##########
  C=rbind(cbind(COV11+NUG1,COV12),cbind(t(COV12),COV22+NUG2))
  print(dim(C))
  test.index = c(1:dim(test.data)[1],
                 (dim(un.grd.train)[1]+dim(test.data)[1]+1):(dim(un.grd.train)[1]+
                                                            dim(test.data)[1]+dim(test.data)[1]))
  
  C.test=C[test.index,test.index]
  C.train=C[-test.index,-test.index]
  C.test.train=C[test.index,-test.index]
  
  print(dim(C.train))
  print(dim(C.test.train))
  
  z=c(un.grd.train$var1 , un.grd.train$var2)
  
  
  
  ## Inverse of C.train
  H = eigen(C.train)
  e_vec = H$vectors
  D = H$values
  D = 1/D
  for(i in 1:length(D)){
    if(D[i] > 1e+05) D[i] = 0
  }
  D.inv = diag(D)
  C.train.inv<-e_vec%*%D.inv%*%t(e_vec)
  ##################
  
  ## Inverse calculation of V
  H = eigen(t(X)%*%C.train.inv%*%X)
  e_vec = H$vectors
  D = H$values
  D = 1/D
  for(i in 1:length(D)){
    if(D[i] > 1e+05) D[i] = 0
  }
  D.inv = diag(D)
  V  = e_vec%*%D.inv%*%t(e_vec)
  #################
  
  beta = V%*%t(X)%*%C.train.inv%*%z
  
  
  prediction=X0%*%beta + C.test.train%*%(C.train.inv)%*%(z - X%*%beta) #Conditional mean of a Multivariate Gaussian
  cond_var_MAT= C.test - C.test.train%*%(C.train.inv)%*%t(C.test.train)
  
 # prediction = prediction + c(rep(mean(un.grd.train$var1),dim(test.data)[1]),
  #                             rep(mean(un.grd.train$var2),dim(test.data)[1]))
  return(list(pred = prediction, conditional_var = cond_var_MAT))
  
  
}

########################################################################
############## FUNCTIONS FOR PARSIMONIOUS MATERN #######################
########################################################################

par_Matern.loglikelihood_allcomponents=function(p,data_train){
  ############# Covariate X ##############
  if(dim(data_train)[2] > 4){
    un.grd.train = data_train[,c(6:9)]
    a = data.matrix(data_train[,c(1:5)])
    b = matrix(0, dim(a)[1], dim(a)[2])
    X = rbind(cbind(a,b),cbind(b,a)) 
  }
  else{
    un.grd.train = data_train[,c(1,2,3,4)]
    X = matrix(c(rep(1,dim(un.grd.train)[1]),rep(0,dim(un.grd.train)[1]),
                 rep(0,dim(un.grd.train)[1]),rep(1,dim(un.grd.train)[1])), 
               nrow = 2*(dim(un.grd.train)[1]),
               ncol = 2, byrow = F)
  }
  
  dist.mat =fields::rdist(un.grd.train[,-c(3,4)])
  z=c(un.grd.train$var1 , un.grd.train$var2)
  theta=p
  a1=theta[1]
  nu1=theta[2]
  sigma1=theta[3]
  a2=theta[4]
  nu2=theta[5]
  sigma2=theta[6]
  nug1=theta[7]
  nug2=theta[8]
  R=theta[9]
  nu12 = (nu1 + nu2) /2
  a12 = (a1 + a2)/2
  ######## Putting hard constraints on the parameters #############
  if( a1<=0 | nu1<=0 | sigma1<=0 | a2<=0 | nu2<=0 | sigma2<=0 |nug1<0|nug2<0)
  {
    return(list(mlv=Inf))
  }
  else
  {
    n=nrow(dist.mat)
    constant = (sqrt(sigma1*sigma2)*R*(gamma(nu12)/sqrt(gamma(nu1)*gamma(nu2)))) /
      ((a12^nu12)/sqrt(a1^nu1 * a2^nu2))
    C11=sigma1*matern(dist.mat,sqrt(a1),nu1)
    C22=sigma2*matern(dist.mat,sqrt(a2),nu2)
    C12=constant*matern(dist.mat,sqrt(a12),nu12)
    
    NUG1=diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    NUG2=diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    
    
    ############## Inverting C11 ##########
    C=rbind(cbind(C11+NUG1,C12),cbind(t(C12),C22+NUG2))
    
    # print("I am here 1")
    
    ## Inverse calculation of C
    H = eigen(C)
    e_vec = H$vectors
    D = H$values
    D = 1/D
    for(i in 1:length(D)){
      if(D[i] > 1e+05) D[i] = 0
    }
    D.inv = diag(D)
    Cinv = e_vec%*%D.inv%*%t(e_vec)
    #################
    
    logD=determinant(C)$modulus
    # print("I am here 2")
    
    ## Inverse calculation of V
    H = eigen(t(X)%*%Cinv%*%X)
    e_vec = H$vectors
    D = H$values
    D = 1/D
    for(i in 1:length(D)){
      if(D[i] > 1e+05) D[i] = 0
    }
    D.inv = diag(D)
    V  = e_vec%*%D.inv%*%t(e_vec)
    #################
    
    beta = V%*%t(X)%*%Cinv%*%z
    # print("I am here 3")
    
    nloglikelihood =(0.5 * logD + 0.5 * t(z - X%*%beta) %*% Cinv %*% (z - X%*%beta)+0.5*length(z)*log(2*pi))
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood = 1e+08}
    return(list(mlv=nloglikelihood,full.cov=C))
  }
}

par_Matern.loglikelihood=function(par)
{
  
  return(par_Matern.loglikelihood_allcomponents(p=par,data_train)$mlv)
}

optim_par_Matern.loglikelihood = function(par){
  optim(par=par,
        fn = par_Matern.loglikelihood,
        hessian=FALSE, 
        control=list(trace=4,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=300))
}

par_Matern.pred.summary=function(estim.par, data_train, test_data)
{
  if(dim(data_train)[2] > 4){
    un.grd.train = data_train[,c(6:9)]
    a = data.matrix(data_train[,c(1:5)])
    b = matrix(0, dim(a)[1], dim(a)[2])
    X = rbind(cbind(a,b),cbind(b,a))    
    test.data = test_data[,c(6:9)]
    a = data.matrix(test_data[,c(1:5)])
    b = matrix(0, dim(a)[1], dim(a)[2])
    X0 = rbind(cbind(a,b),cbind(b,a))
  }
  else{
    un.grd.train = data_train[,c(1,2,3,4)]
    X = matrix(c(rep(1,dim(un.grd.train)[1]),rep(0,dim(un.grd.train)[1]),
                 rep(0,dim(un.grd.train)[1]),rep(1,dim(un.grd.train)[1])), 
               nrow = 2*(dim(un.grd.train)[1]),
               ncol = 2, byrow = F)
    test.data = test_data[,c(1,2,3,4)]
    X0 = matrix(c(rep(1,dim(test.data)[1]),rep(0,dim(test.data)[1]),
                  rep(0,dim(test.data)[1]),rep(1,dim(test.data)[1])), 
                nrow = 2*(dim(test.data)[1]),
                ncol = 2, byrow = F)
  }
  full.data.coords = rbind(as.matrix(test.data[,c(1,2)]),as.matrix(un.grd.train[,c(1,2)]))
  dist.mat=rdist(full.data.coords)
  p=estim.par
  print(dim(full.data.coords))
  
  theta=p
  a1=theta[1]
  nu1=theta[2]
  sigma1=theta[3]
  a2=theta[4]
  nu2=theta[5]
  sigma2=theta[6]
  nug1=theta[7]
  nug2=theta[8]
  R=theta[9]
  nu12 = (nu1 + nu2) /2
  a12 = (a1 + a2)/2
  
  n=nrow(dist.mat)
  constant = (sqrt(sigma1*sigma2)*R*(gamma(nu12)/sqrt(gamma(nu1)*gamma(nu2)))) /
    ((a12^nu12)/sqrt(a1^nu1 * a2^nu2))
  C11=sigma1*matern(dist.mat,sqrt(a1),nu1)
  C22=sigma2*matern(dist.mat,sqrt(a2),nu2)
  C12=constant*matern(dist.mat,sqrt(a12),nu12)
  
  NUG1=diag(nug1,nrow = nrow(C11),ncol=ncol(C11))
  NUG2=diag(nug2,nrow = nrow(C22),ncol=ncol(C22))
  
  ############## Inverting C11 ##########
  C=rbind(cbind(C11+NUG1,C12),cbind(t(C12),C22+NUG2))
  print(dim(C))
  test.index = c(1:dim(test.data)[1],
                 (dim(un.grd.train)[1]+dim(test.data)[1]+1):(dim(un.grd.train)[1]+
                                                               dim(test.data)[1]+dim(test.data)[1]))
  
  C.test=C[test.index,test.index]
  C.train=C[-test.index,-test.index]
  C.test.train=C[test.index,-test.index]
  
  print(dim(C.train))
  print(dim(C.test.train))
  
  z=c(un.grd.train$var1 , un.grd.train$var2)
  
  
  
  ## Inverse of C.train
  H = eigen(C.train)
  e_vec = H$vectors
  D = H$values
  D = 1/D
  for(i in 1:length(D)){
    if(D[i] > 1e+05) D[i] = 0
  }
  D.inv = diag(D)
  C.train.inv<-e_vec%*%D.inv%*%t(e_vec)
  ##################
  
  ## Inverse calculation of V
  H = eigen(t(X)%*%C.train.inv%*%X)
  e_vec = H$vectors
  D = H$values
  D = 1/D
  for(i in 1:length(D)){
    if(D[i] > 1e+05) D[i] = 0
  }
  D.inv = diag(D)
  V  = e_vec%*%D.inv%*%t(e_vec)
  #################
  
  beta = V%*%t(X)%*%C.train.inv%*%z
  
  
  prediction=X0%*%beta + C.test.train%*%(C.train.inv)%*%(z - X%*%beta) #Conditional mean of a Multivariate Gaussian
  cond_var_MAT= C.test - C.test.train%*%(C.train.inv)%*%t(C.test.train)
  
  # prediction = prediction + c(rep(mean(un.grd.train$var1),dim(test.data)[1]),
  #                             rep(mean(un.grd.train$var2),dim(test.data)[1]))
  return(list(pred = prediction, conditional_var = cond_var_MAT))
  
  
}


True.model.estimation=function(test_data,train_data)
{
  g = 0.5
  h = 0.5
  tukey_inverse = inverse(function(x){((exp(g*x)-1)/g) * exp(h*(x^2)/2)}, -20, 20)
  ########## Test set #############
  un.grd.train = train_data[,c(1,2,3,4)]
  test.data = test_data[,c(1,2,3,4)]
  for(i in 1:dim(un.grd.train)[1]){
    un.grd.train[i,3] = tukey_inverse(un.grd.train[i,3])
    un.grd.train[i,4] = tukey_inverse(un.grd.train[i,4])
  }
  
  full.data.coords = rbind(as.matrix(test.data[,c(1,2)]),as.matrix(un.grd.train[,c(1,2)]))
  m = as.matrix(dist(full.data.coords))
  
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
  matern_cov2 = constant*matern(m,sqrt(alpha12),nu12)
  matern_cov4 = s22*matern(m,sqrt(alpha22),nu22)
  C = rbind(cbind(matern_cov1,matern_cov2),cbind(t(matern_cov2),matern_cov4))
  test.index = c(1:dim(test.data)[1],(dim(un.grd.train)[1]+1):(dim(un.grd.train)[1]+dim(test.data)[1]))
  
  C.test=C[test.index,test.index]
  C.train=C[-test.index,-test.index]
  C.test.train=C[test.index,-test.index]
  
  # print(dim(C.train))
  # print(dim(C.test.train))
  
  C.inv = solve(C.train)
  prediction=C.test.train%*%C.inv%*%c(un.grd.train$var1,un.grd.train$var2) #Conditional mean of a Multivariate Gaussian
  cond_var_MAT= C.test - C.test.train%*%C.inv%*%t(C.test.train)
  
  return(list(pred = prediction, conditional_var = cond_var_MAT ))
  
  
}


