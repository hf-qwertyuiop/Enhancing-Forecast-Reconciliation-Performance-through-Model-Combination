setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(forecast)
library(foreach)
library(dplyr)


##############
#compute C matrix, which is used to map base forecasts to forecasts for combination
#parm1: S, summing matrix
compute_C <- function(S) {
  #constract Z matrix
  m=dim(S)[2]
  n=dim(S)[1]
  num.of.agg=n-m
  
  C=matrix(0,nrow = sum(S),ncol = n)
  d=1
  for (i in 1:m) {
    #extract non-zeros rows
    rows_with_one <- S[S[, i] == 1, ]
    row_indices <- which(S[, i] == 1)
    num=length(row_indices)
    for (j in 1:num) {
      C[d,row_indices[j]]=1
      col_inds=which(rows_with_one[j,] != 0)
      #find column indicies excluding the target column 
      col_inds_1=col_inds[col_inds!=i]+num.of.agg
      #
      C[d,col_inds_1]=-1
      d=d+1
    }
    
  }
  return(C)
}

library(MASS)
inverse_matrix <- function(matrix) {
  tryCatch({
    inv_matrix <- solve(matrix)
    inv_matrix  
  }, error = function(e) {
    message("Ordinary inverse failed, using generalized inverse instead.")
    tryCatch({
      inv_matrix <- ginv(matrix)
      inv_matrix 
    }, error = function(e) {
      message("Generalized inverse also failed, returning identity matrix instead.")
      inv_matrix <- diag(nrow(matrix))
      inv_matrix  
    })
  })
}
library(quadprog)

optimize_weights <- function(Sigma) {
  N <- nrow(Sigma)
  
  # Define the objective function: Dmat is 2 * Sigma because we need to solve 1/2 * w' * Dmat * w
  Dmat <- 2 * Sigma
  
  # dvec is the linear term in the objective function; here it is a zero vector
  dvec <- rep(0, N)
  
  # Equality constraint: sum(w) = 1
  A_eq <- matrix(1, nrow = 1, ncol = N)
  b_eq <- 1
  
  # Inequality constraints: w >= 0
  A_ge <- diag(N)
  b_ge <- rep(0, N)
  
  # Combine all constraints
  Amat <- t(rbind(A_eq, A_ge))
  bvec <- c(b_eq, b_ge)
  
  # Number of equality constraints
  meq <- 1
  
  # Solve the quadratic programming problem
  solution <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq)
  
  # Return the optimal weights
  return(matrix(solution$solution,nrow = N,ncol = 1))
}
#compute optimal weights given variance-covariance matrix
compute_weights<-function(k,var_cov){
  inv_matrix=inverse_matrix(var_cov)
  print('dim')
  print(var_cov)
  e=matrix(1,nrow = k,ncol=1)
  w=(inv_matrix%*%e)/(t(e)%*%inv_matrix%*%e)[1,1]
  return(w)
}
#compute variance-covariance matrix of errors
est_var_cov<-function(X,method){
  if (method=='sample-var-cov'){
    S <- cov(X)
    return(S)
  }else if(method=='diag'){
    #
    variances <- apply(X, 2, var)
    # #return(diag(dim(S)[1]))
    return(diag(variances))
  }
}

#compute combination weight matrix Phi, shape(m,r)
#parm1: S, summing matrix
#parm2: var_cov_y_c, variance-covariance matrix of errors of candidate forecasts in combination,shape(r,r)
compute_Phi_lessdf<-function(S,var_cov_y_c){
  #model combination
  g=1
  s=1
  m=dim(S)[2]
  #combination weight matrix
  Phi=matrix(0,nrow = m,ncol = sum(S))
  for (i in 1:m) {
    print(i)
    #k:number of forecasts for combination
    k=sum(S[,i])
    weights=compute_weights(k,var_cov_y_c[g:(g+k-1),g:(g+k-1)])
    g=g+k
    Phi[i,s:(s+k-1)]=weights[,1]
    s=s+k
  }
  return(Phi)
}

#compute forecasts and residuls of all nodes
#param1: aggts_all, observations of all time series, shape(T,n)
compute_forecasts<-function(aggts_all,method){
  library(forecast)
  T=dim(aggts_all)[1]
  base_f=matrix(0,nrow = n,ncol = h)
  residual=matrix(0,nrow = n,ncol=T)
  for (i in 1:n) {
    if (method=='arima'){
      fit=auto.arima(aggts_all[,i])
    }else if (method=='ets'){
      fit=ets(aggts_all[,i])
    }
    
    residual[i,]=fit$residuals
    base_f[i,]=forecast(fit,h=2)$mean
  }
  return(list(base_f=base_f,residual=residual))
}

#forecasting performance evaluation metric
mase_cal <- function(insample, outsample, forecasts) {
  stopifnot(stats::is.ts(insample))
  #Used to estimate MASE
  frq <- stats::frequency(insample)
  forecastsNaiveSD <- rep(NA,frq)
  for (j in (frq+1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
  }
  masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)
  
  outsample <- as.numeric(outsample) ; forecasts <- as.numeric(forecasts)
  mase <- (abs(outsample-forecasts))/masep
  return(mase)
}
