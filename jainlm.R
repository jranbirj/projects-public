jainlm <- function(Y, Xk, mtype='ls', ridged=0){
  n = nrow(Xk)
  p = ncol(Xk)
  X = cbind(1, Xk)
  
  if (length(Y) != n){
    stop("Length of Y must match number of rows in Xk") 
  } 
  if (!mtype %in% c('ls', 'ridge', 'lasso')) {
    stop("mtype must be one of: 'ls', 'ridge', 'lasso'")
  }
  if (p == 0){
    break
  }
  if(is.null(p)){
    p = 1
  }

  if(mtype == 'ls'){
    betahat = solve(t(X)%*%X)%*%t(X)%*%Y
    
  } else if(mtype == 'lasso'){
    library(glmnet)
    lamideal = cv.glmnet(Xk, Y, alpha=1)$lambda.min
    betahat = coef(glmnet(Xk, Y, alpha = 1, lambda = lamideal))
    
  } else if(mtype == 'ridge'){
    betahat = solve(t(X)%*%X + ridged * diag(p + 1))%*%t(X)%*%Y
  
  } else {
    print("Invalid Model Type")
    break
  }
  
  betahat = as.vector(betahat)
  hatmatrix = X%*%solve(t(X)%*%X)%*%t(X)
  
  #model assessment metrics
  yhat = X%*%betahat
  resid = Y - yhat
  
  SSE = sum(resid^2)
  MSE = SSE/(n - p - 1)
  
  SST = var(Y) * (n - 1)
  MST = SST / (n - 1)

  SSM = SST - SSE
  MSM = SSM / p 
  
  RMSE = sqrt(MSE)
  
  Fstat = MSM / MSE
  pval = pf(Fstat, p, n-p-1, lower.tail=F)
  
  r2 = 1 - SSE/SST
  r2adj = 1 - MSE/MST
  
  betahat_std = betahat[-1] * apply(Xk, 2, sd) / sd(Y)
  betahat_star = betahat[-1] * apply(Xk, 2, sd)

  leverage = diag(hatmatrix)
  sresid = (resid - 0) / (RMSE * sqrt(1 - leverage))
  
  VIF = numeric(p)
  MCI = numeric(p)
  # AV = matrix(0, p, p)
  slopes = c()

  for(i in 1:p) {
    X_i = Xk[,i]
    X_others = Xk[,-i]
    X_others = cbind(1, X_others)
    
    beta_others = solve(t(X_others) %*% X_others) %*% t(X_others) %*% X_i
    fitted_i = X_others %*% beta_others
    residuals_i = X_i - fitted_i
    SSE_i = sum(residuals_i^2)
    SST_i = sum((X_i - mean(X_i))^2)
    r2_i = 1 - SSE_i/SST_i
    
    VIF[i] = 1/(1-r2_i)
    MCI[i] = sqrt(VIF[i])
    
    X_kc = X[, -c(i+1)]
    X_k = X[, c(i+1)]
    beta_Xk = solve(t(X_kc) %*% X_kc) %*% t(X_kc) %*% X_k
    X_fitted = X_kc %*% beta_Xk
    res_x = X_k - X_fitted
    SSEk = sum(res_x^2)
    MSEk = SSEk / (n - p - 2)
    Hk = X_kc %*% solve(t(X_kc) %*% X_kc) %*% t(X_kc)
    leverages_k = diag(Hk)
    std_resid_Xk = res_x / (sqrt(MSEk) * sqrt(1 - leverages_k))
    
    beta_Yk = solve(t(X_kc) %*% X_kc) %*% t(X_kc) %*% Y
    Y_fitted = X_kc %*% beta_Yk
    res_y = Y - Y_fitted
    SSE_Y = sum(res_y^2)
    MSE_Y = SSE_Y / (n - p - 2)
    std_resid_Y = res_y / (sqrt(MSE_Y) * sqrt(1 - leverages_k))

    beta_av = solve(t(std_resid_Xk)%*%std_resid_Xk)%*%t(std_resid_Xk)%*%std_resid_Y
    slopes = c(slopes, beta_av)
  }
  
  R = cor(cbind(Y, Xk))
  det_X = sqrt(det(t(Xk) %*% Xk))
  
  results = list('y_hat'      = yhat, 
                 'bhat'       = betahat, 
                 'bhatstd'    = betahat_std,
                 'bhatstar'   = betahat_star, 
                 'sse'        = SSE, 
                 'rmse'       = RMSE,
                 'ssm'        = SSM,
                 'sst'        = SST, 
                 'mse'        = MSE, 
                 'fstat'      = Fstat, 
                 'pval'       = pval, 
                 'r2'         = r2, 
                 'r2adj'      = r2adj, 
                 'hatmatrix'  = hatmatrix, 
                 'sresid'     = sresid, 
                 'leverage'   = leverage,
                 'vif'        = VIF, 
                 'mci'        = MCI, 
                 'AV'         = slopes,
                 'det_X'      = det_X,
                 'R'          = R)
  return(results)
}