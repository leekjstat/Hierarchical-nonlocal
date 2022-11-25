
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # setting the workspace to source file location
source("ind_pmom.R")
source("H_ind_pmom.R")
source("SSS.R")
source("H_RSSS.R")
source("result.R")

library(MASS)
library(truncnorm)
library(BoomSpikeSlab)
library(ncvreg)
library(EBglmnet)

n.rep = 10

n = 100 # the number of sample size
p = 100 # the number of variables
t = 3 # the number of nonzero coefficients


# Hyperparameters
lambda1 = 1
lambda2 = 10^{2} * n^{-1/3} * p^{(2 + 0.001)*2/3}

# True coefficients (value)
# beta.val = 1 # setting 1
beta.val = 2 # setting 2

# True coefficients (index)
b = rep(beta.val,t)
Bc = c(b, array(0, p-t)) 

# Covariance matrix for design matrix
rho = 0 # cov.set 1
# rho = 0.3 # cov.set 2

Sigma <- diag(p)
for(j in 1:(p-1)){
  for(k in (j+1):p)
    Sigma[j,k] = Sigma[k,j] = rho^abs(j-k)
}

Selection_list = list()
for(seed.num.i in 1:n.rep){
  Selection_mat = matrix(NA, nrow=6, ncol=9)
  colnames(Selection_mat) <- c("Precision", "Recall", "Sensitivity", "Specificity", "MCC", "TP", "TN", "FP", "FN")
  rownames(Selection_mat) <- c("H-nonlocal", "nonlocal", "SS", "LASSO", "EBlasso", "SCAD")
  Selection_list[[seed.num.i]] = Selection_mat
}
MSPE_mat = matrix(NA, nrow=6, ncol=n.rep)
rownames(MSPE_mat) <- c("H-nonlocal", "nonlocal", "SS", "LASSO", "EBlasso", "SCAD")

# Start simulation study
for(seed.num in 1:n.rep){
  
  # Generating data
  set.seed(seed.num)
  X <- mvrnorm(n, rep(0, p), Sigma=Sigma)
  X <- scale(X)
  y = rlogis(n, location = X %*% Bc)
  y = ifelse(y>0, 1,0) # Logistic response E in the paper
  
  
  # Run SSS for nonlocal prior
  set.seed(seed.num)
  fit_de_SSS = SSS(X, y, ind_fun=pmom_laplace, N=150, C0=1, verbose=FALSE)
  
  # Run RSSS for hierarchical nonlocal prior
  set.seed(seed.num)
  H_fit_de_RSSS = H_RSSS(X, y, lambda1, lambda2, ind_fun=H_pmom_laplace, N=150, C0=1, verbose=FALSE) # use corr(X, y)
  # H_fit_de_RSSS = H_RSSS(X, y, lambda1, lambda2, ind_fun=H_pmom_laplace, N=150, C0=1, verbose=FALSE, method="working") # use working residual instead of y
  
  res_de_SSS = result(fit_de_SSS)
  H_res_de_RSSS = result(H_fit_de_RSSS)
  
  ################ Comparision of variable selection methods ###################################
  X_train = X
  Y_train = y
  
  # Hierarchical nonlocal & nonlocal priors
  nonzerobetaid = (res_de_SSS$hppm)
  H_nonzerobetaid = (H_res_de_RSSS$hppm)
  H_nonlocalgamma <- nonlocalgamma <- rep(0, p)
  nonlocalgamma[nonzerobetaid] = 1
  H_nonlocalgamma[H_nonzerobetaid] = 1
  true_model = c(rep(1, t), rep(0, p - t))
  
  # spike and slab
  dataSS = as.data.frame(cbind(X_train, Y_train))
  model <- logit.spike(dataSS$Y_train ~ . -1, data = as.data.frame(dataSS), niter = 2000, ping = 1000, seed=seed.num)
  nonzerobeta_ss <- which(abs(colMeans(model$beta)) > 0.01)
  ssgamma = rep(0, p)
  ssgamma[nonzerobeta_ss] = 1
  
  # LASSO
  fit_lasso <- cv.ncvreg(X_train, Y_train, family="binomial", penalty="lasso", nfolds=5, seed=seed.num)
  lasso <- as.vector(coef(fit_lasso, lambda=fit_lasso$lambda[fit_lasso$min]))
  lasso <- lasso[-1] # delete intercept
  
  # Empirical Bayesian Lasso (EBlasso)
  set.seed(seed.num)
  ebglmnet <- cv.EBglmnet(as.matrix(X_train), Y_train, family=c("binomial"),prior= c("lassoNEG"), nfolds=5,Epis = FALSE, group = FALSE, verbose = 0)
  # ebglmnet <- EBglmnet(as.matrix(X_train), Y_train, family=c("binomial"),prior= c("lassoNEG"), hyperparameters=c(-0.3, 0.1), Epis = FALSE,group = FALSE, verbose = 0)
  ebfit <- as.vector(ebglmnet$fit[,2])
  ebfitbeta <- as.vector(ebglmnet$fit[,3])
  eb <- rep(0, p)
  eb[ebfit] <- 1
  
  # SCAD
  fit_SCAD <- cv.ncvreg(X_train, Y_train, family="binomial", penalty="SCAD", nfolds=5, seed=seed.num)
  scad <- as.vector(coef(fit_SCAD, lambda=fit_SCAD$lambda[fit_SCAD$min]))
  scad <- scad[-1] # delete intercept
  
  
  # Evaluation variable selection performance
  Selection_list[[seed.num]][1,] = Evaluation(true_model, H_nonlocalgamma) # Hierarchical nonlocal
  Selection_list[[seed.num]][2,] = Evaluation(true_model, nonlocalgamma) # nonlocal
  Selection_list[[seed.num]][3,] = Evaluation(true_model, ssgamma) # spike and slab
  Selection_list[[seed.num]][4,] = Evaluation(true_model, 1*(lasso != 0)) # LASSO
  Selection_list[[seed.num]][5,] = Evaluation(true_model, eb) # EBlasso
  Selection_list[[seed.num]][6,] = Evaluation(true_model, 1*(scad != 0)) # SCAD
  
  
  ################ Calculating MSPE ################################################
  # Testing set
  set.seed(seed.num)
  n_test = 50
  X_test = mvrnorm(n_test,rep(0, p), Sigma=Sigma)
  Y_test = rlogis(n_test, location = X_test %*% Bc)
  Y_test = ifelse(Y_test>0, 1,0) 
  
  # Hierarchical nonlocal
  beta.est.H = H_pmom_laplace_beta(X_train[,H_nonzerobetaid], Y_train, r=1, lambda1, lambda2)
  phi.H.nonlocal = X_test[,H_nonzerobetaid, drop=FALSE] %*% beta.est.H
  H_predfitted = 1/(1 + exp(-phi.H.nonlocal))
  MSPE_H_nonlocal = round(mean((Y_test - H_predfitted)^2), digits = 4)
  
  # nonlocal priors
  beta.est = pmom_laplace_beta(X_train[,nonzerobetaid], Y_train, r=1, tau = p^(2+0.05)/sqrt(n))
  phi.nonlocal = X_test[,nonzerobetaid, drop=FALSE] %*% beta.est
  predfitted = 1/(1 + exp(-phi.nonlocal))
  MSPE_nonlocal = round(mean((Y_test - predfitted)^2), digits = 4)
  
  # spike and slab
  beta.est.ss = as.vector(colMeans(model$beta)[nonzerobeta_ss])
  phi.nonlocal.ss = X_test[,as.vector(nonzerobeta_ss), drop=FALSE] %*% beta.est.ss
  predfitted.ss = 1/(1 + exp(-phi.nonlocal.ss))
  MSPE_ss = round(mean((Y_test - predfitted.ss)^2), digits = 4)
  
  # LASSO
  Y_hat_lasso <- as.vector(predict(fit_lasso, as.matrix(X_test), type="response", lambda=0.08))
  MSPE_glasso <- round(mean((Y_test - Y_hat_lasso)^2), digits = 4)
  
  # EBlasso
  Y_hat_eb <- exp(X_test[,ebfit]*ebfitbeta)/(1+exp(X_test[,ebfit]*ebfitbeta))
  MSPE_elasso <- round(mean((Y_test - Y_hat_eb)^2), digits = 4)
  
  # SCAD
  Y_hat_scad <- as.vector(predict(fit_SCAD, as.matrix(X_test), type="response", lambda=0.058))
  MSPE_scad <- round(mean((Y_test - Y_hat_scad)^2), digits = 4)
  
  
  MSPE_mat[1, seed.num] = MSPE_H_nonlocal # Hierarchical nonlocal
  MSPE_mat[2, seed.num] = MSPE_nonlocal # nonlocal
  MSPE_mat[3, seed.num] = MSPE_ss # spike and slab
  MSPE_mat[4, seed.num] = MSPE_glasso # LASSO
  MSPE_mat[5, seed.num] = MSPE_elasso # EBlasso
  MSPE_mat[6, seed.num] = MSPE_scad # SCAD
  
  cat(seed.num,"th simulation is completed! \n")
}

Selection_mat = matrix(0, nrow=6, ncol=9)
for(seed.num.i in 1:n.rep){
  Selection_mat = Selection_mat + Selection_list[[seed.num.i]]
}
Selection_mat = Selection_mat/n.rep

round(cbind(Selection_mat[,-c(2, 6:9)]), digits = 3)
round(rowMeans(MSPE_mat), digits = 3)


# filename = paste("n",n,"p",p,"t",t,"_betaval",beta.val,"_rho",rho,".RData", sep = "") 
# save(list=ls(), file=filename)

