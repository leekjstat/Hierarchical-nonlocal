
## Laplace approximation for (\beta_k)

H_pmom_laplace <- function(X.ind, y, r, lambda1, lambda2, display = FALSE) {
  n_k <- ncol(X.ind)
  if(length(n_k) == 0) return(-100000)
  
  # modified by KL
  obj_func <- function(beta) {
    inner <- X.ind %*% beta
    obj <- sum(inner * y - log(1 + exp(inner))) - n_k * log(prod(seq(1, 2*r - 1, 2))) - 
      n_k/2 * log(2*pi) - (r*n_k + 0.5*n_k + 0.5) * log(lambda2 + 0.5*sum(beta^2)) +
      2*r*sum(log(abs(beta))) + lgamma(r*n_k + 0.5*n_k + lambda1) + 
      lambda1 * log(lambda2) - lgamma(lambda1)
    return(-obj)
  }
  
  # modified by KL
  obj_grad <- function(beta) {
    inner <- X.ind %*% beta
    grad <- colSums(X.ind * as.vector(y - exp(inner) / (1 + exp(inner)))) + 
      2*r / beta - (r*n_k + 0.5*n_k + 0.5) * beta / (lambda2 + 0.5*sum(beta^2))
    return(-grad)
  }
  
  # modified by KL
  # res <- optim(rep(0.5, ncol(X.ind) + 1), obj_func, obj_grad, method = "BFGS")
  Error.ind = 0
  term.ind = -1
  tryCatch(res <- optim(rnorm(n_k, mean=0, sd=2), obj_func, obj_grad, method = "BFGS"), 
           error = function(e){ Error.ind <<- 1; term.ind <<- Error.ind })
  # if optim produces an error, try different initial value
  while(Error.ind == term.ind){
    cat("Error occurred", Error.ind, "times. \n")
    tryCatch(res <- optim(rnorm(n_k, mean=0, sd=2), obj_func, obj_grad, method = "BFGS"), 
             error = function(e){ Error.ind <<- Error.ind + 1; term.ind <<- Error.ind })
    if(Error.ind >= 10){
      res = list()
      res$par = rep(1000, n_k)
      print("Error!")
      break
    } 
  }
  # res$par
  
  # modified by KL
  get_V <- function(beta) {
    inner <- X.ind %*% beta
    V = matrix(0, n_k, n_k)
    V <- - t(X.ind) %*% (X.ind * as.vector(exp(inner) /((1 + exp(inner))^2) )) - 
      diag(2*r / (beta^2)) - 
      (r*n_k + 0.5*n_k + 0.5) * ( diag(n_k)/(lambda2 + 0.5*sum(beta^2)) - beta%*%t(beta)/(lambda2 + 0.5*sum(beta^2))^2 )
    return(V)
  }
  
  # modified by KL
  get_log_marginal <- function(beta) {
    V <- get_V(beta)
    if(sum(is.na(V)) != 0){
      log_marginal <- -100000
    }else{
      log_marginal <- 0.5*n_k* log(2*pi) - obj_func(beta) - 0.5 * sum(log(abs(eigen(V)$values)))
    }
    return(log_marginal)
  }
  
  if(display == TRUE) print(res$par)
  return(get_log_marginal(res$par))
}


H_pmom_laplace_beta <- function(X.ind, y, r, lambda1, lambda2) {
  n_k <- ncol(X.ind)
  if(length(n_k) == 0) return(-100000)
  
  # modified by KL
  obj_func <- function(beta) {
    inner <- X.ind %*% beta
    obj <- sum(inner * y - log(1 + exp(inner))) - n_k * log(prod(seq(1, 2*r - 1, 2))) - 
      n_k/2 * log(2*pi) - (r*n_k + 0.5*n_k + 0.5) * log(lambda2 + 0.5*sum(beta^2)) +
      2*r*sum(log(abs(beta))) + lgamma(r*n_k + 0.5*n_k + lambda1) + 
      lambda1 * log(lambda2) - lgamma(lambda1)
    return(-obj)
  }
  
  # modified by KL
  obj_grad <- function(beta) {
    inner <- X.ind %*% beta
    grad <- colSums(X.ind * as.vector(y - exp(inner) / (1 + exp(inner)))) + 
      2*r / beta - (r*n_k + 0.5*n_k + 0.5) * beta / (lambda2 + 0.5*sum(beta^2))
    return(-grad)
  }
  
  # modified by KL
  # res <- optim(rep(0.5, ncol(X.ind) + 1), obj_func, obj_grad, method = "BFGS")
  Error.ind = 0
  term.ind = -1
  tryCatch(res <- optim(rnorm(n_k, mean=0, sd=2), obj_func, obj_grad, method = "BFGS"), 
           error = function(e){ Error.ind <<- 1; term.ind <<- Error.ind })
  # if optim produces an error, try different initial value
  while(Error.ind == term.ind){
    # print("Error occurred")
    tryCatch(res <- optim(rnorm(n_k, mean=0, sd=2), obj_func, obj_grad, method = "BFGS"), 
             error = function(e){ Error.ind <<- Error.ind + 1; term.ind <<- Error.ind })
    if(Error.ind >= 10){
      res = list()
      res$par = rep(1000, n_k)
      print("Error!")
      break
    } 
  }
  
  return(res$par)
}

