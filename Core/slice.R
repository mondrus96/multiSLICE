library(clime)
library(glasso)
library(huge)
library(RSpectra)
library(Matrix)

# Main slice estimator
slice <- function(Sigma, rho, rank, Sest = "glasso",
                  tol = 1e-3, maxiter = 100){
  # Sigma = the input covariance matrix
  # rho = regularization parameter for clime/graphical lasso
  # rank = rank
  # Sest = sparse estimator
  
  if(!Sest %in% c("glasso", "clime", "gscad", "huge_ct", "huge_glasso")){
    stop(paste(Sest, "is not a valid sparse model"))
  }
  
  p <- ncol(Sigma) # Make Sigma PD
  Sigma <- makePD(Sigma)
  invSigma <- Matrix::chol2inv(Matrix::chol(Sigma))
  
  L <- 0 # zero initialization L
  E <- invSigma - L # Expectation
  S <- 0 # Empty S
  
  deltaS <- deltaL <- deltalogL <- c() 
  for(i in 1:maxiter){
    if(!isPD(E)){
      E <- makePD(E) # Make expectation PD
    }
    
    # Sparse step
    Sold <- S
    if(Sest == "glasso"){
      S <- glasso(Matrix::chol2inv(Matrix::chol(E)), rho, 
                  thr = tol, maxit = maxiter)$wi
    } else if(Sest == "clime"){
      S <- clime(Matrix::chol2inv(Matrix::chol(E)), rho, 
                 sigma = TRUE, linsolver = "simplex")$Omegalist[[1]]
      S[abs(S) < tol] <- 0
    } else if(Sest == "gscad"){
      S <- gscad(Matrix::chol2inv(Matrix::chol(E)), rho)
    } else if(Sest == "huge_glasso"){
      S <- huge(Matrix::chol2inv(Matrix::chol(E)), rho, 
                method = "glasso", verbose = FALSE)$icov[[1]]
    }
    S <- (S + t(S))/2
    
    # Latent step
    Lold <- L
    tsvdL <- svds(invSigma - S, rank)
    L <- tsvdL$v %*% diag(tsvdL$d) %*% t(tsvdL$v)
    L <- (L + t(L))/2
    
    # New expectation
    E <- invSigma - L
    E <- (E + t(E))/2
    
    deltaS <- c(deltaS, sqrt(sum((S - Sold)^2))); deltaL <- c(deltaL, sqrt(sum((L - Lold)^2))) # Convergence check
    if(i == 1){
      deltalogL <- c(deltalogL, suppressWarnings(abs(logL(Sigma, S + L) - logL(Sigma, Diagonal(p, x = tol)))))
    } else{
      deltalogL <- c(deltalogL, suppressWarnings(abs(logL(Sigma, S + L) - logL(Sigma, Sold + Lold)))) 
    }
    if((deltaS[i] < tol) && (deltaL[i] < tol) || 
       ifelse(is.na(deltalogL[i] < tol), FALSE, deltalogL[i] < tol)){
      break
    }
  }
  return(list(S = S, L = L, rho = rho, rank = rank, 
              misc = list(converged = ifelse(i != maxiter, TRUE, FALSE), iters = i,
              deltaS = deltaS, deltaL = deltaL, deltalogL = deltalogL)))
}