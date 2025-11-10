library(RSpectra)

# Main multislice estimator
multislice <- function(Sigmas, rhos, rank, Sest = "glasso",
                  tol = 1e-3, maxiter = 100){
  # Sigma = the input covariance matrices, in a list format
  # rho = regularization parameter for clime/graphical lasso
  # rank = rank
  # Sest = sparse estimator
  
  if(!Sest %in% c("glasso", "clime", "gscad", "huge_ct", "huge_glasso")){
    stop(paste(Sest, "is not a valid sparse model"))
  }
  
  if(length(rhos) == 1){
    rhos = rep(rhos, length(Sigmas))
  }
  
  # Apply SLICE
  labs <- lapply(Sigmas, nrow)
  labs <- Map(seq, cumsum(c(1, head(labs, -1))), cumsum(labs))
  slices <- mapply(slice, Sigmas, rhos, MoreArgs = list(rank = rank,
                                                        Sest = Sest,
                                                        tol = tol,
                                                        maxiter = maxiter), SIMPLIFY = FALSE) # Apply independent slice models
  Ls <- lapply(slices, "[[", "L"); Ss <- lapply(slices, "[[", "S") # Get Ls and Ss
  Ss <- Map(function(matrix, labels){ # For relabeling Ss with actual labs
    colnames(matrix) <- labels
    rownames(matrix) <- labels
    matrix
  }, Ss, labs)
  
  L <- list2mat(Ls, labs) # Bring together Ls
  L <- matcomp(L, rank, labs)
  
  # Return result
  slices_misc <- lapply(slices, "[[", "misc")
  result <- list(S = Ss, L = L, rhos = rhos, rank = rank,
                 misc = list(converged =
                               sapply(slices_misc, "[[", "converged",
                                        simplify = FALSE),
                             iters = sapply(slices_misc, "[[", "iters", 
                                              simplify = FALSE),
                             deltaS = sapply(slices_misc, "[[", "deltaS",
                                             simplify = FALSE),
                             deltaL = sapply(slices_misc, "[[", "deltaL",
                                             simplify = FALSE)))
  return(result)
}

# Matrix completion step
matcomp <- function(L, rank, labs) {
  # One‐pass truncated SVD per block to build factor H
  p <- nrow(L)
  H <- matrix(0, p, rank)
  for (idx in labs) {
    X  <- L[idx, idx]
    sv <- svds(X, rank)
    H[idx, ] <- sv$u[, 1:rank] %*% diag(sqrt(sv$d[1:rank]), rank, rank)
  }
  
  # Reconstruct low‐rank completion
  return(H %*% t(H))
}