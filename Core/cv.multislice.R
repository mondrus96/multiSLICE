library(pracma)

cv.multislice = function(Xs, folds = 3, rhos = logseq(1e-5, 0.1, 5), rs = c(2:6), Sest = "glasso"){
  # Xs = input data matrices
  # k = number of folds to perform for CV
  # rhos = list of rhos values to try
  # rs = list of rs to try
  # sest = sparse estimator - either glasso, clime, or gscad
  
  ns <- lapply(Xs, nrow) # number of samples
  cvmat <- matrix(NA, length(rs), length(rhos)); rownames(cvmat) <- rs; colnames(cvmat) <- rhos
  
  # Go over grid of rhos and rs
  for(i in 1:length(rs)){
    print(paste0("rank: ", rs[i]))
    for(j in 1:length(rhos)){
      print(paste0("rho: ", rhos[j]))
      
      inds <- mapply(sample, replicate(length(Xs), 1:folds, simplify = FALSE), ns, 
                     MoreArgs = list(replace = TRUE), SIMPLIFY = FALSE) # Define indices
      mulogL <- c()
      for(k in 1:folds){
        train <- mapply(function(x, ind) x[ind != k,], Xs, inds, SIMPLIFY = FALSE) # List of data for train and test
        test <- mapply(function(x, ind) x[ind == k,], Xs, inds, SIMPLIFY = FALSE)
        
        train <- lapply(train, cov); test <- lapply(test, cov) # Define covariance
        
        out <- multislice(train, rhos[j], rs[i], Sest) # Run method
        Ss <- out$S; L <- out$L
        
        labs <- lapply(Xs, colnames) # Labels from each dimension
        Ls <- mat2list(L, labs)

        likl <- mapply(logL, test, Map("+", Ss, Ls)) # Append to mulogL
        mulogL <- c(mulogL, sum(likl))
      }
      cvmat[i, j] <- mean(mulogL)
    }
  }
  best <- which(cvmat == max(cvmat, na.rm=TRUE), arr.ind = TRUE)
  if(nrow(best) > 2){
    best <- best[1,]
  }
  
  return(list(cvmat = cvmat, maxlogL = max(cvmat), 
              rho = rhos[best[2]], r = rs[best[1]]))
}