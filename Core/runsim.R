library(MASS)
library(BANS)
library(R.matlab)
library(matrixStats)
library(coglasso)
library(qgraph)

runsim <- function(simtype, pobs, plat = NULL, ns, iters, Sest = "glasso", outname,
                   maxiter = 100, tol = 1e-3, rhos = NULL, rs = NULL, Kmax = 5){
  S_hat_out <- L_hat_out <- S_star_out <- L_star_out <- z_star_out <- cv_rho_out <- cv_rank_out <- vector("list", length(iters))
  # Loop through the iterations of the simulation
  for(i in 1:length(iters)){
    T_start <- Sys.time()
    print(paste0("SIM ITER ", iters[i]))
    set.seed(123*iters[i])
    l <- length(pobs)
    
    labs <- unique(unlist(pobs)) # Labels
    if (simtype == "rand"){
      Lout <- Lrand(length(labs), plat, 1.5)
    } else if (simtype == "cres"){
      Lout <- Lcres(length(labs), 0.1)
      plat <- 2
    } else if (simtype == "spir"){
      Lout <- Lspir(length(labs), 0.05)
      plat <- 2
    }
    
    L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels
    colnames(L_star) <- rownames(L_star) <- labs
    L_stars <- mat2list(L_star, pobs)
    
    Sigma_stars <- Xs <- Sigmas <- S_stars <- vector("list", l)
    for(j in 1:l){
      S_stars[[j]] <- Smat(length(pobs[[j]]), 2, 1.5)
      S_stars[[j]][S_stars[[j]] < 0.01] <- 0 # True sparse component
      
      Sigma_stars[[j]] <- solve(S_stars[[j]] + L_stars[[j]]) # True Sigma
      Xs[[j]] <- mvrnorm(ns[[j]], rep(0, length(pobs[[j]])), 
                         Sigma = Sigma_stars[[j]]) # Finite sample data
      Sigmas[[j]] <- cov(Xs[[j]]) # Sample Sigma 
    }
    
    cvsli <- cv.multislice(Xs, Sest = Sest, 
                           rhos = seq(4, 18, length.out = 5)*
                             sqrt(log(ncol(Xs[[1]])))/nrow(Xs[[1]]), 
                           rs = rs[[i]])
    sli <- multislice(Sigmas, cvsli$rho, cvsli$r, Sest, maxiter=maxiter, tol=tol)
    Ss <- sli$S; L <- sli$L
    
    # Append to all outputs
    if(K == Kmax){
      S_hat_out[(i - K + 1):i] <- Ss; L_hat_out[(i - K + 1):i] <- Ls
      # Reset K and datalist
      K <- 0
      datalist <- vector("list", Kmax)}
    
    cv_rho_out[[i]] <- sli$rho; cv_rank_out[[i]] <- sli$rank # Append CV information as well 
    
    # Save as rda
    save(S_hat_out, L_hat_out, S_star_out, L_star_out, z_star_out, cv_rho_out, cv_rank_out, file = paste0(outname, ".rda"))
    
    T_end <- Sys.time()
    print(T_end - T_start)
  }
}

set.mat.names <- function(matrices, labs) {
  Map(function(matrix, labs) {
    colnames(matrix) <- labs
    rownames(matrix) <- labs
    matrix
  }, matrices, labs)
}