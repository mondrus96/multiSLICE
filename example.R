# ---- Setup ----
# Required packages
required_packages <- c("pracma", 
                       "glasso", 
                       "RSpectra",
                       "clime",
                       "huge",
                       "Matrix")
installed <- required_packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(required_packages[!installed], repos = "https://cloud.r-project.org/")
}
lapply(required_packages, library, character.only = TRUE)

# Core files
core_files <- list.files("Core", pattern = "\\.R$", full.names = TRUE)
lapply(core_files, source)

# ---- Simulate data ----
set.seed(123)
plat <- 2
p <- 50
l <- 2
n <- 500
pobs <- makeseq(p, l)

# Generate L matrix
labs <- unique(unlist(pobs)) # Labels
Lout <- Lrand(length(labs), plat, 1.5)
ns <- rep(n, l)
L_star <- Lout$L + 0.1; z_star <- Lout$z # True latent component; True cluster labels
colnames(L_star) <- rownames(L_star) <- labs
L_stars <- mat2list(L_star, pobs)

# Generate Sigmas for each measurement modality
Sigma_stars <- Xs <- Sigmas <- S_stars <- vector("list", l)
for(j in 1:l){
  S_stars[[j]] <- Smat(length(pobs[[j]]), 2, 1.5)
  S_stars[[j]][S_stars[[j]] < 0.01] <- 0 # True sparse component
  
  Sigma_stars[[j]] <- solve(S_stars[[j]] + L_stars[[j]]) # True Sigma
  Xs[[j]] <- mvrnorm(ns[[j]], rep(0, length(pobs[[j]])), 
                     Sigma = Sigma_stars[[j]]) # Finite sample data
  Sigmas[[j]] <- cov(Xs[[j]]) # Sample Sigma 
}

# ---- Estimate via multiSLICE ----
cvout <- cv.multislice(Xs) # hyperparameter selection via CV
out <- multislice(Sigmas, cvout$rho, cvout$r) # estimation

# ---- Plot results ----
# Plot L star (full matrix, including observed and unobserved)
heatmap(L_star, Rowv = NA, Colv = NA)

# Plot observed portions of L star (ie., layerwise L stars)
L_obs <- list2mat(L_stars, pobs)
heatmap(L_obs, Rowv = NA, Colv = NA)

# Plot recovered L hat (ie., the full multilayer network)
heatmap(out$L, Rowv = NA, Colv = NA)

# Plot recovered sparse matrices (ie., layerwise S hats)
Ss <- list2mat(out$S, pobs)
heatmap(1*(Ss != 0), Rowv = NA, Colv = NA)