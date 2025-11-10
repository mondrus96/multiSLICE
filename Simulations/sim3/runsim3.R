# Showing relationship between number of layers l, sample size n
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  iter <- as.numeric(args[1])
  n <- as.numeric(args[2])
  p <- as.numeric(args[2])
  Sest <- as.character(args[3])
  print(paste("iter:", iter, "n:", n, ", p:", p, ", Sest:", Sest))
}

# Load functions
sapply((paste0("../../Core/", list.files("../../Core/", pattern = "\\.R$"))), source)
Rcpp::sourceCpp('../../Core/GemBag-algo.cpp')

# Load parameters if clime or gscad
if(Sest %in% c("clime", "gscad")){
  file <- paste0("glasso_", n, "samples_iter", iter, ".rda")
  load(file)
  rhos <- cv_rho_out; rs <- cv_rank_out
} else{
  rhos <- NULL; rs <- NULL
}

# Auxiliary parameters
simtype <- "rand"
plat <- l <- 2
p <- 100
iters <- ((iter-1)*25+1):(iter*25)

# Run the simulation
outname <- paste0(Sest, "_", n, "samples_iter", iter)
pobs <- makeseq(p, l)
ns <- rep(n, l)
rs <- as.list(rep(plat, length(iters)))
runsim(simtype, pobs, plat, ns, iters, Sest, outname, rhos = rhos, r = rs)