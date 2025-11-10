# Showing relationship between amount of overlap o, sample size n
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  iter <- as.numeric(args[1])
  n <- as.numeric(args[2])
  plat <- as.numeric(args[3])
  Sest <- as.character(args[4])
  print(paste("iter:", iter, "n:", n, ", plat:", plat, ", Sest:", Sest))
}

# Load functions
sapply((paste0("../../Core/", list.files("../../Core/", pattern = "\\.R$"))), source)
Rcpp::sourceCpp('../../Core/GemBag-algo.cpp')

# Load parameters if clime or gscad
if(Sest %in% c("clime", "gscad")){
  file <- paste0("glasso_", l, "layers_", n, "samples_iter", iter, ".rda")
  load(file)
  rhos <- cv_rho_out; rs <- cv_rank_out
} else{
  rhos <- NULL; rs <- NULL
}

# Auxiliary parameters
simtype <- "rand"
p <- 100
l <- 2
iters <- ((iter-1)*25+1):(iter*25)

# Run the simulation
outname <- paste0(Sest, "_", plat, "plat_", n, "samples_iter", iter)
pobs <- makeseq(p, l)
ns <- rep(n, l)
runsim(simtype, pobs, plat, ns, iters, Sest, outname, rhos = rhos, r = rs)