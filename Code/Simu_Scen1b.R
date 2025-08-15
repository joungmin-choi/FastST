# Simulation for Scenario 1b (directionality inference) with fully simulated microbiome data
suppressMessages(suppressWarnings(library(gtools)))
suppressMessages(suppressWarnings(library(npreg)))

setwd("./") # make sure the FastST folder is working directory
source("./Code/funST.R")

res_dirname <- "./Result"
if (!dir.exists(res_dirname)){
 dir.create(res_dirname)
}

# Set parameters
N <- 500; C <- 1e5
K.vec <- c(2, 5, 10, 50, 100); nk <- length(K.vec)
nmajor.vec <- c(2, 5); nn <- length(nmajor.vec)
nsimu <- 100 #1000

accr.mat <- matrix(NA, nk, nn)
runt.mat <- matrix(NA, nk, nn)
set.seed(123)
for (i in 1 : nk)
{
  K <- K.vec[i]
  for (j in 1 : nn)
  {
    nmajor <- nmajor.vec[j]
    sinkid.est.vec <- rep(NA, nsimu)
    if (nmajor <= K)
    {
      kk <- 0
      runt <- system.time({
        while(kk < nsimu)
        {
          res <- gendata(N, K, C, nmajor = nmajor, pmajor = 0.9, lb.amajor = 0.1)
          x.mat <- as.vector(res$x.vec)
          gammaSS.mat <- rbind(res$beta.vec, res$gamma.mat[-1, ])
          for (l in 2 : (K + 1))
          {
            x.mat <- rbind(x.mat, as.vector(rmultinom(1, C, res$gamma.mat[l, ])))
          }
          kk <- kk + 1
          suppressWarnings({sinkid.est.vec[kk] <- DIRinfer(x.mat, gammaSS.mat)})
        }
      })
      accr.mat[i, j] <- mean(sinkid.est.vec == 1)
      runt.mat[i, j] <- runt[3]
      print(paste0("Number of observed sources: ", K, ", number of major sources: ", nmajor))
      print(paste0("Accuracy rate: ", accr.mat[i, j]))
      print(paste0("running time: ", runt.mat[i, j]))
    }
  }
}

save(accr.mat, runt.mat, file = "./Result/Simu_Scen1b.RData")
