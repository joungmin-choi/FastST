# Simulation for Scenario 1a (proportion estimation) with fully simulated microbiome data
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
nsimu <- 1000

amae.observed.mat <- matrix(NA, nk, nn)
amae.unobserved.mat <- matrix(NA, nk, nn)
amae.major.mat <- matrix(NA, nk, nn)
runt.mat <- matrix(NA, nk, nn)
set.seed(123)
for (i in 1 : nk)
{
  K <- K.vec[i]
  for (j in 1 : nn)
  {
    nmajor <- nmajor.vec[j]
    mae.observed.vec <- rep(NA, nsimu)
    mae.unobserved.vec <- rep(NA, nsimu)
    mae.major.vec <- rep(NA, nsimu)
    if (nmajor <= K)
    {
      kk <- 0
      runt <- system.time({
        while(kk < nsimu)
        {
          res <- gendata(N, K, C, nmajor = nmajor, pmajor = 0.9, lb.amajor = 0.1)
          # x.vec <- as.vector(res$x.vec)
          # yobs.mat <- res$y.mat[-1, ]
          # data.mat <- rbind(x.vec, yobs.mat)
          # write.table(cbind(1 : N, t(data.mat)), file = paste0("./Result/Scen2_data/data_K", K, "_nmajor", nmajor, "_simu", kk, ".txt"), quote = F, sep = "\t", row.names = F, col.names = c("out_id", "SNK", paste0("SRC", 1 : K)))
          suppressWarnings({estres <- glsest(res$x.vec, res$y.mat[-1, ])})
          alpha.est.vec <- estres[[1]]
          if (sum(is.na(alpha.est.vec)) == 0)
          {
            kk = kk + 1
            mae.observed.vec[kk] <- mean(abs(alpha.est.vec[-1] - res$alpha.vec[-1]))
            mae.unobserved.vec[kk] <- abs(alpha.est.vec[1] - res$alpha.vec[1] * mean(res$gamma.mat[1, ]))
            mae.major.vec[kk] <- mean(abs(alpha.est.vec[2 : (nmajor + 1)] - res$alpha.vec[2 : (nmajor + 1)]))
          }
        }
      })
      amae.observed.mat[i, j] <- mean(mae.observed.vec)
      amae.unobserved.mat[i, j] <- mean(mae.unobserved.vec)
      amae.major.mat[i, j] <- mean(mae.major.vec)
      runt.mat[i, j] <- runt[3]
      print(paste0("Number of observed sources: ", K, ", number of major sources: ", nmajor))
      print(paste0("amae.observed: ", round(amae.observed.mat[i, j], 4), ", amae.major: ", round(amae.major.mat[i, j], 4), ", amae.unobserved: ", round(amae.unobserved.mat[i, j], 4)))
      print(paste0("running time: ", runt.mat[i, j]))
    }
  }
}

save(amae.observed.mat, amae.unobserved.mat, amae.major.mat, runt.mat, file = "./Result/Simu_Scen1.RData")
