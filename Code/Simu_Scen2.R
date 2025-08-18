# Simulation for Scenario 2 with semi-synthetic data generated from the Knights et al. dataset
suppressMessages(suppressWarnings(library(gtools)))
suppressMessages(suppressWarnings(library(npreg)))

setwd("C:/Users/Xiaowei/Documents/Work/FastST") # make sure the FastST folder is working directory
source("./Code/funST.R")

gensemidata <- function(file, K, nmajor = 2, pmajor = 0.9, lb.amajor = 0.1)
{
  rawdata <- as.matrix(read.table(filename, sep = "\t"))
  tmpdata <- rawdata[, -1]
  srcdata <- tmpdata[, (305 - 180 + 1) : 305]
  npos.vec <- colSums(srcdata > 0)
  npos.srt <- sort(npos.vec, decreasing = T, index.return = T)
  
  taxaix.vec <- NULL # index of taxa
  for (i in 1 : (K + 1))
  {
    taxaix.vec <- union(taxaix.vec, which(srcdata[, npos.srt$ix[i]] > 0))
  }
  sampleix.vec <- npos.srt$ix[1 : (K + 1)] # index of less sparse sources
  N <- length(taxaix.vec)
  y.mat <- t(srcdata[taxaix.vec, sampleix.vec])
  gamma.mat <- y.mat / matrix(rep(rowSums(y.mat), N), ncol = N)
  amajor.vec <- rdirichlet(1, rep(1, nmajor)) * pmajor
  while(min(amajor.vec) < lb.amajor)
  {
    amajor.vec <- rdirichlet(1, rep(1, nmajor)) * pmajor
  }
  nminor <- K + 1 - nmajor
  pminor <- 1 - pmajor
  aminor.vec <- rdirichlet(1, rep(1, nminor)) * pminor
  alpha.vec <- c(amajor.vec, aminor.vec)
  alpha.vec <- c(alpha.vec[nmajor + 1], alpha.vec[-(nmajor + 1)]) # switch order so that the 1st element is for the unobserved minor source
  beta.vec <- alpha.vec %*% gamma.mat
  C <- round(mean(rowSums(y.mat)))
  # C <- round(mean(rowSums(y.mat)) / 100)
  x.vec <- rmultinom(1, C, beta.vec)
  return(list(x.vec = x.vec, beta.vec = beta.vec, gamma.mat = gamma.mat, alpha.vec = alpha.vec, y.mat = y.mat))
}

# Set parameters
K.vec <- c(2, 5, 10, 50, 100); nk <- length(K.vec)
nmajor.vec <- c(2, 5); nn <- length(nmajor.vec)
nsimu <- 100

filename <- "./Real/otus.txt"
pmajor <- 0.9; lb.amajor <- 0.1
amae.observed.mat <- matrix(NA, nk, nn)
amae.unobserved.mat <- matrix(NA, nk, nn)
amae.major.mat <- matrix(NA, nk, nn)
acon.major.mat <- matrix(NA, nk, nn)
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
    con.major.vec <- rep(NA, nsimu)
    if (nmajor <= K)
    {
      kk <- 0
      runt <- system.time({
        while(kk < nsimu)
        {
          res <- gensemidata(filename, K, nmajor, pmajor, lb.amajor)
          x.vec <- res$x.vec
          beta.vec <- res$beta.vec
          gamma.mat <- res$gamma.mat
          alpha.vec <- res$alpha.vec
          y.mat <- res$y.mat
          yobs.mat <- y.mat[-1, ]
          suppressWarnings({estres <- glsest(x.vec, yobs.mat)})
          alpha.est.vec <- estres[[1]]
          alpha.se.vec <- estres[[2]]
          if (sum(is.na(alpha.est.vec)) == 0)
          {
            kk = kk + 1
            mae.observed.vec[kk] <- mean(abs(alpha.est.vec[-1] - alpha.vec[-1]))
            mae.unobserved.vec[kk] <- abs(alpha.est.vec[1] - alpha.vec[1] * mean(gamma.mat[1, ]))
            mae.major.vec[kk] <- mean(abs(alpha.est.vec[2 : (nmajor + 1)] - alpha.vec[2 : (nmajor + 1)]))
            alpha.srt <- sort(alpha.vec, decreasing = T, index.return = T)
            alpha.est.srt <- sort(alpha.est.vec, decreasing = T, index.return = T)
            con.major.vec[kk] <- length(intersect(alpha.srt$ix[1 : nmajor], alpha.est.srt$ix[1 : nmajor])) / nmajor
          }
        }
      })
      amae.observed.mat[i, j] <- mean(mae.observed.vec)
      amae.unobserved.mat[i, j] <- mean(mae.unobserved.vec)
      amae.major.mat[i, j] <- mean(mae.major.vec)
      acon.major.mat[i, j] <- mean(con.major.vec)
      runt.mat[i, j] <- runt[3]
      print(paste0("Number of known sources: ", K, ", number of major sources: ", nmajor))
      # print(paste0("amae.observed: ", amae.observed.mat[i, j], ", amae.major: ", amae.major.mat[i, j], ", amae.unobserved: ", amae.unobserved.mat[i, j]))
      print(paste0("amae.observed: ", round(amae.observed.mat[i, j], 4), ", amae.major: ", round(amae.major.mat[i, j], 4), ", amae.unobserved: ", round(amae.unobserved.mat[i, j], 4)))
      print(paste0("running time: ", runt.mat[i, j]))
    }
  }
}

save(amae.observed.mat, amae.unobserved.mat, amae.major.mat, acon.major.mat, runt.mat, file = "./Result/Simu_Scen2.RData")