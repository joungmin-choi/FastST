# library(gtools) # for function rdirichlet
# library(npreg) # for function psolve

gendata <- function(N, K, C, nmajor = 2, pmajor = 0.9, lb.amajor = 0.1, C.vec = NULL)
{
  # generate sink data by using mixed-proportion model, for simulation purpose
  # note, the unknown source is always not in major sources
  # Input: 
  #         N: scalar, number of taxa groups in source/sink
  #         K: scalar, number of KNOWN sources
  #         C: scalar, total taxa count in sink
  #         nmajor: scalar, number of major sources
  #         pmajor: scalar, sum of major proportions
  #         lb.amajor: scalar, lower bound of major proportions
  #         C.vec: vector, total taxa counts in K sources
  # Output: 
  #         a list with the following components:
  #           x.vec: vector of length N, sink data (check: sum up to C)
  #           beta.vec: vector of length N, abundance of sink (check: sum up to 1)
  #           gamma.mat: matrix of size (K+1) by N, abundance of KNOWN+UNKNOWN sources (note: sum up each row to 1)
  #                      note, the 1st row is reserved for the unobserved source
  #           alpha.vec: vector of length K+1, proportion of KNOWN+UNKNOWN sources (note: sum up to 1)
  #                      note, unobserved source is always put in the 1st element
  #           y.mat: matrix of size (K+1) by N, data of KNOWN+UNKNOWN sources
  #                      note, the 1st row is reserved for the unobserved source
  
  if (nmajor > K)
  {
    stop("Too many major sources! Need to increase K or decrease nmajor!")
  }
  if (nmajor * lb.amajor > pmajor)
  {
    stop("Bad configuration for major sources! Need to adjust nmajor, lb.amajor, and pmajor!")
  }
  gamma.mat <- rdirichlet(K + 1, rep(1, N))
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
  x.vec <- rmultinom(1, C, beta.vec)
  if (is.null(C.vec))
  {
    C.vec <- rep(C, K + 1)
  }
  y.mat <- round(gamma.mat * matrix(rep(C.vec, N), nrow = K + 1))
  return(list(x.vec = x.vec, beta.vec = beta.vec, gamma.mat = gamma.mat, alpha.vec = alpha.vec, y.mat = y.mat))
}

gendata_s <- function(N, K, C, nmajor = 2, pmajor = 0.9, lb.amajor = 0.1, C.vec = NULL)
{
  # generate sink data by using mixed-proportion model, for simulation purpose (major sources have identical taxa distributions)
  # note, the unknown source is always not in major sources
  # Input: 
  #         N: scalar, number of taxa groups in source/sink
  #         K: scalar, number of KNOWN sources
  #         C: scalar, total taxa count in sink
  #         nmajor: scalar, number of major sources
  #         pmajor: scalar, sum of major proportions
  #         lb.amajor: scalar, lower bound of major proportions
  #         C.vec: vector, total taxa counts in K sources
  # Output: 
  #         a list with the following components:
  #           x.vec: vector of length N, sink data (check: sum up to C)
  #           beta.vec: vector of length N, abundance of sink (check: sum up to 1)
  #           gamma.mat: matrix of size (K+1) by N, abundance of KNOWN+UNKNOWN sources (note: sum up each row to 1)
  #                      note, the 1st row is reserved for the unobserved source
  #           alpha.vec: vector of length K+1, proportion of KNOWN+UNKNOWN sources (note: sum up to 1)
  #                      note, unobserved source is always put in the 1st element
  #           y.mat: matrix of size (K+1) by N, data of KNOWN+UNKNOWN sources
  #                      note, the 1st row is reserved for the unobserved source
  
  if (nmajor > K)
  {
    stop("Too many major sources! Need to increase K or decrease nmajor!")
  }
  if (nmajor * lb.amajor > pmajor)
  {
    stop("Bad configuration for major sources! Need to adjust nmajor, lb.amajor, and pmajor!")
  }
  gamma1.vec <- rdirichlet(1, rep(1, N))
  gamma1.mat <- matrix(rep(gamma1.vec, each = nmajor), nrow = nmajor)
  amajor.vec <- rdirichlet(1, rep(1, nmajor)) * pmajor
  while(min(amajor.vec) < lb.amajor)
  {
    amajor.vec <- rdirichlet(1, rep(1, nmajor)) * pmajor
  }
  nminor <- K + 1 - nmajor
  gamma2.mat <- rdirichlet(nminor, rep(1, N))
  pminor <- 1 - pmajor
  aminor.vec <- rdirichlet(1, rep(1, nminor)) * pminor
  alpha.vec <- c(amajor.vec, aminor.vec)
  alpha.vec <- c(alpha.vec[nmajor + 1], alpha.vec[-(nmajor + 1)]) # switch order so that the 1st element is for the unobserved minor source
  gamma.mat <- rbind(gamma1.mat, gamma2.mat)
  gamma.mat <- rbind(gamma.mat[nmajor + 1, ], gamma.mat[-(nmajor + 1), ]) # switch order so that the 1st element is for the unobserved minor source
  beta.vec <- alpha.vec %*% gamma.mat
  x.vec <- rmultinom(1, C, beta.vec)
  if (is.null(C.vec))
  {
    C.vec <- rep(C, K + 1)
  }
  y.mat <- round(gamma.mat * matrix(rep(C.vec, N), nrow = K + 1))
  return(list(x.vec = x.vec, beta.vec = beta.vec, gamma.mat = gamma.mat, alpha.vec = alpha.vec, y.mat = y.mat))
}

gendata_ss <- function(N, K, C, nmajor = 2, pmajor = 0.9, lb.amajor = 0.1, C.vec = NULL)
{
  # generate sink data by using mixed-proportion model, for simulation purpose (major sources have similar taxa distributions)
  # note, the unknown source is always not in major sources
  # Input: 
  #         N: scalar, number of taxa groups in source/sink
  #         K: scalar, number of KNOWN sources
  #         C: scalar, total taxa count in sink
  #         nmajor: scalar, number of major sources
  #         pmajor: scalar, sum of major proportions
  #         lb.amajor: scalar, lower bound of major proportions
  #         C.vec: vector, total taxa counts in K sources
  # Output: 
  #         a list with the following components:
  #           x.vec: vector of length N, sink data (check: sum up to C)
  #           beta.vec: vector of length N, abundance of sink (check: sum up to 1)
  #           gamma.mat: matrix of size (K+1) by N, abundance of KNOWN+UNKNOWN sources (note: sum up each row to 1)
  #                      note, the 1st row is reserved for the unobserved source
  #           alpha.vec: vector of length K+1, proportion of KNOWN+UNKNOWN sources (note: sum up to 1)
  #                      note, unobserved source is always put in the 1st element
  #           y.mat: matrix of size (K+1) by N, data of KNOWN+UNKNOWN sources
  #                      note, the 1st row is reserved for the unobserved source
  
  if (nmajor > K)
  {
    stop("Too many major sources! Need to increase K or decrease nmajor!")
  }
  if (nmajor * lb.amajor > pmajor)
  {
    stop("Bad configuration for major sources! Need to adjust nmajor, lb.amajor, and pmajor!")
  }
  gamma1.mat <- matrix(NA, nmajor, N)
  gamma1.vec <- rdirichlet(1, rep(1, N))
  gamma1.mat[1, ] <- gamma1.vec
  for (i in 2 : nmajor)
  {
    gamma1.mat[i, ] <- c(gamma1.vec[i], gamma1.vec[-i])
  }
  amajor.vec <- rdirichlet(1, rep(1, nmajor)) * pmajor
  while(min(amajor.vec) < lb.amajor)
  {
    amajor.vec <- rdirichlet(1, rep(1, nmajor)) * pmajor
  }
  nminor <- K + 1 - nmajor
  gamma2.mat <- rdirichlet(nminor, rep(1, N))
  pminor <- 1 - pmajor
  aminor.vec <- rdirichlet(1, rep(1, nminor)) * pminor
  alpha.vec <- c(amajor.vec, aminor.vec)
  alpha.vec <- c(alpha.vec[nmajor + 1], alpha.vec[-(nmajor + 1)]) # switch order so that the 1st element is for the unobserved minor source
  gamma.mat <- rbind(gamma1.mat, gamma2.mat)
  gamma.mat <- rbind(gamma.mat[nmajor + 1, ], gamma.mat[-(nmajor + 1), ]) # switch order so that the 1st element is for the unobserved minor source
  beta.vec <- alpha.vec %*% gamma.mat
  x.vec <- rmultinom(1, C, beta.vec)
  if (is.null(C.vec))
  {
    C.vec <- rep(C, K + 1)
  }
  y.mat <- round(gamma.mat * matrix(rep(C.vec, N), nrow = K + 1))
  return(list(x.vec = x.vec, beta.vec = beta.vec, gamma.mat = gamma.mat, alpha.vec = alpha.vec, y.mat = y.mat))
}

glsest = function(x.vec, yobs.mat)
{
  # estimate proportions of sources using GLS for multiple linear regression
  # Input: 
  #         x.vec: vector of length N, sink data
  #         yobs.mat: matrix of size K by N, data of KNOWN sources
  # Output: 
  #         a vector of estimated proportions for ALL K + 1 sources, note negative estimates are shrunk to 0
  #         note, the 1st element is for unknown source, it is the estimate of alpha0 * gamma0
  
  N <- length(x.vec)
  if (ncol(yobs.mat) != N)
  {
    stop("Dimension not match!")
  }
  K <- nrow(yobs.mat)
  C.vec <- rowSums(yobs.mat)
  C.mat <- matrix(rep(C.vec, N), nrow = K)
  gammao.mat <- yobs.mat / C.mat
  C <- sum(x.vec)

  gamma.mat <- rbind(1, gammao.mat) ### add intercept, this is for the UNKNOWN source
  beta.est.vec <- matrix(x.vec / C, ncol = 1)
  S.mat <- diag(as.vector(beta.est.vec)) - beta.est.vec %*% t(beta.est.vec)
  Sinv.mat <- psolve(S.mat)
  alpha.var.vec <- psolve(gamma.mat %*% Sinv.mat %*% t(gamma.mat))
  alpha.est.vec <- alpha.var.vec %*% (gamma.mat %*% Sinv.mat %*% x.vec)
  alpha.est.vec[alpha.est.vec < 0] <- 0
  alpha.est.vec <- as.vector(alpha.est.vec) / C
  alpha.se.vec <- sqrt(diag(alpha.var.vec)) / C
  return(list(alpha.est.vec, alpha.se.vec))
}

glsest1 = function(x.vec, yobs.mat)
{
  # estimate proportions of sources using GLS for multiple linear regression
  # Input: 
  #         x.vec: vector of length N, sink data
  #         yobs.mat: matrix of size K by N, data of KNOWN sources
  # Output: 
  #         a vector of estimated proportions for ALL K + 1 sources, note negative estimates are shrunk to 0
  #         note, the 1st element is for unknown source, it is the estimate of alpha0 * gamma0
  
  N <- length(x.vec)
  if (ncol(yobs.mat) != N)
  {
    stop("Dimension not match!")
  }
  K <- nrow(yobs.mat)
  C.vec <- rowSums(yobs.mat)
  C.mat <- matrix(rep(C.vec, N), nrow = K)
  gammao.mat <- yobs.mat / C.mat
  C = sum(x.vec)
  
  gamma.mat = rbind(1, gammao.mat) ### add intercept, this is for the UNKNOWN source
  beta.est.vec = matrix(x.vec / C, ncol = 1)
  S.mat = diag(as.vector(beta.est.vec)) - beta.est.vec %*% t(beta.est.vec)
  Sinv.mat = psolve(S.mat)
  mat = gamma.mat %*% Sinv.mat %*% t(gamma.mat)
  invmat = try(solve(mat), silent = T)
  
  alpha.est.vec = psolve(gamma.mat %*% Sinv.mat %*% t(gamma.mat)) %*% (gamma.mat %*% Sinv.mat %*% x.vec)
  alpha.est.vec[alpha.est.vec < 0] = 0
  alpha.est.vec = as.vector(alpha.est.vec) / C
  alpha.est.vec = alpha.est.vec / sum(alpha.est.vec)
  return(alpha.est.vec)
}

nll = function(alpha.vec, x.vec, gamma.mat)
{
  # calculate negative log-likelihood
  # Input: 
  #         alpha.vec: vector of length K+1, (estimated) proportion of ALL sources (no check: sum up to 1), note, the 1st element = alpha0 * gamma0 is for unknown source
  #         x.vec: vector of length N, sink data
  #         gamma.mat: matrix of size K by N, abundance of KNOWN sources (check: sum up each row to 1), note, abundance of unknown source should not be included
  # Output: 
  #         scalar, negative log-likelihood of observing sink data
  
  N = length(x.vec)
  K = nrow(gamma.mat)
  if (ncol(gamma.mat) != N || length(alpha.vec) != K + 1)
  {
    stop("Dimension not match!")
  }
  if (!all.equal(sum(gamma.mat), K))
  {
    stop("Check gamma.mat!")
  }
  gamma.mat = rbind(1, gamma.mat) ### add intercept, this is for the UNKNOWN source
  return(-(dmultinom(x.vec, sum(x.vec), as.vector(alpha.vec %*% gamma.mat), log = T)))
}

nll_s = function(alpha.vec, x.vec, gamma.mat)
{
  # calculate negative log-likelihood
  # Input: 
  #         alpha.vec: vector of length K+1, (estimated) proportion of ALL sources (no check: sum up to 1), note, the 1st element = alpha0 * gamma0 is for unknown source
  #         x.vec: vector of length N, sink data
  #         gamma.mat: matrix of size K+1 by N, abundance of ALL sources (check: sum up each row to 1), note, abundance of unknown source should not be included
  # Output: 
  #         scalar, negative log-likelihood of observing sink data
  
  if (length(x.vec) != ncol(gamma.mat) || length(alpha.vec) != nrow(gamma.mat))
  {
    stop("Dimension not match!")
  }
  if (!all.equal(sum(gamma.mat), length(alpha.vec)))
  {
    stop("Check gamma.mat!")
  }
  return(-(dmultinom(x.vec, sum(x.vec), as.vector(alpha.vec %*% gamma.mat), log = T)))
}

glsest_0 = function(x.vec, gammao.mat)
{
  # estimate proportions of sources using multiple linear regression
  # Input: 
  #         x.vec: vector of length N, sink data
  #         gammao.mat: matrix of size K by N, abundance of KNOWN sources (check: sum up each row to 1)
  # Output: 
  #         a vector of estimated proportions for ALL K + 1 sources, note negative estimates are shrunk to 0
  #         note, the 1st element is for unknown source, it is the estimate of alpha0 * gamma0
  
  N <- length(x.vec)
  if (ncol(gammao.mat) != N)
  {
    stop("Dimension not match!")
  }
  K <- nrow(gammao.mat)
  if (!all.equal(sum(gammao.mat), K))
  {
    stop("Check gamma.mat!")
  }
  C <- sum(x.vec)

  gamma.mat <- rbind(1, gammao.mat) ### add intercept, this is for the UNKNOWN source
  beta.est.vec <- matrix(x.vec / C, ncol = 1)
  S.mat <- diag(as.vector(beta.est.vec)) - beta.est.vec %*% t(beta.est.vec)
  Sinv.mat <- psolve(S.mat)
  alpha.var.vec <- psolve(gamma.mat %*% Sinv.mat %*% t(gamma.mat))
  alpha.est.vec <- alpha.var.vec %*% (gamma.mat %*% Sinv.mat %*% x.vec)
  alpha.est.vec[alpha.est.vec < 0] <- 0
  alpha.est.vec <- as.vector(alpha.est.vec) / C
  alpha.se.vec <- sqrt(diag(alpha.var.vec)) / C
  return(list(alpha.est.vec, alpha.se.vec))
}

DIRinfer = function(x.mat, gammaSS.mat)
{
  # infer direction from source-sink data
  # Input: 
  #         x.mat: matrix of size K+1 by N, source-sink data (note: the 1st row, ground truth for sink, sums up to C and each of the rest rows sums up to Ci)
  #         gammaSS.mat: matrix of size K+1 by N, abundance of KNOWN source and sink (check: sum up each row to 1)
  # Output: 
  #         scalar, index of sink from 1 : (K+1)
  
  K <- nrow(x.mat) - 1
  N <- ncol(x.mat)
  if (nrow(gammaSS.mat) != K + 1 || ncol(gammaSS.mat) != N)
  {
    stop("Dimension not match!")
  }
  if (!all.equal(sum(gammaSS.mat), K + 1))
  {
    stop("Check gammaSS.mat!")
  }
  njll.vec <- rep(NA, K + 1)
  for (i in 1 : (K + 1))
  {
    x.vec <- x.mat[i, ]
    gamma.mat <- gammaSS.mat[-i, ]
    estres <- glsest_0(x.vec, gamma.mat)
    alpha.est.vec <- estres[[1]]
    if (sum(is.na(alpha.est.vec)) == 0)
    {
      njll.vec[i] <- nll(alpha.est.vec, x.vec, gamma.mat)
      l <- 0
      for (j in setdiff(1 : (K + 1), i))
      {
        y.vec <- x.mat[j, ]
        l <- l + 1
        p.vec <- gamma.mat[l, ]
        njll.vec[i] <- njll.vec[i] - dmultinom(y.vec, sum(y.vec), p.vec, log = T)
      }
    }
  }
  ix <- which.min(njll.vec)
  return(ix)
}

error.bar <- function(x, y, upper, lower = upper, length = 0.1, ...)
{
  arrows(x, y + upper, x, y - lower, angle = 90, code = 3, length = length, ...)
}
