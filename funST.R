library(gtools) # for function rdirichlet
library(npreg) # for function psolve

gendata = function(N, K, C, nmajor = 2, pmajor = 0.9, lb.amajor = 0.1)
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
  # Output: 
  #         a list with the following components:
  #           x.vec: vector of length N, sink data 
  #           beta.vec: vector of length N, abundance of sink 
  #           gamma.mat: matrix of size (K+1) by N, abundance of KNOWN+UNKNOWN sources (note: sum up each row to 1)
  #                      note, the 1st row is reserved for the unobserved source
  #           alpha.vec: vector of length K+1, proportion of KNOWN+UNKNOWN sources (note: sum up to 1)
  #                      note, unobserved source is always put in the 1st element
  
  if (nmajor > K)
  {
    stop("Too many major sources! Need to increase K or decrease nmajor!")
  }
  if (nmajor * lb.amajor > pmajor)
  {
    stop("Bad configuration for major sources! Need to adjust nmajor, lb.amajor, and pmajor!")
  }
  gamma.mat = rdirichlet(K + 1, rep(1, N))
  amajor.vec = rdirichlet(1, rep(1, nmajor)) * pmajor
  while(min(amajor.vec) < lb.amajor)
  {
    amajor.vec = rdirichlet(1, rep(1, nmajor)) * pmajor
  }
  nminor = K + 1 - nmajor
  pminor = 1 - pmajor
  aminor.vec = rdirichlet(1, rep(1, nminor)) * pminor
  alpha.vec = c(amajor.vec, aminor.vec)
  alpha.vec = c(alpha.vec[nmajor + 1], alpha.vec[-(nmajor + 1)]) # switch order so that the 1st element is for the unobserved minor source
  beta.vec = alpha.vec %*% gamma.mat
  x.vec = rmultinom(1, C, beta.vec)
  source.vec = rmultinom(1, C, alpha.vec)
  return(list(x.vec = x.vec, beta.vec = beta.vec, gamma.mat = gamma.mat, alpha.vec = alpha.vec))
}

glsest = function(x.vec, gamma.mat, useGLS = T)
{
  # estimate proportions of sources using multiple linear regression
  # Input: 
  #         x.vec: vector of length N, sink data
  #         gamma.mat: matrix of size K by N, abundance of KNOWN sources 
  #         useGLS: scalar, whether use GLS (TRUE) or OLS (FALSE)
  # Output: 
  #         a vector of estimated proportions for ALL K + 1 sources, note negative estimates are shrunk to 0
  #         note, the 1st element is for unknown source, it is the estimate of alpha0 * gamma0
  
  if (length(x.vec) != ncol(gamma.mat))
  {
    stop("Dimension not match!")
  }
  K = nrow(gamma.mat)
  if (!all.equal(sum(gamma.mat), K))
  {
    stop("Check gamma.mat!")
  }
  C = sum(x.vec)
  
  gamma.mat = rbind(1, gamma.mat) ### add intercept, this is for the UNKNOWN source
  if (useGLS)
  {
    beta.est.vec = matrix(x.vec / C, ncol = 1)
    S.mat = diag(as.vector(beta.est.vec)) - beta.est.vec %*% t(beta.est.vec)
    Sinv.mat = psolve(S.mat)
    alpha.est.vec = try(solve(gamma.mat %*% Sinv.mat %*% t(gamma.mat)) %*% (gamma.mat %*% Sinv.mat %*% x.vec), silent = T)
    if(inherits(alpha.est.vec, "try-error") | sum(alpha.est.vec < 0) >= K)
    {
      alpha.est.vec = rep(NA, K + 1)
    } else
    {
      alpha.est.vec = as.vector(alpha.est.vec) / C
      alpha.est.vec[alpha.est.vec < 0] = 0
    }
  } else
  {
    m = lm(x.vec ~ 0 + t(gamma.mat))
    alpha.est.vec = as.vector(coef(m)) / C
    alpha.est.vec[alpha.est.vec < 0] = 0
  }
  return(alpha.est.vec)
  
  # x.vec = res$x.vec
  # beta.est.vec = matrix(x.vec / C, ncol = 1)
  # S.mat = diag(as.vector(beta.est.vec)) - beta.est.vec %*% t(beta.est.vec)
  # tmp = chol(S.mat, pivot = T)
  # tmp1 = t(tmp)
  # S.res = tmp1 %*% tmp
  # S.res1 = tmp %*% tmp1
  # 
  # Sred.mat = S.mat[-1, -1]
  # system.time({Sredinv.mat = solve(Sred.mat)})
  # tmp = chol(S.mat, pivot = T)
  # 
  # data = data.frame(psink.est.vec, t(psource.mat))
  # mlr = lm(psink.est.vec ~ 0 + ., data = data)
  # # ah.raw.vec = summary(mlr)$coef[, 1]
  # ah.raw.vec = coef(mlr)
  # ah.raw.vec[is.na(ah.raw.vec)] = 0
  # ah.raw.vec = c(ah.raw.vec, sum(mlr$residuals)) ### newly added
  # ah.raw.vec[ah.raw.vec < 0] = 0
  # ah.vec = as.vector(ah.raw.vec / sum(ah.raw.vec))
  # return(ah.vec)
}

# N = 500; K = 100; C = 1e5
# nsimu = 1000
# mae.major.vec = rep(NA, nsimu)
# mae1.major.vec = rep(NA, nsimu)
# mae.all.vec = rep(NA, nsimu)
# mae1.all.vec = rep(NA, nsimu)
# set.seed(0)
# i = 0
# system.time({
# while(i < nsimu)
# {
#   res = gendata(N, K, C, nmajor = 2, pmajor = 0.9, lb.amajor = 0.1)
#   alpha.est.vec = glsest(res$x.vec, res$gamma.mat, useGLS = T)
#   if (sum(is.na(alpha.est.vec)) == 0)
#   {
#     i = i + 1
#     mae.major.vec[i] = mean(abs(alpha.est.vec - res$alpha.vec)[2 : 3])
#     mae.all.vec[i] = mean(abs(alpha.est.vec - res$alpha.vec))
#   }
# }
# })
# mean(mae.major.vec)
# mean(mae.all.vec)
# set.seed(0)
# system.time({
# for (i in 1 : nsimu)
# {
#   res = gendata(N, K, C, nmajor = 2, pmajor = 0.9, lb.amajor = 0.1)
#   alpha.est.vec = glsest(res$x.vec, res$gamma.mat, useGLS = F)
#   mae1.major.vec[i] = mean(abs(alpha.est.vec - res$alpha.vec)[2 : 3])
#   mae1.all.vec[i] = mean(abs(alpha.est.vec - res$alpha.vec))
# }
# })
# mean(mae1.major.vec)
# mean(mae1.all.vec)

# min(res$gamma.mat)
# sum(res$x.vec == 0)

nll = function(alpha.vec, x.vec, gamma.mat)
{
  # calculate negative log-likelihood
  # Input: 
  #         alpha.vec: vector of length K+1, (estimated) proportion of ALL sources, note, the 1st element = alpha0 * gamma0 is for unknown source
  #         x.vec: vector of length N, sink data
  #         gamma.mat: matrix of size K by N, abundance of KNOWN sources, note, abundance of unknown source should not be included
  # Output: 
  #         scalar, negative log-likelihood of observing sink data
  
  N = length(x.vec)
  K = nrow(gamma.mat)
  if (ncol(gamma.mat) != N || length(alpha.vec) != K + 1)
  {
    stop("Dimension not match!")
  }
  # if (sum(alpha.vec) != 1)
  # {
  #   stop("Check alpha.vec!")
  # }
  if (!all.equal(sum(gamma.mat), K))
  {
    stop("Check gamma.mat!")
  }
  gamma.mat = rbind(1, gamma.mat) ### add intercept, this is for the UNKNOWN source
  return(-(dmultinom(x.vec, sum(x.vec), as.vector(alpha.vec %*% gamma.mat), log = T)))
}

DIRinfer = function(x.mat, gammaSS.mat)
{
  # infer direction from source-sink data
  # Input: 
  #         x.mat: matrix of size K+1 by N, source-sink data (note: the 1st row, ground truth for sink, sums up to C and each of the rest rows sums up to Ci)
  #         gammaSS.mat: matrix of size K+1 by N, abundance of KNOWN source and sink 
  # Output: 
  #         scalar, index of sink from 1 : (K+1)
  
  K = nrow(x.mat) - 1
  N = ncol(x.mat)
  if (nrow(gammaSS.mat) != K + 1 || ncol(gammaSS.mat) != N)
  {
    stop("Dimension not match!")
  }
  if (!all.equal(sum(gammaSS.mat), K + 1))
  {
    stop("Check gammaSS.mat!")
  }
  njll.vec = rep(NA, K + 1)
  for (i in 1 : (K + 1))
  {
    x.vec = x.mat[i, ]
    gamma.mat = gammaSS.mat[-i, ]
    alpha.est.vec = glsest(x.vec, gamma.mat, useGLS = T)
    if (sum(is.na(alpha.est.vec)) == 0)
    {
      njll.vec[i] = nll(alpha.est.vec, x.vec, gamma.mat)
      l = 0
      for (j in setdiff(1 : (K + 1), i))
      {
        y.vec = x.mat[j, ]
        l = l + 1
        p.vec = gamma.mat[l, ]
        njll.vec[i] = njll.vec[i] - dmultinom(y.vec, sum(y.vec), p.vec, log = T)
      }
    }
  }
  ix = which.min(njll.vec)
  # print(njll.vec)
  return(ix)
}
