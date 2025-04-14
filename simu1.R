workpath = "./"
setwd(workpath)
source(paste(workpath, "funST.R", sep = ""))


N = 500; C = 1e5
K.vec = c(2, 5, 10, 50, 100)

nk = length(K.vec)
nmajor.vec = c(2, 5)
nn = length(nmajor.vec)
nsimu = 1000
amae1.observed.mat = matrix(NA, nk, nn) # GLS
amae1.unobserved.mat = matrix(NA, nk, nn)
amae1.major.mat = matrix(NA, nk, nn)
runt1.mat = matrix(NA, nk, nn)
set.seed(0)

for (i in 1 : nk){
  K = K.vec[i]
  for (j in 1 : nn)
  {
    nmajor = nmajor.vec[j]
    mae.observed.vec = rep(NA, nsimu)
    mae.unobserved.vec = rep(NA, nsimu)
    mae.major.vec = rep(NA, nsimu)
    if (nmajor <= K)
    {
      kk = 0
      runt = system.time({
        while(kk < nsimu)
        {
          res = gendata(N, K, C, nmajor = nmajor, pmajor = 0.9, lb.amajor = 0.1)        
          alpha.est.vec = glsest(res$x.vec, res$gamma.mat[-1, ], useGLS = T)
          if (sum(is.na(alpha.est.vec)) == 0)
          {
            mae.observed.vec[kk] = mean(abs(alpha.est.vec[-1] - res$alpha.vec[-1]))
            mae.unobserved.vec[kk] = abs(alpha.est.vec[1] - res$alpha.vec[1] * mean(res$gamma.mat[1, ]))
            mae.major.vec[kk] = mean(abs(alpha.est.vec[2 : 3] - res$alpha.vec[2 : 3]))
          }
        }
      })
      amae1.observed.mat[i, j] = mean(mae.observed.vec)
      amae1.unobserved.mat[i, j] = mean(mae.unobserved.vec)
      amae1.major.mat[i, j] = mean(mae.major.vec)
      runt1.mat[i, j] = runt[3]
    }
  }
}

