workpath = "./"
source(paste(workpath, "funST.R", sep = ""))

N = 500; C = 1e5
K.vec = c(2, 5, 10, 50, 100)
nk = length(K.vec)
nmajor.vec = c(2, 5)
nn = length(nmajor.vec)
nsimu = 1000
accr.mat = matrix(NA, nk, nn)
runt.mat = matrix(NA, nk, nn)

set.seed(0)
for (i in 1 : nk)
{
  K = K.vec[i]
  for (j in 1 : nn)
  {
    nmajor = nmajor.vec[j]
    sinkid.est.vec = rep(NA, nsimu)
    if (nmajor <= K)
    {
      kk = 0
      runt = system.time({
        while(kk < nsimu)
        {
          res = gendata(N, K, C, nmajor = nmajor, pmajor = 0.9, lb.amajor = 0.1)
          x.mat = as.vector(res$x.vec)
          gammaSS.mat = rbind(res$beta.vec, res$gamma.mat[-1, ])
          for (l in 2 : (K + 1))
          {
            x.mat = rbind(x.mat, as.vector(rmultinom(1, C, res$gamma.mat[l, ])))
          }
          kk = kk + 1
          sinkid.est.vec[kk] = DIRinfer(x.mat, gammaSS.mat)
        }
      })
      accr.mat[i, j] = mean(sinkid.est.vec == 1)
      runt.mat[i, j] = runt[3]
      print(paste("K = ", K, ", nmajor = ", nmajor, ", time = ", runt[3], sep = ""))
    }
  }
}
