library(FEAST)

K.vec = c(10)
nk = length(K.vec)

nsimu = 1000
nmajor.vec = c(2, 5)
nn = length(nmajor.vec)

work_dir <- "/data2/project/joungmin/MST/exp_v7_020325/dataset"
for (i in 1 : nk){
  K = K.vec[i]
  meta_file <- paste0("/data2/project/joungmin/MST/exp_v7_020325/dataset/metadata_feast/metadata_K_", K, ".txt")
  metadata <- Load_metadata(metadata_path = meta_file)
  for (j in 1 : nn)
  {
    nmajor = nmajor.vec[j]
    res_dirname <- paste0("/data2/project/joungmin/MST/exp_v7_020325/dataset/results_FEAST/K_", K, "_nmajor_", nmajor)
    if (!dir.exists(res_dirname)){
      dir.create(res_dirname)
    }
    if (nmajor <= K){
      kk = 0
      runt = system.time({
      while(kk < nsimu){
        otu_file <- paste0("/data2/project/joungmin/MST/exp_v7_020325/dataset/simul_data/K_", K, "_nmajor_", nmajor, "/simul_", kk+1, ".txt")
        print(otu_file)
        otus <- Load_CountMatrix(otu_file)
        res_file <- paste0('result_simul_', kk+1)
        FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, 
                                                  dir_path = paste0(res_dirname,"/"), outfile=res_file)
        kk = kk + 1
      }
      })
      runt.df <- data.frame(v1 = runt[3])
      runt_file <- paste0('/data2/project/joungmin/MST/exp_v7_020325/dataset/runtime_feast_K_',K,'_nmajor_', nmajor, '.csv')
      write.csv(runt.df, runt_file, row.names = F, quote = F)
    }
  }
}


  
