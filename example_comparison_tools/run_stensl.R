library(FEAST)

K.vec = c(10)
nk = length(K.vec)

nsimu = 1000
MAX_ITER = 1000
nmajor.vec = c(2, 5)
nn = length(nmajor.vec)

work_dir <- "/home/joungmin/mst/"
for (i in 1 : nk){
  setwd(work_dir)
  K = K.vec[i]
  meta_file <- paste0("/home/joungmin/mst/metadata_feast/metadata_K_", K, ".txt")
  metadata <- read.table(meta_file, sep='\t', header=T)
  rownames(metadata) = metadata$SampleID

  for (j in 1 : nn)
  {
    setwd(work_dir)
    nmajor = nmajor.vec[j]
    res_dirname <- paste0("/home/joungmin/mst/STENSL/results_STENSL/K_", K, "_nmajor_", nmajor, "/")
    if (!dir.exists(res_dirname)){
      dir.create(res_dirname)
    }
    if (nmajor <= K){
      kk = 0
      runt = system.time({
      while(kk < nsimu){
        setwd(work_dir)
        otu_file <- paste0("/home/joungmin/mst/simul_data/K_", K, "_nmajor_", nmajor, "/simul_", kk+1, ".txt")
        print(otu_file)
        otus = read.table(otu_file, sep='\t', header=T)
        sink_sample <- "Sink"
        source_samples <- setdiff(colnames(otus), sink_sample)

        
        samples_to_check <- c(sink_sample, source_samples)
        col_depths <- colSums(otus[, samples_to_check])
        COVERAGE_DEPTH <- min(col_depths)
        otus <- t(otus)

        res_file <- paste0(res_dirname, 'result_simul_', kk+1, ".txt")
        result <- STENSL(C=as.matrix(otus), metadata=metadata, EM_iterations=MAX_ITER, COVERAGE=COVERAGE_DEPTH)
        FEAST_output <- result$proportions_mat
        write.table(FEAST_output, res_file, row.names = T, quote = F, sep = '\t')
        kk = kk + 1
      }
      })
      runt.df <- data.frame(v1 = runt[3])
      setwd(work_dir)
      runt_file <- paste0('/home/joungmin/mst/STENSL/results_STENSL/runtime_feast_K_',K,'_nmajor_', nmajor, '.csv')
      write.csv(runt.df, runt_file, row.names = F, quote = F)
    }
  }
}


  
