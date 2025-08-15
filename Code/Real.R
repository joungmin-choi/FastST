# Real data analysis using the Knights et al. dataset
suppressMessages(suppressWarnings(library(npreg)))
suppressMessages(suppressWarnings(library(ggplot2)))
                 
setwd("./") # make sure the FastST folder is working directory
source("./Code/funST.R")

res_dirname <- "./Result"
if (!dir.exists(res_dirname)){
 dir.create(res_dirname)
}

# Read data
sinkfile <- "./Real/otus_sink.csv"
sourcefile <- "./Real/otus_source.csv"
metafile <- "./Real/metadata-subset.txt"
rawsink <- read.table(sinkfile, header = F, sep = ",")
sink <- rawsink[, -1]
nsink <- ncol(sink)
sink.vec <- c("Lab1", "Lab2", "NICU", "Office")

rawsource <- read.table(sourcefile, header = F, sep = ",")
source <- rawsource[, -1]

rawmeta <- read.table(metafile, header = F, sep = "\t", skip = 5)
class.vec <- unlist(unique(rawmeta[3]))
nclass <- length(class.vec)
ix.lst <- vector("list", nclass)
for (i in 1 : nclass)
{
  ix.lst[[i]] <- which(rawmeta[3] == class.vec[i])
}

# Estimate proportions
est.mat <- data.frame(matrix(NA, nsink * nclass, 3))
colnames(est.mat) <- c("Proportion", "Env", "Sample")
est.mat$Env <- rep(class.vec, nsink)
est.mat$Sample <- rep(sink.vec, each = nclass)
yobs.mat <- t(source)
runt.vec <- rep(NA, nsink)
for (i in 1 : nsink)
{
  x.vec <- sink[, i]
  runt <- system.time({
    suppressWarnings({estres <- glsest(x.vec, yobs.mat)})
  })
  alpha.est.vec <- estres[[1]]
  for (j in 1 : nclass)
  {
    temp <- sum(alpha.est.vec[ix.lst[[j]]])
    est.mat$Proportion[(i - 1) * nclass + j] <- temp / sum(alpha.est.vec)
  }
  runt.vec[i] <- runt[3]
  # print(paste0("Sink ", i, " done! Time = ", runt[3]))
}

# Save results
resfile <- "./Result/real_result.txt"
cat("Source proportion estimation result:\n")
print(est.mat)
write.table(est.mat, file = resfile, sep = "\t", row.names = F, col.names = T)
cat(paste0("\nResult saved. File name: ", resfile), "\n")
cat("Elapsed time:", sum(runt.vec), "seconds \n\n")

# Save plot
figfile <- "./Result/real.jpg"
gp <- ggplot(data = est.mat, aes(x = Sample, y = Proportion, fill = Env)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  labs(title = "FastST", x = "Sink") +
  theme(axis.text.x = element_text(), plot.title = element_text(hjust = 0.5)) #, hjust = 1
#angle = 45, hjust = 1
# x11()
# postscript(file = figfile)
jpeg(filename = figfile)
plot(gp)
trash <- dev.off()
cat(paste0("Figure generated. File name: ", figfile), "\n")
