#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(gtools)))
suppressMessages(suppressWarnings(library(npreg)))

setwd("./") # make sure the FastST folder is working directory
source("./Code/funST.R")

res_dirname <- "./Result"
if (!dir.exists(res_dirname)){
 dir.create(res_dirname)
}

# Set parameters
N <- 500; K <- 10; C <- 1e5; nmajor <- 5; pmajor <- 0.9; lb.amajor <- 0.1

# Define command-line options
option_list <- list(
  make_option(c("-O", "--filename"), type = "character", default = "demo",
              help = "File name of simulated data and result [default = %default]")
)

# Create parser object and parse arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
datafile <- paste0(args$filename, "_data.txt")
resfile <- paste0(args$filename, "_result.txt")
figfile <- paste0(args$filename, ".jpg")

cat("Simulating sink data using the following settings:\n")
cat("-- Number of taxa:", N, "\n")
cat("-- Number of sources:", K, "\n")
cat("-- Total taxa counts:", C, "\n")
cat(paste0("-- Number of major sources: ", nmajor, "; From #", 2, " to #", nmajor + 1, "\n"))
cat("-- Sum of major proportions:", pmajor, "\n")
cat("-- Lower bound of major proportions:", lb.amajor, "\n\n")

# Generate sink data
set.seed(123)
res <- gendata(N, K, C, nmajor = nmajor, pmajor = pmajor, lb.amajor = lb.amajor)
x.vec <- res$x.vec
beta.vec <- res$beta.vec
gamma.mat <- res$gamma.mat
alpha.vec <- res$alpha.vec
y.mat <- res$y.mat
yobs.mat <- y.mat[-1, ]

data.mat <- cbind(x.vec, t(yobs.mat)) # 1st: sink, rest: observed sources
write.table(data.mat, file = paste0("./Result/", datafile), sep = "\t", row.names = F, col.names = F)
cat(paste0("Data simulated. File name: ./Result/", datafile), "\n\n")

cat("Estimating source proportions...\n")

runt <- system.time(estres <- glsest(x.vec, yobs.mat))
alpha.est.vec <- estres[[1]]
alpha.se.vec <- estres[[2]]
prop.mat <- cbind(alpha.vec, alpha.est.vec, alpha.se.vec)
cat("Source proportion estimation result:\n")
print(prop.mat)
write.table(prop.mat, file = paste0("./Result/", resfile), sep = "\t", row.names = F, col.names = F)
cat(paste0("\nResult saved. File name: ./Result/", resfile), "\n")
cat("Elapsed time:", runt[3], "seconds \n\n")

jpeg(filename = paste0("./Result/", figfile))
alpha.lst <- data.frame(name = 1 : K, alpha = alpha.vec[-1], alpha.est = alpha.est.vec[-1])
ub <- 1.25 * max(alpha.lst$alpha)
bardata <- t(alpha.lst[, 2 : 3])
barplot(bardata, axes = F, beside = T, legend.text = c("True proportions", "Estimated proportions"), col = c("blue", "skyblue"), ylim = c(-0.02, ub), ylab = "Proportion of source contributions")
axis(1, at = seq(2, 3 * K, by = 3), labels = 1 : K, line = 1)
axis(2)
lines(c(1, 3 * nmajor), c(-0.01, -0.01), col = "black")
mtext("major sources", side = 1, line = -0.5, at = 3 * nmajor / 2)
suppressWarnings({error.bar(3 * (1 : K) - 0.5, alpha.est.vec[-1], alpha.se.vec[-1])})
trash <- dev.off()
cat(paste0("Figure generated. File name: ./Result/", figfile), "\n")
