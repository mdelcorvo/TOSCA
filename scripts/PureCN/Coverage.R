args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])
if (!requireNamespace("PureCN", quietly = TRUE))
    BiocManager::install("PureCN")    

library(PureCN)

interval.file<-args[1]
bam.file<-args[2]
coverage.file<-args[3]
coverage.loess.file<- args[4]

calculateBamCoverageByInterval(bam.file = bam.file,interval.file = interval.file, output.file = coverage.file)
correctCoverageBias(coverage.file, interval.file,output.file = coverage.loess.file, plot.bias = F)

