args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])
if (!requireNamespace("PureCN", quietly = TRUE))
    BiocManager::install("PureCN")    

library(PureCN)

vcf= args[1]
list= args[2]
rds= args[3]
rds_bias= args[4]
ref= args[5]

l1<-read.table(list)
normalDB <- createNormalDatabase(l1$V1)
saveRDS(normalDB, file = rds)
bias <- calculateMappingBiasVcf(vcf, genome = ref)
saveRDS(bias, rds_bias)


