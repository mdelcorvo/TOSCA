args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])
if (!requireNamespace("PureCN", quietly = TRUE))
    BiocManager::install("PureCN")    
if (!requireNamespace("data.table", quietly = TRUE))
    BiocManager::install("data.table")
    
library(PureCN)
library(data.table)

#input file
rds= args[1]
cov_tum= args[2]
vcf= args[3]
interval= args[4]
rds_bias= args[5]
#params
ref= args[6]
sample= args[7]
#output file
rds_output= args[8]
pS_output= args[9]
plot_output= args[10]

normalDB <- readRDS(rds)
pool <- calculateTangentNormal(cov_tum, normalDB)

ret <-try(runAbsoluteCN(normal.coverage.file = pool,tumor.coverage.file = cov_tum, 
 vcf.file = vcf,genome = ref, sampleid = sample,interval.file = interval,
  args.setMappingBiasVcf=list(mapping.bias.file = rds_bias),
   normalDB = normalDB, post.optimize = TRUE, plot.cnv = F, verbose = TRUE,model.homozygous=TRUE))

if (class(ret) == "try-error"){
normalDB<-NULL
ret <-runAbsoluteCN(normal.coverage.file = pool,tumor.coverage.file = cov_tum, 
 vcf.file = vcf,genome = ref, sampleid = sample,interval.file = interval,
  args.setMappingBiasVcf=list(mapping.bias.file = rds_bias),
   normalDB = normalDB, post.optimize = TRUE, plot.cnv = F, verbose = TRUE,model.homozygous=TRUE)
} 
    
seqlevelsStyle(ret$input$centromeres)<-"NCBI"   
pS= predictSomatic(ret); pS$Sample<-sample
pdf(plot_output, width = 16, height = 11)
plotAbs(ret, type = "overview")
plotAbs(ret, 1, type = "hist")
plotAbs(ret, 1, type = "AF")
plotAbs(ret, 1, type = "BAF")
dev.off()

fwrite(pS,file=pS_output, row.names=F, col.names=T,quote=F, sep='\t')
saveRDS(ret, rds_output)
