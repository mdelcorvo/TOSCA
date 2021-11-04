args <- (commandArgs(trailingOnly = TRUE))
isdata_table <- "data.table" %in% installed.packages()[, 1]
isRutils <- "R.utils" %in% installed.packages()[, 1]

if (!isdata_table) {install.packages("data.table", repos = "http://cran.rstudio.com/",lib = .libPaths()[1])}
if (!isRutils) {install.packages("R.utils", repos = "http://cran.rstudio.com/",lib = .libPaths()[1])}

library(data.table)
options(warn = -1)
vcf=args[1]
db=args[2]
outputfile=args[3]

cosmic<-fread(db,skip=grep('#CHROM',readLines(db))-1)
raw_vcf<-fread(vcf,skip=grep('#CHROM',readLines(vcf))-1)
colnames(cosmic)[1]<-'CHROM';colnames(raw_vcf)[1]<-'CHROM'
cosmic$id<-paste(cosmic$CHROM,cosmic$POS,sep='_');raw_vcf$id<-paste(raw_vcf$CHROM,raw_vcf$POS,sep='_')

cosmic<-cosmic[!duplicated(cosmic$id),]
raw_vcf$ID<- cosmic[match(raw_vcf$id,cosmic$id),]$ID

fwrite(raw_vcf,file=outputfile, row.names=F, col.names=T,quote=F, sep='\t')
