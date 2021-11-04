args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])    

library(data.table)

gtf=args[1]
rd=args[2]

gtf<-fread(gtf)
colnames(gtf)[1]<-'V1'
gtf<-subset(gtf,V3=='exon')
gtf$V6<-as.numeric(gsub('.*exon_number "|"[;].*','',gtf$V9))
gtf$V7<-gsub('.*gene_name "|"[;].*','',gtf$V9)
gtf<-gtf[,c(1,4:7)]
colnames(gtf)<-c('chrom','start','end','exon_rank','gene_name')
save(gtf,file=rd)


