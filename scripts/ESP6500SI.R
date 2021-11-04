args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])

library(data.table)
options(warn = -1)

build=args[1]

setwd(paste('resources/database/',build,sep=''))
fileNames <- Sys.glob("*.vcf")
for (i in 1:length(fileNames)) {	
nex <- fread(fileNames[i],skip=as.numeric(grep("CHROM", readLines(fileNames[i]))),sep='\t')
nex$V3<- as.numeric(gsub('.*[,]','',gsub('.*MAF=|[;].*','',nex$V8)))/100
	if (build=='GRCh38') {
	nex$V1<-gsub('.*GRCh38_POSITION=|[:].*','',nex$V8)
	nex$V2<-gsub('.*GRCh38_POSITION=|.*[:]','',nex$V8)
	}
colnames(nex)<-c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
nex$QUAL<-'.';nex$FILTER<-'.';nex$INFO<-'.';
fwrite(nex,file=fileNames[i], row.names=F, col.names=T,quote=F, sep='\t')
}

