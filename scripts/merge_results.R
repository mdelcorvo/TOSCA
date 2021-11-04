args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])
if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])       

library(data.table)
options(warn = -1)

list1=args[1]
list1_ps=args[3]
outputfile=args[2]

if (!grepl('list',list1_ps)) {
	
	l1<-as.data.frame(fread(list1,header=F))
	res<-data.frame()
	for (i in l1$V1) {
		tmp<-fread(i)
		tmp$Sample<-gsub('.txt','',i)
		tmp<-tmp[,c(24,1:23)]
		res<-rbind(res,tmp)
	}
fwrite(res,outputfile, row.names=F, col.names=T,quote=F, sep='\t')

} else {
	
	l1<-as.data.frame(fread(list1,header=F))
	l2<-as.data.frame(fread(list1_ps,header=F))
	res<-data.frame()
	for (i in 1:length(l1$V1)) {
		tmp<-fread(l1$V1[i])
		tmp1<-fread(l2$V1[i])
		tmp$Sample<-gsub('.txt','',l1$V1[i])
		tmp$id<-paste(tmp$CHROM,tmp$POS,sep='_')
		tmp1$id<-paste(tmp1$chr,tmp1$start,sep='_')
		tmp2<-merge(tmp,tmp1[,c(26,27,53)],by.x='id',by.y='id',all.x=T)
		tmp3<-tmp2[,c(25,2:5,26,27,6:24)]
		colnames(tmp3)[6:7]<-c('PureCN_Prediction','PureCN_posterior probability')
		tmp3$PureCN_Prediction<-ifelse(is.na(tmp3$PureCN_Prediction),'germline',ifelse(tmp3$PureCN_Prediction==FALSE,'germline','somatic'))
		res<-rbind(res,tmp3)
	}
	fwrite(res,outputfile, row.names=F, col.names=T,quote=F, sep='\t')
	
}	
