args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])
if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])       
if (!requireNamespace("GenomicRanges", quietly = TRUE))
    BiocManager::install("GenomicRanges")
    
library(data.table)
library(GenomicRanges)
options(warn = -1)
#Function to assemble SnpEff output
snpeff_f <- function(vcf) {
tmp1<-strsplit(vcf$INFO,'|',fixed = TRUE)
tmp2<- tmp1[grepl('HIGH|MODERATE|LOW|MODIFIER',tmp1)]
list1<- lapply(tmp2,function(x) c(length(grep('HIGH',  x)),length(grep('MODERATE',  x)),length(grep('LOW',  x)),length(grep('MODIFIER',  x))))
data1<- do.call(rbind.data.frame, list1)
colnames(data1)<-c('HIGH','MODERATE','LOW','MODIFIER')
data1<-cbind(vcf[grepl('HIGH|MODERATE|LOW|MODIFIER',tmp1),1:5],data1)
data1$Tot_transcript <- rowSums(data1[,6:9])
data1$Main_impact<-ifelse(data1$HIGH>0,'HIGH',ifelse(data1$MODERATE>0,'MODERATE',ifelse(data1$LOW>0,'LOW','MODIFIER')))
data1$Effect<-unlist(lapply(tmp2,function(x) {tmp<-c(x[grep('HIGH',  x)-1][1],x[grep('MODERATE',  x)-1][1],x[grep('LOW',  x)-1][1],x[grep('MODIFIER',  x)-1][1]);tmp<-tmp[!is.na(tmp)][1];return(tmp)}))
data1$Gene<-unlist(lapply(tmp2,function(x) {tmp<-c(x[grep('HIGH',  x)+1][1],x[grep('MODERATE',  x)+1][1],x[grep('LOW',  x)+1][1],x[grep('MODIFIER',  x)+1][1]);tmp<-tmp[!is.na(tmp)][1];return(tmp)}))
data1$Ensembl_id<-unlist(lapply(tmp2,function(x) {tmp<-c(x[grep('HIGH',  x)+2][1],x[grep('MODERATE',  x)+2][1],x[grep('LOW',  x)+2][1],x[grep('MODIFIER',  x)+2][1]);tmp<-tmp[!is.na(tmp)][1];return(tmp)}))
data1$DNA_type<-unlist(lapply(tmp2,function(x) {tmp<-c(grep('HIGH',  x)[1],grep('MODERATE',  x)[1],grep('LOW',  x)[1],grep('MODIFIER',  x)[1]);tmp<-tmp[!is.na(tmp)][1];tmp<-grep(paste(c('A>G','G>A','A>T','T>A','A>C','C>A','G>T','T>G','G>C','C>G','C>T','T>C','c\\.'),collapse='|'),x[tmp:(tmp+15)],value=T)[1];return(tmp)}))
data1$Protein_type<-unlist(lapply(tmp2,function(x) {tmp<-c(grep('HIGH',  x)[1],grep('MODERATE',  x)[1],grep('LOW',  x)[1],grep('MODIFIER',  x)[1]);tmp<-tmp[!is.na(tmp)][1];tmp<-grep(paste(c('Gly','Ala','Val','Leu','Ile','Met','Cys','Pro','Phe','Trp','Tyr','Thr','Ser','Asn','Gln','Asp','Glu','His','Lys','Arg','p\\.'),collapse='|'),x[tmp:(tmp+15)],value=T)[1];return(tmp)}))
colnames(data1)[1]<-'CHROM'
data1$ID<-NULL
return(data1)
}
#Function to check the perfect match between raw and annotation files
fix_vcf <- function(vcf,tmp) {
if (nrow(vcf)!=nrow(tmp)) {
	colnames(vcf)[1]<-'chrom'
	vcf$id<-paste(vcf$chrom,vcf$POS,sep='_')
	colnames(tmp)[1]<-'chrom'
	tmp$id<-paste(tmp$chrom,tmp$POS,sep='_')
	c1<-as.numeric(row.names(vcf[is.na(match(vcf$id,tmp$id)),]))
	 for (i in c1) {
		 if (i == nrow(vcf)) {
			tmp<- rbind(tmp[1:(i-1),],vcf[i,])	 
		 } else {
			tmp<- rbind(tmp[1:(i-1),],vcf[i,],tmp[i:nrow(tmp),])	
	     } 
	 }
return(tmp)	 
} else {
return(tmp)	
}		 
}	

vcf=args[1]
snpeff=args[2]
gen1k=args[3]
esp=args[4]
exac=args[5]
dbsnp=args[6]
cosmic=args[7]
clinvar=args[8]
cov=args[9]
af=args[10]
gtf=args[11]
ref=args[12]
depth=as.numeric(args[13])
vaf=as.numeric(args[14])
outputfile=args[15]

raw_res<- as.data.frame(fread(vcf,skip=grep('#CHROM',readLines(vcf))-1))
snpeff<-fix_vcf(raw_res,fread(snpeff,skip=grep('#CHROM',readLines(snpeff))-1))
gen1k<-fix_vcf(raw_res,fread(gen1k,skip=grep('#CHROM',readLines(gen1k))-1))
esp<-fix_vcf(raw_res,fread(esp,skip=grep('#CHROM',readLines(esp))-1))
exac<-fix_vcf(raw_res,fread(exac,skip=grep('#CHROM',readLines(exac))-1))
dbsnp<-fix_vcf(raw_res,fread(dbsnp,skip=grep('#CHROM',readLines(dbsnp))-1))
cosmic<-fix_vcf(raw_res,fread(cosmic,skip=grep('CHROM',readLines(cosmic))-1))
clinvar<-fix_vcf(raw_res,fread(clinvar,skip=grep('#CHROM',readLines(clinvar))-1))
cov<-fread(cov)
af<-fread(af)
load(gtf)

####################################################################
#The tumor-only filtration algorithm, inspired by Sukhai et al. 2019
#https://www.sciencedirect.com/science/article/pii/S1525157817305986
####################################################################

##############################################################################################################################################
#Step 1 
# 1.1 Retained by criteria of quality pass, read depth > cutoff (250x suggested) and VAF > cutoff (5% suggested)
# 1.2 Retained by criteria of variant type (non-synonymous) and in regions of interest (exons or first 2 base pairs of introns flanking exons)
##############################################################################################################################################
#Data recovery for step 1
res<-snpeff_f(snpeff)
res$Depth<-cov$avgcov_0
res$Type<-af[,3]
res$AF<-af[,4]
if (sum(grepl(',',res$AF))>0) {
	res$AF<-unlist(lapply(strsplit(res$AF,',',fixed = TRUE),function(x) round(mean(as.numeric(x)),digits=3)))
}

###################################################################################################################################
#Step 2 
# - Retained by criteria of not present in germline population variant databases (PVDs) with MAF >1%, or present in COSMIC database
###################################################################################################################################
#Data recovery for step 2
res$dbSNP<-ifelse(dbsnp$ID!='.',dbsnp$ID,NA)
res$ESP<-ifelse(esp$ID!='.',round(as.numeric(esp$ID),digits=4),NA)
res$ExAC<-ifelse(exac$ID!='.',round(as.numeric(exac$ID),digits=4),NA)
res$Genome1k<-ifelse(gen1k$ID!='.',round(as.numeric(gsub('.*MAF=|[;].*','',gen1k$INFO)),digits=4),NA)
res$COSMIC<-cosmic$ID

####################################################################################################################
#Step 3 
# - Retained by criteria of not present in either or ClinVar databases with a benign or likely benign classification
####################################################################################################################
#Data recovery for step 3
res$ClinVarDisease_name<-ifelse(clinvar$ID!='.',gsub('.*CLNDN=|;.*','',clinvar$INFO),NA)
res$ClinVarClinical_significance<-ifelse(clinvar$ID!='.',gsub('.*CLNSIG=|;.*','',clinvar$INFO),NA)

#######################
#Application of 3 steps
#######################

res$Prediction<- ifelse(res$AF < vaf | res$Depth < depth, 'germline', #Step 1.1
					ifelse(grepl(paste(c("synonymous_variant","start_retained","stop_retained_variant"),collapse="|"),res$Effect) | grepl('MODIFIER',res$Main_impact),'germline', #Step 1.2
						ifelse((((!is.na(res$ESP) & res$ESP>0.01) | (!is.na(res$ExAC) & res$ExAC>0.01) | (!is.na(res$Genome1k) & res$Genome1k > 0.01)) & res$COSMIC=='') |  #Step 2.1
							((!is.na(res$Genome1k) & res$Genome1k > 0.0001) & (!is.na(res$ESP) & res$ESP>0.0001) & (!is.na(res$ExAC) & res$ExAC>0.0001)) |  #Step 2.2
								((!is.na(res$dbSNP)) & (!is.na(res$ESP) & res$ESP>0.0001) & (!is.na(res$ExAC) & res$ExAC>0.0001))  |  #Step 2.3
									((!is.na(res$dbSNP)) & (!is.na(res$Genome1k) & res$Genome1k>0.0001) & (!is.na(res$ExAC) & res$ExAC>0.0001)) |  #Step 2.4
										((!is.na(res$dbSNP)) & (!is.na(res$Genome1k) & res$Genome1k>0.0001) & (!is.na(res$ESP) & res$ESP>0.0001)),'germline', #Step 2.5
							ifelse(grepl('Benign|Likely_benign',res$ClinVarClinical_significance),'germline','somatic')))) #Step 3

#######################
res$Contamination<-raw_res$FILTER
res$refGenome<-ref
res <- res[with(res, order(res$CHROM, res$POS)), ]


#############
#exon ranking
#############
res<-as.data.frame(res)
res$exon_rank<-NA
exon <- GRanges(gtf$chrom, ranges =IRanges(gtf$start,gtf$end))
tmp<-GRanges(res$CHROM, ranges =IRanges(res$POS,res$POS))
exon_match = findOverlaps(query = exon, subject = tmp)
res[subjectHits(exon_match),]$exon_rank<- gtf[queryHits(exon_match),]$exon_rank
res<-res[,c('Gene','Ensembl_id','Type','CHROM','POS','REF','ALT','refGenome','Prediction','Contamination','exon_rank','AF','Depth','Effect','dbSNP','ESP','ExAC','Genome1k','COSMIC','ClinVarDisease_name','ClinVarClinical_significance','Main_impact','DNA_type','Protein_type')]
colnames(res)[18]<-'1000Genome'

res<-res[order(res$Prediction,decreasing = T),]
fwrite(res,outputfile, row.names=F, col.names=T,quote=F, sep='\t')
