args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])
if (!requireNamespace("PureCN", quietly = TRUE))
    BiocManager::install("PureCN")    
if (!requireNamespace("data.table", quietly = TRUE))
    BiocManager::install("data.table")  
    
library(PureCN)
library(data.table)

baits=args[1]
reference=args[2]
mappability=args[3]
gtf= args[4]
min_target=as.numeric(args[5])
min_off_target=as.numeric(args[6])
intervals=args[7]

baits <- import(baits)
mappability <- import(mappability)

interval_tmp <- as.data.frame(preprocessIntervals(baits, reference, mappability = mappability, off.target=T, min.target.width = min_target, min.off.target.width= min_off_target))

#annotate Targets with GTF file from Ensembl...
gtf<-fread(gtf)
colnames(gtf)[1]<-'V1'
gene<-subset(gtf,V3=='gene')
gene <- gene[gene$V1 %in% c(1:22,'X','Y'),]
gene$V9<-gsub('gene_id ["]','',gsub('.*gene_name "|"[;].*','',gene$V9))

gr_interval<-GRanges(interval_tmp$seqnames, ranges =IRanges(interval_tmp$start,interval_tmp$end))
gr_gene <- GRanges(gene$V1, ranges =IRanges(gene$V4,gene$V5))
overlap_gr = findOverlaps(query = gr_gene, subject = gr_interval)

interval_tmp<-data.frame(interval_tmp[subjectHits(overlap_gr),], gene[queryHits(overlap_gr),])
interval_tmp$Target<-paste(paste(interval_tmp$seqnames,interval_tmp$start,sep=':'),interval_tmp$end,sep='-')
interval_tmp<-interval_tmp[!duplicated(interval_tmp$Target),]
interval_tmp$seqnames<- factor(interval_tmp$seqnames,levels =c(1:22,'X','Y'))
interval_tmp <- interval_tmp[with(interval_tmp, order(interval_tmp$seqnames,interval_tmp$start)), ]
interval_tmp$Gene <- interval_tmp$V9

res <- interval_tmp[,c(20,9,7:8,10,6)]
colnames(res)[6]<- 'on_target'
fwrite(res,file=intervals, row.names=F, col.names=T,quote=F, sep='\t')
