args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])       

library(data.table)

inputfile=args[1]
outputfile=args[2]

map<-fread(inputfile,header=F)
map$V1<-gsub('chr','',map$V1)
map <- map[map$V1 %in% c(1:22,'X','Y'),]
fwrite(map,outputfile, row.names=F, col.names=F,quote=F, sep='\t')
