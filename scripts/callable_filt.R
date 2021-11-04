args <- (commandArgs(trailingOnly = TRUE))

if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])

library(data.table)
options(warn = -1)

bed=args[1]
outputfile=args[2]

bed <- fread(bed)
bed <- bed[bed$V1 %in% c(1:22,'X','Y','MT'),]

fwrite(bed,file=outputfile, row.names=F, col.names=F,quote=F, sep='\t')
