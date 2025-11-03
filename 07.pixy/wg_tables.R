#! /sw/apps/R/4.1.1/rackham/bin/Rscript
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

prefix <- args[1]
outname <- args[2]

files<-list.files(pattern=prefix)
tables<-lapply(files,read.table, header=T)
dxy<-rbindlist(tables,use.names=TRUE)
write.table(dxy,file=outname,col.names = T,row.names=F,quote=F,sep="\t")