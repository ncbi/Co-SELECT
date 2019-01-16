#!/usr/bin/env Rscript
library(vcdExtra)

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outfile = args[2]

df = read.table(infile)
names(df) = c('value', 'count')
df$count <- round(df$count*1.0/min(df$count),0)
exdf = expand.dft(df, freq="count")
write.csv(exdf, file=outfile, row.names=F)
