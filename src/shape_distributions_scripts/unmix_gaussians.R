#!/usr/bin/env Rscript
library(mixtools)

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
num_comp = as.integer(args[2])
outfile = args[3]

dist = read.csv(infile)
mixmdl = normalmixEM(dist$value, k=num_comp, maxrestarts=100)
save(mixmdl, file=outfile)
