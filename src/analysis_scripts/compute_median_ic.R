#!/usr/bin/env Rscript

infile = '../data/homeodomain/NKX2-3/TGGAAT20NGA/4.non.CACT.d2.bg.MGW.S4.000M4.600H5.600X.mer.pwm'
outfile = 'haha.scv'

args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  infile = args[1]
  outfile = args[2]
}

df = read.csv(infile)

require(plyr)
require(dplyr)
require(tidyr)

all_cols = names(df)
mat_cols = all_cols[grepl("X\\d+", all_cols)]
group_cols = setdiff(all_cols, c('base', mat_cols))

getPWM <- function(df) {
    mat = as.matrix(df[, mat_cols])
    rownames(mat) = df$base
    return(mat)
}

computeBits <- function(pfm, N=4, Nseqs=NULL){
  pwm = apply(pfm, 2, function(x) x / sum(x, na.rm=T))
  Nseqs = attr(pwm, 'nongapped') # TODO
  H_i = - apply(pwm, 2, function(col) sum(col * log2(col), na.rm=T))
  e_n = 0
  if(!is.null(Nseqs)) e_n = (1/logb(2)) * (N-1)/(2*Nseqs) 
 
  R_i = log2(N) - (H_i  + e_n)
  # Set any negatives to 0
  R_i = pmax(R_i, 0)
  return(median(R_i))
}

df <- df %>% group_by_(.dots = group_cols) %>% nest()
          
df$data <- lapply(df$data, getPWM)
df$medianic <- lapply(df$data, computeBits)

df <- df %>% select(-data) %>% unnest()

write.csv(df, file = outfile, row.names=F)

