#!/usr/bin/env Rscript
require(plyr, quietly=TRUE)
source("results_scripts/tf_utils.R")

top_data_dir <- '.'
top_res_dir <- '.'
#SPDEF,TCATTG20NCG,ETS,,1;2;3;4,Yes,GGAT;GNAA,GGAA,1,1
#ARNTL,TCAAAA20NCG,bHLH,,1;2;3;4,No,,CACGTG,1,1

tf <- 'ARNTL'
primer <- 'TCAAAA20NCG'
motif <- 'CACGTG'

tf <- 'SPDEF'
primer <- 'TCATTG20NCG'
motif <- 'GGAA'


tf <- 'MAX'
primer <- 'TGACCT20NGA'
motif <- 'CACGTG'

lflank <- 1
rflank <- 1
cycle <- 4
shape <- 'MGW'
dist <- 2

tf <- 'ALX3'
primer <- 'TGCAAG20NGA'
motif <- 'TAAT'

tf <- 'USF1'
primer <- 'TGACGA20NGCA'
motif <- 'CACGTG'
dist <- 3

tf <- 'Elk3'
primer <- 'TGAGTG20NTGA'
motif <- 'GGAA'
dist <- 2

en_th = '1.20'

shape_levels = 'S4.250M5.400H5.700X'

# <tf> <primer> <cycle> <shape> <motif> <dist> <lflank> <rflank> <top-res-dir> <shape-levels>

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  tf = args[1]
}
if (length(args) > 1) {
  primer = args[2]
}
if (length(args) > 2) {
  cycle = as.numeric(args[3])
}
if (length(args) > 3) {
  shape = args[4]
}
if (length(args) > 4) {
  motif = args[5]
}
if (length(args) > 5) {
  dist = args[6]
}
if (length(args) > 6) {
  lflank = as.numeric(args[7])
}
if (length(args) > 7) {
  rflank = as.numeric(args[8])
}
if (length(args) > 8) {
  en_th = args[9]
}
if (length(args) > 9) {
  top_res_dir =args[10]
}
if (length(args) > 10) {
  shape_levels =args[11]
}
if (length(args) > 11) {
  family =args[12]
}
if (length(args) > 12) {
  top_data_dir =args[13]
}

ENRICHMENT_THRESHOLD <<- as.numeric(en_th)


makeBoth <- function(x) {
  if (length(unique(x)) == 1) {return(x[1])}
  return('both')
}


makeBoth2 <- function(x, y) {
  t <- na.omit(c(x,y,'both'))
  return(na.omit(c(t[3],t[1]))[1])
}

getCrossSummary(top_data_dir, top_res_dir, tf, primer, family, motif, dist, cycle, tf, primer, family, motif, dist, cycle, shape, lflank, rflank, shape_levels)

