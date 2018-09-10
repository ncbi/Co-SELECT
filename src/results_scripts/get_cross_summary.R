#!/usr/bin/env Rscript
require(plyr, quietly=TRUE)
source("results_scripts/tf_utils.R")
#source("results_scripts/params_config.R")

top_data_dir <- '../data'
top_res_dir <- '../results'
#SPDEF,TCATTG20NCG,ETS,,1;2;3;4,Yes,GGAT;GNAA,GGAA,1,1
#ARNTL,TCAAAA20NCG,bHLH,,1;2;3;4,No,,CACGTG,1,1

tf1 <- 'SPDEF'
primer1 <- 'TCATTG20NCG'
motif1 <- 'GGAA'


tf2 <- 'ARNTL'
primer2 <- 'TCAAAA20NCG'
motif2 <- 'CACGTG'

#232  ALX3  TGCAAG20NGA  homeodomain       NaN  1;2;3;4        No   NaN  TAAT   
#17   ELF1  TATGGG20NCG          ETS       NaN  1;2;3;4        No   NaN  GGAA 

tf1 <- 'ALX3'
primer1 <- 'TGCAAG20NGA'
motif1 <- 'TAAT'


tf2 <- 'ELF1'
primer2 <- 'TATGGG20NCG'
motif2 <- 'GGAA'


lflank <- 1
rflank <- 1
cycle <- 4
shape <- 'MGW'
shape_levels <- 'haha'

# <tf1> <primer1> <motif1> <tf2> <primer2> <motif2> <cycle> <shape> <lflank> <rflank> 

args = commandArgs(trailingOnly=TRUE)

tf1 = args[1]
primer1 = args[2]
family1 = args[3]
motif1 = args[4]
dist1 = args[5]
cycle1 = as.numeric(args[6])

tf2 = args[7]
primer2 = args[8]
family2 = args[9]
motif2 = args[10]
dist2 = args[11]
cycle2 = as.numeric(args[12])

shape = args[13]
lflank = as.numeric(args[14])
rflank = as.numeric(args[15])

en_th = args[16]
top_res_dir = args[17]
shape_levels = args[18]

ENRICHMENT_THRESHOLD <<- as.numeric(en_th)
  

makeBoth <- function(x) {
  if (length(x) == 1) {return(x[1])}
  return('both')
}

dir.create(top_res_dir, showWarnings = FALSE, recursive = TRUE)

#tfs <- read.csv('../inventory/tf_use.csv', comment.char='#', stringsAsFactors=FALSE)
#tfs <- tfs[tfs$run == 1,]
#
#tmp <- tfs[, c('family', 'tf', 'primer', 'motif')]
#
#cross <- merge(tmp, tmp, by=NULL)
#print(head(cross))
#print(nrow(tfs))
#print(nrow(cross))
#
#cross <- cross[cross$family.x != cross$family.y, ]
#print(head(cross))
#print(nrow(tfs))
#print(nrow(cross))
#
#
##cross <- cross[cross$tf.x < cross$tf.y, ]
##print(head(cross))
##print(nrow(tfs))
##print(nrow(cross))

getCrossSummary(top_data_dir, top_res_dir, tf1, primer1, family1, motif1, dist1, cycle1, tf2, primer2, family2, motif2, dist2, cycle2, shape, lflank, rflank, shape_levels)
