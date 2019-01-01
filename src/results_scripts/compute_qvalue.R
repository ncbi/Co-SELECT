#!/usr/bin/env Rscript
library(qvalue)
#source("results_scripts/tf_utils.R")

getQvalues <- function(df, text) {
  #qobj <- qvalue(p = df$pvalue, lambda = seq(0.0125, 0.10, 0.0125), pi0.method="bootstrap")
  #qobj <- qvalue(p = df$pvalue, pi0.method="bootstrap")
  #qobj <- qvalue(p = df$pvalue, pi0.method="bootstrap")
  #qobj <- qvalue(p = df$pvalue)
  #print(qobj)
  #qobj <- tryCatch(qvalue(p = df$pvalue, pi0.method="bootstrap"), error=function(e) { NA })
  qobj <- tryCatch(qvalue(p = df$pvalue), error=function(e) { NA })
  #summary(qobj)
  #  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  #  text(x = 0.5, y = 0.5, text, cex = 1.6, col = "black")
  #print(plot(qobj, title='haha'))
  #print(qobj)
  qvalue <- rep(NA, nrow(df))
  pi0 <- 1.0
  if (all(!is.na(qobj))) {
    pi0 <- qobj$pi0           # TODO: HACK for the case qvalue complains
    qvalue = qobj$qvalues
  } else if (max(df$pvalue) <= 0.75) {
    pi0 <- 0.0
    qvalue <- df$pvalue
  }
  return(data.frame(df, qvalue=qvalue, nullratio=pi0))
}


args = commandArgs(trailingOnly=TRUE)

shape = args[1]
family = args[2]
infile = args[3]
outfile = args[4]

sdf <- read.csv(infile, stringsAsFactors=F)
sdf <- sdf[sdf$family.x == family, ]
df <- getQvalues(sdf[, c('tf.x', 'primer.x', 'tf.y', 'primer.y', 'pvalue')])
df$rvalue <- p.adjust(df$pvalue, method='fdr')

write.csv(df, file=outfile, row.names=FALSE)
