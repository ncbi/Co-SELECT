#!/usr/bin/env Rscript
require(plyr)
require(data.table)
require(dplyr)
require(optparse)
#source("tf_utils.R")


en_th <- "1.20"
shape <- 'MGW'

option_list <- list( 
    make_option(c('-s', "--shape"), default=shape, 
        help = "Shape feature to use [default \"%default\"]"),
    make_option(c('-t', "--enrich_threshold"), default=en_th, 
        help = "Enrichment threshold '1.10', '1.20' etc. [default \"%default\"]"),
    make_option(c('-i', "--in_file"), default="denriched.csv", 
        help = "Input file. [default \"%default\"]"),
    make_option(c('-n', "--tfs_per_family_file"), default="tfs_per_family.csv", 
        help = "File showing number of TFs in family. [default \"%default\"]"),
    make_option(c('-o', "--out_file"), default="promiscuous_shapemers.csv", 
        help = "Output file. [default \"%default\"]")
    )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

shape = opt$shape
en_th = opt$enrich_threshold
infile = opt$in_file
outfile = opt$out_file
tfs_per_family_file = opt$tfs_per_family_file

print(shape)
print(en_th)
print(infile)
print(outfile)


computePromiscuity <- function(df, tf_per_family) {
  df <- df[,c('tf.x', 'primer.x', 'family.x', 'kmer', 'label')]
  print(head(df))
  names(df) <- c('tf', 'primer', 'family', 'kmer', 'ctx')
  print(head(df))
  df$ctx[df$ctx=='bg'] <- 'distant'
  df$ctx[df$ctx=='both'] <- 'distant'
  df <- df[df$ctx == 'distant',]
  print(head(df))

  df <- ddply(df, .(kmer, family), summarize, count=length(family))
  print(df)

  df <- merge(df, tf_per_family, all.x=TRUE)
  print(df)

  df$frac <- df$count/df$num_tfs
  df <- ddply(df, .(kmer), summarize, num_family=length(family), promiscuity=min(frac))
  # discard those shapemers which are not present in all considered families
  df <- df[df$num_family >= nrow(tf_per_family),]
  df <- df[order(-df$promiscuity), c('kmer', 'promiscuity')]

  return(df)
}

all <- read.csv(infile, stringsAsFactors = F)
tfs_per_family <- read.csv(tfs_per_family_file, stringsAsFactors=F)
print(tfs_per_family)


df <- computePromiscuity(all, tfs_per_family)

print(head(df, n=20))

write.csv(df, file=outfile, row.names=FALSE)
