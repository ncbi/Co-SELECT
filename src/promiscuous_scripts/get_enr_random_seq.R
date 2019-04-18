#!/usr/bin/env Rscript
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressPackageStartupMessages(library("optparse"))


source("results_scripts/tf_utils.R")
#source("params_config.R")

top_res_dir <- '/panfs/pan1/dnashape/FinalResults/d0/publish/1.20'

lflank <- 1
rflank <- 1
cycle <- 4
EN_TH <- 1.2
num_iter <- 2
shp = 'MGW'

map_file <- fread(paste0('/panfs/pan1/aptax/seqmer2shapemer/map_seq2shape_', shp, '.csv'), stringsAsFactors=FALSE)
print(map_file)

promiscuous_file = paste0(top_res_dir, "/promiscuous_same.", shp,'.', cycle, ".l1.r1.csv") 
out_file <- paste(top_res_dir, paste0("promiscuous_significance_", shp, ".csv"), sep='/')
bg_file <- paste(top_res_dir, paste0("bg_summary_", shp, ".csv"), sep='/')

option_list <- list( 
    make_option(c("-n", "--num_iter"), type="integer", default=5, 
        help="Number of random normals to generate [default %default]",
        metavar="number"),
    make_option(c('-m', "--map_file"), default=map_file, 
        help = "Shapemer to sequence mer mapping file [default \"%default\"]"),
    make_option(c('-p', "--promiscuous_file"), default=promiscuous_file, 
        help = "Promiscuous shapemer file [default \"%default\"]"),
    make_option(c('-b', "--bg_file"), default=promiscuous_file, 
        help = "BG shapemers enrichment info file [default \"%default\"]"),
    make_option(c('-o', "--out_file"), default=out_file, 
        help = "Output file [default \"%default\"]")
    )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

map_file = opt$map_file
num_iter = opt$num_iter
map_file = opt$map_file
promiscuous_file = opt$promiscuous_file
bg_file = opt$bg_file
out_file = opt$out_file


all_promiscuous <- fread(promiscuous_file)
print(all_promiscuous)

#all_promiscuous <- all_promiscuous[promiscuous == TRUE, .N, keyby = .(shape, kmer)]
#all_promiscuous <- all_promiscuous[, minfrac, keyby = .(shape, kmer)]
#print(all_promiscuous)

num_promis <- 2
num_promis <- nrow(all_promiscuous)


seq2shape = read.csv(map_file, stringsAsFactors=F)
num_sequences <- merge(all_promiscuous, seq2shape, by.x='kmer', by.y='shapemer')[,.(numseq=.N), by=.(kmer)]
print(num_sequences)


tfs <- fread(bg_file)
print(tfs)

computePromiscuity <- function(dt, new_col_name) {
  #print(dt)
  dt <- dt[, list(count=sum(enrichment >= EN_TH, na.rm=TRUE), frac=mean(enrichment >= EN_TH, na.rm=TRUE)), by=list(kmer, family)]
  dt <- dt[, list(minfrac=min(frac)), by=list(kmer)]
  setnames(dt, 'minfrac', new_col_name)
  return(dt)
}


results <- all_promiscuous[, .(promiscuous = promiscuity), by=kmer]

for (iter in 1:num_iter) {
  print(paste('iter', iter))
  new_results <- data.table()
  for (promiscuous in all_promiscuous$kmer[1:num_promis]) {
    skip_shapemers <- c(promiscuous)
    rdf <- seq2shape[!(shapemer %in% skip_shapemers)]
    num_seq <- seq2shape[shapemer == promiscuous, .N]
    new_random <- rdf[sample(.N, num_seq),list(count=.N),keyby=shapemer]
    new_random <- merge(new_random, tfs, by.x='shapemer', by.y='kmer')
    new_random <- new_random[, .(enrichment = sum(bg4.frac, na.rm=TRUE)/sum(bg0.estimated, na.rm=TRUE), kmer=promiscuous), by=list(tf, primer,family)]
    new_results <- rbind(new_results, new_random)
  }
  new_col <- paste('en.rand', iter, sep='.')
  results <- merge(results, computePromiscuity(new_results, new_col))
}
results <- melt(results, id.vars=c('kmer', 'promiscuous'))
results <- results[, list(promiscuity=min(promiscuous), pvalue = mean(value >= promiscuous)), by=list(kmer)]
results <- merge(results, num_sequences, by='kmer')
setorder(results, -promiscuity, kmer)
print(results)
fwrite(results, out_file)
