#!/usr/bin/env Rscript
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressPackageStartupMessages(library("optparse"))

source("results_scripts/tf_utils.R")
#source("params_config.R")

top_data_dir <- '../data'
top_prob_dir <- '../data/simple_predict'
top_res_dir <- '../results/d0/publish/1.20'

lflank <- 1
rflank <- 1
cycle <- 4
shape_type <- 'MGW'

outfile <- paste(top_res_dir, paste0("bg_summary_", shape_type, ".csv"), sep='/')
print(outfile)

option_list <- list( 
    make_option(c('-s', "--shape"), default="MGW", 
        help = "Shape feature to use [default \"%default\"]"),
    make_option(c('-o', "--outfile"), default=outfile, 
        help = "Output file [default \"%default\"]")
    )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

shape_type = opt$shape
outfile = opt$outfile

print(shape_type)
print(outfile)

shape_levels <- fread('shape_levels.csv')[shape == shape_type & levels_type == 'publish', levels_str]
print(shape_levels)


getTfs <- function() {
  tfs <- read.csv('tf_run_coselect.csv', comment.char='#', stringsAsFactors=FALSE)
  details <- read.csv('tf_inventory_jolma_ronshamir.csv', comment.char='#', stringsAsFactors=FALSE)
  motifs <- read.csv('tf_coremotif.csv', comment.char='#', stringsAsFactors=FALSE)
  tfs = merge(tfs, details)
  tfs = merge(tfs, motifs)
  #tfs <- tfs[tfs$run == 1,]
  tfs$distance = nchar(tfs$motif) - 2
  return(tfs)
}

getData <-function(top_data_dir, top_prob_dir, dt, cycle, shape_type, lflank, rflank, shape_levels) {

  print(paste(dt[1,tf], dt[1,primer]))

  df <- getEnrichmentInfo(top_data_dir, top_prob_dir, dt[1,tf], dt[1,primer], dt[1,family], cycle, 'bg', shape_type, dt[1,motif], dt[1,distance], lflank, rflank, shape_levels) 
  df <- df[, c('kmer', 'bg4.frac', 'bg0.estimated')]
  df$tf <- dt[1,tf]
  df$primer <- dt[1,primer]
  df$family <- dt[1,family]
  
  return(as.list(df))

}

tfs <- data.table(getTfs() %>% select(tf, primer, family, motif, distance), keep.rownames=TRUE)

tfs <- rbindlist(lapply(1:nrow(tfs), function(x) { getData(top_data_dir, top_prob_dir, tfs[x], cycle, shape_type, lflank, rflank, shape_levels)}))

print(head(tfs))
print(nrow(tfs))


fwrite(tfs, outfile)

