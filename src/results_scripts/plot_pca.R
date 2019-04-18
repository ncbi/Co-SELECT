#!/usr/bin/env Rscript
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
script.dir <- dirname(thisFile())
src_file = paste(script.dir, "plot_pca_func.R", sep='/')
print(src_file)
source(src_file)

args = commandArgs(trailingOnly=TRUE)

pca <- read.csv(args[1], stringsAsFactors=F)
pca$label <- ifelse(pca$tf == 'ISX', paste(pca$tf, pca$barcode, sep='.'), pca$tf)

p <- getPCAplot(pca)

pdf(args[2], width=9, height=6)

print(p)
