#!/usr/bin/env Rscript
suppressMessages(require(optparse))
source("seqlogo_scripts/common_funcs.R")


option_list <- list( 
    make_option(c('-i', "--in_file"), default="SEQLOGO/seqlogo_info_PITX3_TGCATC20NGA_NEW.csv", 
        help = "Input file. [default \"%default\"]"),
    make_option(c('-o', "--out_file"), default="logos.pdf", 
        help = "Output file. [default \"%default\"]")
    )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

print(opt$in_file)

pwms = read.csv(opt$in_file, stringsAsFactors=F)
pwms = pwms[pwms$promiscuity == "high", ]
pwms$promiscuity <- NULL
pwms$enrichment <- NULL

res <- dlply(pwms, c("shapemer"))
files <- lapply(names(res), function(shapemer){
 df <- res[[shapemer]]
 if (all(is.na(df))) {
   text = paste0("No enriched shapemers in motif-free sequences for\n", opt$title, ".")
   print(paste(text, "Quiting ..."))
   grb=textGrob(text)
   pdf(opt$out_file)
   grid.arrange(grb)
   dev.off()
   quit()
 }
 
 experiments <- unique(paste(df$tf, df$primer, sep='.'))
 cycles <- unique(df$cycle)
 contexts <- unique(df$context)
 
 p <- getLogoPlot(df)
 p <- p + facet_grid(family + tf + primer + motif ~ context + cycle,
                     labeller = labeller(cycle = label_both))
 p <- p + ggtitle(opt$title)
 
 height = 1 + 2.5*length(experiments)
 width = 1 + 3*length(contexts)*length(cycles)
 
 out_file = paste(tools::file_path_sans_ext(opt$out_file), shapemer, 'pdf', sep='.') 
 pdf(out_file, height=height, width=width)
 print(p)

 out_file
})

df = do.call(rbind, files)
rownames(df) = names(res)
colnames(df) = c("pdf")
print(df)

write.csv(df, file = opt$out_file)

