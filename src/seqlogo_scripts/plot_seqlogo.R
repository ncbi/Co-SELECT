#!/usr/bin/env Rscript
suppressMessages(require(optparse))
source("seqlogo_scripts/common_funcs.R")


option_list <- list( 
    make_option(c('-i', "--in_file"), default="SEQLOGO/seqlogo_info_PITX3_TGCATC20NGA_NEW.csv", 
        help = "Input file. [default \"%default\"]"),
    make_option(c('-t', "--title"), default="PITX3_TGCATC20NGA", 
        help = "Output file. [default \"%default\"]"),
    make_option(c('-o', "--out_file"), default="logos.pdf", 
        help = "Output file. [default \"%default\"]")
    )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

print(opt$in_file)

pwms = read.csv(opt$in_file, stringsAsFactors=F)
pwms = pwms[pwms$enrichment != "", ]
pwms$promiscuity <- NULL
if (all(is.na(pwms))) {
  text = paste0("No enriched shapemers in motif-free sequences for\n", opt$title, ".")
  print(paste(text, "Quiting ..."))
  grb=textGrob(text)
  pdf(opt$out_file)
  grid.arrange(grb)
  dev.off()
  quit()
}

pwms$enrichment <- ifelse(pwms$enrichment == 'bg', 'motif-free', pwms$enrichment)
pwms <- with(pwms, pwms[order(enrichment, shapemer), ])
pwms$shapemer <- factor(pwms$shapemer, levels = unique(pwms$shapemer), ordered=T)

shapemers <- unique(pwms$shapemer)
cycles <- unique(pwms$cycle)
contexts <- unique(pwms$context)

p <- getLogoPlot(pwms)
p <- p + facet_grid(shapemer + enrichment ~ context + cycle,
                    labeller = labeller(enrichment = label_both, cycle = label_both))
p <- p + ggtitle(opt$title)

height = 1 + 2.5*length(shapemers)
width = 1 + 3*length(contexts)*length(cycles)

pdf(opt$out_file, height=height, width=width)

print(p)

