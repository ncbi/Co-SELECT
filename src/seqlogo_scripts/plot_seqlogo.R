#!/usr/bin/env Rscript
suppressMessages(require(optparse))
suppressMessages(require(ggthemes))
suppressMessages(require(tidyr))
suppressMessages(require(dplyr))
suppressMessages(require(spseqlogo))
suppressMessages(require(ggplot2))
#source("seqlogo_scripts/common_funcs.R")

tf = 'HMX2'
primer = 'TAGTGG20NCG'

enrich_file = paste0('../results/d0/publish/1.20/homeodomain/', tf, '/', primer, '/enriched.', tf, '.', primer, '.MGW.4.TAAT.1.1.csv')
in_file = paste0('../data/homeodomain/', tf, '/', primer, '/1111.non.all.d0.combined.MGW.S4.000M4.600H5.600X.mer.pwm')

option_list <- list( 
    make_option(c('-i', "--in_file"), default=in_file, 
        help = "Input file. [default \"%default\"]"),
    make_option(c('-e', "--enrich_file"), default=enrich_file, 
        help = "Enrichment file. [default \"%default\"]"),
    make_option(c('-t', "--title"), default=paste(tf, primer, sep = '_'), 
        help = "Output file. [default \"%default\"]"),
    make_option(c('-o', "--out_file"), default="logos.pdf", 
        help = "Output file. [default \"%default\"]")
    )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

print(opt$in_file)
print(opt$enrich_file)
print(opt$out_file)




enrich = read.csv(opt$enrich_file, stringsAsFactors=F)
enrich = enrich[,c('kmer', 'label')]
enrich = enrich[enrich$label == 'both', ]
names(enrich) <- c('shapemer', 'label')
enrich2 = enrich
enrich2$shapemer = stringi::stri_reverse(enrich$shapemer)
enrich = unique(rbind(enrich, enrich2))
print("Enriched shapemers")
print(enrich)

if (all(is.na(enrich))) {
  text = paste0("No enriched shapemers in motif-free sequences for\n", opt$title, ".")
  print(paste(text, "Quiting ..."))
suppressMessages(require(grid))
suppressMessages(require(gridExtra))
  grb=textGrob(text)
  pdf(opt$out_file, title="no logo")
  grid.arrange(grb)
  dev.off()
  quit()
}


cs = make_col_scheme(
       chars = c('A', 'C', 'G', 'T', 'U', 'X', 'B', 'D', 'E', 'F', 'H', 'I', 'J'),
       groups = c('A', 'C', 'G', 'T', 'U', 'full', 'fgfull', 'motif-containing', 'motif-free', 'fg', 'bg', 'motif-containing-total', 'motif-containing-enriched'),
       cols = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839', 'white', 'green', 'blue', 'red', 'blue', 'red', 'blue', 'blue'))

getPwm <- function(df) {
  mat <- as.matrix(df[, !(names(df) %in% c('base'))])
  if (all(mat == 0)) {
    mat <- cbind(mat, c(1,0,0,0))
  }
  row.names(mat) <- df$base
  mat
}


getLogoPlot <- function(pwms) {
  all_cols = names(pwms)
  mat_cols = c('base', all_cols[grepl("X\\d+", all_cols)])
  group_cols = setdiff(all_cols, mat_cols)
  seq_length = length(mat_cols) - 1
  
  pwms <- pwms %>%
            group_by_(.dots = group_cols) %>%
            nest()
  
  pwms$data <- lapply(pwms$data, getPwm)
  pwms$nseq <- lapply(pwms$data, function(x) max(colSums(x[, 1:seq_length])))
  pwms <- pwms %>% unnest(nseq)

  pwms$context <- ifelse(pwms$context == 'fg', 'motif-containing', 'motif-free')

  print(pwms)
  
  p <- ggplot() + 
         geom_logo(data = pwms, col_scheme = cs) +
         theme_logo() +
         theme_bw() +
         theme(panel.background = element_rect(fill = NA, color = "black")) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         theme(legend.position='bottom', legend.box = "horizontal") +
         guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
         #geom_rect(data = pwms, aes(fill = context), xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1) +
         scale_x_continuous(breaks=1:seq_length, lim=c(0.5, seq_length + 0.5)) +
         #geom_label(data = pwms, aes(label = nseq), x=-Inf, y=Inf, hjust = -0.1, vjust = +1.25) +
         theme(axis.line = element_blank()) +
         ylim(0,2)
}



pwms = read.csv(opt$in_file, stringsAsFactors=F)
print(head(pwms))

pwms = merge(pwms, enrich)

shapemers = unique(pwms[, c('shapemer', 'cycle', 'context')])
print("Shapemers present in some cycle and context")
print(shapemers)

shapemers = shapemers %>% group_by(shapemer, cycle) %>% print(n=Inf) %>% summarize(num_ctx = length(context)) %>% print(n=Inf) %>% filter(num_ctx >= 2) %>% select(shapemer) %>% unique() %>% print() %>% unnest()
print(shapemers)
shapemers = unique(shapemers$shapemer)

print(shapemers)


pwms = pwms[pwms$shapemer %in% shapemers, ]

#pwms$enrichment <- ifelse(pwms$enrichment == 'bg', 'motif-free', pwms$enrichment)
#pwms <- with(pwms, pwms[order(enrichment, shapemer), ])
pwms$shapemer <- factor(pwms$shapemer, levels = unique(pwms$shapemer), ordered=T)

shapemers <- unique(pwms$shapemer)
cycles <- unique(pwms$cycle)
contexts <- unique(pwms$context)

pwms$context <- ifelse(pwms$context == 'fg', 'Motif-containing', 'Motif-free')
pwms = pwms[order(pwms$context, pwms$cycle, pwms$shapemer), ]
pwms$header <- paste0(pwms$context, ', Round ', pwms$cycle)
pwms$header <- factor(pwms$header, levels=unique(pwms$header), ordered=T)

 
  shapemers = unique(pwms[, c('shapemer'), drop=F])
  print("Unique shapemers")
print(shapemers)

  shapemers = shapemers %>% group_by(shapemer) %>% mutate(label = strsplit(as.character(shapemer), split='')) %>% print(n=Inf) %>% mutate(x=strsplit(as.character('3,4,5,6,7,8'), split=',')) %>% print() %>%  unnest()
print(shapemers)
  
  shapemers$x = as.integer(shapemers$x)

  print(shapemers)
  
print(pwms)

p <- getLogoPlot(pwms)
p <- p + facet_grid(shapemer ~ header) 
p <- p + ggtitle(opt$title) +
     theme_tufte() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
     theme(strip.text.y = element_blank()) +
     theme(strip.text.x = element_text(size=14)) +
     theme(panel.spacing.y = unit(0.5, "lines")) +
     theme(panel.border = element_blank()) +
     theme(plot.margin = margin(1, 1, 1, 1, "mm")) +
         guides(fill=FALSE)

 
 
  p <- p + geom_label(data=shapemers, aes(label=label, x=x), y = -0.3, size=3) +
         scale_y_continuous(breaks=seq(0,2,0.5), lim=c(-0.4, 2)) +
             coord_cartesian(clip = 'off')

height = 1 + 3.5*length(shapemers)
width = 1 + 2.5*length(contexts)*length(cycles)

print(opt$out_file)

pdf(opt$out_file, height=height, width=width)

print(p)

