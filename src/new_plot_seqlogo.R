#!/usr/bin/env Rscript
suppressMessages(require(optparse))
suppressMessages(require(cowplot))
source("seqlogo_scripts/common_funcs.R")
library(ggthemes)


#option_list <- list( 
#    make_option(c('-i', "--in_file"), default="SEQLOGO/seqlogo_info_PITX3_TGCATC20NGA_NEW.csv", 
#        help = "Input file. [default \"%default\"]"),
#    make_option(c('-t', "--title"), default="PITX3_TGCATC20NGA", 
#        help = "Output file. [default \"%default\"]"),
#    make_option(c('-o', "--out_file"), default="logos.pdf", 
#        help = "Output file. [default \"%default\"]")
#    )
#
## get command line options, if help option encountered print help and exit,
## otherwise if options not found on command line then set defaults, 
#opt <- parse_args(OptionParser(option_list=option_list))
#
#print(opt$in_file)
#print(opt$out_file)


tf = 'PITX3'
primer = 'TGCATC20NGA'


tf = 'PHOX2B'
primer = 'TGGTCT20NGA'


getLogoPlot <- function(pwms, seq_length) {
    p <- ggplot() + 
         geom_logo(data = pwms, col_scheme = cs) +
         theme_logo() +
         theme_bw() +
         theme(panel.background = element_rect(fill = NA, color = "black")) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         theme(legend.position='bottom', legend.box = "horizontal") +
         guides(fill=FALSE) + #guide_legend(nrow=1,byrow=TRUE)) +
         geom_rect(data = pwms, aes(fill = context), xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1) +
         scale_x_continuous(breaks=1:seq_length, lim=c(0.5, seq_length + 0.5)) +
         #geom_label(data = pwms, aes(label = nseq), x=-Inf, y=Inf, hjust = -0.1, vjust = +1.25) +
     theme_tufte() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
     #theme(strip.text.y = element_blank()) +
     theme(panel.spacing.y = unit(0.5, "lines")) +
     theme(panel.border = element_blank()) +
     theme(axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")) +
     theme(plot.margin = margin(1, 1, 1, 1, "mm")) +
         ylim(0,2)
    p
}



getSeqLogoTF <- function(tf, primer) {

opt = list(in_file = paste0('../data/homeodomain/', tf, '/', primer, '/1111.non.all.d0.combined.MGW.S4.000M4.600H5.600X.mer.pwm'),
           title = 'hihi',
           out_file = 'hihi')



enrich = read.csv(paste0('../results/d0/publish/1.20/homeodomain/', tf, '/', primer, '/enriched.', tf, '.', primer, '.MGW.4.TAAT.1.1.csv'), stringsAsFactors=F)
enrich = enrich[,c('kmer', 'label')]
names(enrich) <- c('shapemer', 'label')
enrich2 = enrich
enrich2$shapemer = stringi::stri_reverse(enrich$shapemer)
enrich = unique(rbind(enrich, enrich2))
print(enrich)

pwms = read.csv(opt$in_file, stringsAsFactors=F)
print(pwms)

all_cols = names(pwms)
mat_cols = c('base', all_cols[grepl("X\\d+", all_cols)])
group_cols = setdiff(all_cols, mat_cols)
seq_length = length(mat_cols) - 1

df <- pwms %>%
          group_by_(.dots = group_cols) %>%
          nest()

df$data <- lapply(df$data, getPwm)


print(df)

df = df %>% left_join(enrich)
print(df)

fg0 = df %>% filter(context == 'fg' & cycle == 0)
total = Reduce('+', fg0$data)
print(total)


fg4enriched <- df %>% filter(context == 'fg' & cycle == 4 & label %in% c('fg', 'both'))
enriched2 = Reduce('+', fg4enriched$data)
print(enriched2)


df = df %>% filter(label == 'both')
print(df)

left_panel = tibble(context='fg-total', shapemer = 'AAAAAA', cycle=0, data=list(total)) %>%
              bind_rows(tibble(context='fg-enriched', shapemer = 'AAAAAA', cycle=4, data=list(enriched2)))

pwms = df %>% select(-c(motif, label))

pwms$context = factor(pwms$context, levels=c('fg-total', 'fg-enriched', 'fg', 'bg'), ordered=T)
pwms$nseq <- lapply(pwms$data, function(x) max(colSums(x[, 1:seq_length])))
left_panel$nseq <- lapply(left_panel$data, function(x) max(colSums(x[, 1:seq_length])))
left_panel$context = factor(left_panel$context, levels=c('fg-total', 'fg-enriched', 'fg', 'bg'), ordered=T)

pwms = pwms %>% filter(!((context == 'fg') & (cycle == 0)))

pwms$type <- substr(pwms$context, 1, 2)
pwms = pwms %>% group_by(shapemer, type) %>% tally() %>% group_by(shapemer) %>% summarize(fg = sum(type == 'fg')) %>% filter(fg > 0) %>% inner_join(pwms)

#pdf(paste0(tf, '_', primer, '.pdf'), width=20, height = 2*length(unique(pwms$shapemer)))

p1 <- getLogoPlot(left_panel, seq_length)
p2 <- getLogoPlot(pwms, seq_length)

p1 <- p1 + facet_grid( ~ context + cycle,
                    labeller = labeller(enrichment = label_both, cycle = label_both))

p2 <- p2 + facet_grid(shapemer ~ context + cycle,
                    labeller = labeller(enrichment = label_both, cycle = label_both))

p <- plot_grid(p1, p2, rel_widths=c(2,6))
}


p1 <- getSeqLogoTF('PITX3', 'TGCATC20NGA')
p2 <- getSeqLogoTF('PHOX2B', 'TGGTCT20NGA')


p <- plot_grid(p1, p2, rel_heights=c(1+2,1+1), labels = c("PITX3", "PHOX2B"), ncol=1)

pdf('seqlogo_pitx3_phox2b.pdf', width=16, height = 6.5)

print(p)

quit()
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

print(pwms)

p <- getLogoPlot(pwms)
p <- p + facet_grid(shapemer + enrichment ~ context + cycle,
                    labeller = labeller(enrichment = label_both, cycle = label_both))
p <- p + ggtitle(opt$title)

height = 1 + 2.5*length(shapemers)
width = 1 + 3*length(contexts)*length(cycles)

print(opt$out_file)

pdf(opt$out_file, height=height, width=width)

print(p)

