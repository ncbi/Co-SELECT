#!/usr/bin/env Rscript
suppressMessages(require(optparse))
suppressMessages(require(cowplot))
source("seqlogo_scripts/common_funcs.R")
source("results_scripts/plot_pca_func.R")
library(ggthemes)
library(gtable)
library(ggrepel)


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


promiscuous = c('HHHMHH', 'HHHHMH')


getLogoPlot <- function(pwms, seq_length) {
    p <- ggplot() + 
         geom_logo(data = pwms, col_scheme = cs) +
         theme_logo() +
         theme_bw() +
         theme(panel.background = element_rect(fill = NA, color = "black")) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         theme(legend.position='bottom', legend.box = "horizontal") +
         guides(fill=FALSE) + #guide_legend(nrow=1,byrow=TRUE)) +
         #geom_rect(data = pwms, aes(fill = context), xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1) +
         scale_x_continuous(breaks=1:seq_length, lim=c(0.5, seq_length + 0.5)) +
         #geom_label(data = pwms, aes(label = nseq), x=-Inf, y=Inf, hjust = -0.1, vjust = +1.25) +
         theme_tufte() + 
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         #theme(strip.text.y = element_blank()) +
         theme(panel.spacing.y = unit(0.5, "lines")) +
         theme(panel.border = element_blank()) +
         theme(axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")) +
         theme(plot.margin = margin(1, 1, 1, 1, "mm")) +
         scale_y_continuous(breaks=seq(0,2,0.5), lim=c(-0.4, 2))
    p
}



getSeqLogoTF <- function(tf, primer, shapemer) {

opt = list(in_file = paste0('../data/homeodomain/', tf, '/', primer, '/1111.non.all.d0.combined.MGW.S4.000M4.600H5.600X.mer.pwm'),
           title = 'hihi',
           out_file = 'hihi')


path = paste0('../results/d0/publish/1.20/homeodomain/', tf, '/', primer, '/enriched.', tf, '.', primer, '.MGW.4.TAAT.1.1.csv')
print(path)

enrich = read.csv(path, stringsAsFactors=F)
enrich = enrich[,c('kmer', 'label')]
names(enrich) <- c('shapemer', 'label')
enrich2 = enrich
enrich2$shapemer = stringi::stri_reverse(enrich$shapemer)
enrich = unique(rbind(enrich, enrich2))
print(enrich)

pwms = read.csv(opt$in_file, stringsAsFactors=F)
#print(pwms)



all_cols = names(pwms)
mat_cols = c('base', all_cols[grepl("X\\d+", all_cols)])
group_cols = setdiff(all_cols, mat_cols)
seq_length = length(mat_cols) - 1

df <- pwms %>%
          group_by_(.dots = group_cols) %>%
          nest()

df$data <- lapply(df$data, getPwm)


#print(df)

df = df %>% left_join(enrich)
#print(df)

fg0 = df %>% filter(context == 'fg' & cycle == 0)
total = Reduce('+', fg0$data)
#print(total)


fg4enriched <- df %>% filter(context == 'fg' & cycle == 4 & label %in% c('fg', 'both'))
enriched2 = Reduce('+', fg4enriched$data)
#print(enriched2)


df = df %>% filter(label == 'both')
#print(df)

left_panel = tibble(context='fg-total', shapemer = 'AAAAAA', cycle=0, data=list(total)) # %>%
             # bind_rows(tibble(context='fg-enriched', shapemer = 'AAAAAA', cycle=4, data=list(enriched2)))

pwms = df %>% select(-c(motif, label))

pwms$context = factor(pwms$context, levels=c('fg-total', 'fg-enriched', 'fg', 'bg'), ordered=T)
pwms$nseq <- lapply(pwms$data, function(x) max(colSums(x[, 1:seq_length])))
left_panel$nseq <- lapply(left_panel$data, function(x) max(colSums(x[, 1:seq_length])))
left_panel$context = factor(left_panel$context, levels=c('fg-total', 'fg-enriched', 'fg', 'bg'), ordered=T)

pwms = pwms %>% filter(!((context == 'fg') & (cycle == 0)))
pwms = pwms %>% filter(!((context == 'bg') & (cycle %in% c(1,2,3))))

pwms$type <- substr(pwms$context, 1, 2)
pwms = pwms %>% group_by(shapemer, type) %>% tally() %>% group_by(shapemer) %>% summarize(fg = sum(type == 'fg')) %>% filter(fg > 0) %>% inner_join(pwms)

#pdf(paste0(tf, '_', primer, '.pdf'), width=20, height = 2*length(unique(pwms$shapemer)))

#print(names(left_panel))
#print(names(pwms))

pwms <- pwms[, names(left_panel)]

#shapemer = setdiff(unique(pwms$shapemer), promiscuous)[1]

pwms = pwms[pwms$shapemer == shapemer, ]
#print(pwms)

pwms <- rbind(left_panel, pwms)

pwms$shapemer = shapemer

data.frame(pwms, tf=tf)
}

niceLabel <- function(cycle, context) {
  if (context == 'fg-total') {
    part1 = "Motif containing sequences"
  } else if (context == 'fg') {
    part1 = "Motif containing sequences with preferred shapemer"
  } else if (context == 'bg') {
    part1 = "Motif free sequences with preferred shapemer"
  } else {
    part1 = ''
  }
  if (cycle == 0) {
    part2 = 'in initial pool'
  } else if (cycle == 4) {
    part2 = 'in final round'
  } else {
    part2 = ''
  }
  paste(part1, part2)
}

getLogos <- function() {
  pwms1 <- getSeqLogoTF('PITX3', 'TGCATC20NGA', 'HXHHMH')
  #pwms2 <- getSeqLogoTF('PHOX2B', 'TGGTCT20NGA')
  pwms2 <- getSeqLogoTF('HMX2', 'TAGTGG20NCG', 'HXHMMH')
  
  pwms <- rbind(pwms1, pwms2)
  print(pwms)
  pwms$rowlabel = mapply(niceLabel, pwms$cycle, pwms$context)
  print(pwms[, c('context', 'cycle', 'rowlabel')])
  pwms$context = factor(pwms$context, levels=c('fg-total', 'fg-enriched', 'fg', 'bg'), labels=c('motif-containing-total', 'motif-containing-enriched', 'motif-containing', 'motif-free'), ordered=T)
  pwms = pwms[order(pwms$context, pwms$cycle), ]
  pwms$rowlabel = factor(pwms$rowlabel, levels=unique(pwms$rowlabel), ordered=T)
  
  seq_length=10
  
  #print(pwms)
  
  
  shapemers = unique(pwms[, c('tf', 'shapemer', 'context', 'rowlabel')])
  shapemers = shapemers[shapemers$context %in% c('motif-containing', 'motif-free'), ]
  #print(shapemers)
  
  shapemers = shapemers %>% mutate(label = strsplit(as.character(shapemer), split=''), x=strsplit(as.character('3,4,5,6,7,8'), split=','),
                                  color  = strsplit(as.character('S,S,S,D,S,S'), split=',')) %>% unnest()
  #print(shapemers)
  
  shapemers$x = as.integer(shapemers$x)
  
  p2 <- getLogoPlot(pwms, seq_length)
  p2 <- p2 + facet_grid(rowlabel ~ tf, labeller = label_wrap_gen(16)) +
             theme(strip.text.y = element_text(angle = 0, hjust = 0,   vjust = 0.6, size = 12, face='plain')) +
             theme(strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0,   size = 14, face='plain', margin = margin(0,0,0.6,0, "cm")))
  
  
  p2 <- p2 + geom_label(data=shapemers, aes(label=label, x=x, color=color), y = -0.25, size=3) +
             scale_color_manual(guide="none", values=c('S'='black', 'D'='red')) +
             coord_cartesian(clip = 'off')
  return(p2)
}

getNmrPlot <- function() {
  seqmotif = unlist(strsplit('GCTCTAATCCCCG', split=''))
  
  #fileConn<-file("output.txt")
  #writeLines(paste(seqmotif, collapse=''), fileConn)
  #close(fileConn)
  #pred <- getShape("output.txt", shapeType='MGW')$MGW[1,]
  pred = c(NA,NA,5.04,5.41,6.00,5.57,4.81,4.19,4.81,4.76,4.94,NA,NA)
  print(pred)
  
  
  df = data.frame(base=seqmotif, mgw=pred, pos=1:length(seqmotif))
  print(df)
  
  mgw_cutoffs = c(3, 4.0, 4.65, 5.67, 7)
  ypos = (mgw_cutoffs[1:(length(mgw_cutoffs)-1)] + mgw_cutoffs[2:(length(mgw_cutoffs))])/2
  ylabs = c('S', 'M', 'H', 'X')
  lab = data.frame(pos = ypos, label = ylabs)
  
  shapemer = data.frame(x=4:9, label=c('H','X','H','H','M','H'))
  
  p1 <- ggplot(df, aes(x=pos, y=mgw)) +
         geom_point() +
         geom_line() +
         scale_x_continuous(breaks=df$pos, labels=df$base) +
         geom_hline(yintercept = mgw_cutoffs[2:4], linetype = 'dotted') +
         ylim(mgw_cutoffs[1], mgw_cutoffs[5]) +
         geom_text(data=lab, aes(label = label, y=pos), x=1) +
         geom_label(data=shapemer, aes(label = label, x=x), y=3.15, size=3) +
         labs(x="", y=expression(paste("MGW [", ring(A), "]")))
  
  p2 <- ggdraw() + draw_image("2lkx_screenshot.png", scale = 1.0)
  p <- plot_grid(p2, p1, labels = c("A", "B"), ncol=1)
}


pca_csv = '../results/d0/pca_shapemers_publish.cycle4.l1.r1.csv'
pca <- read.csv(pca_csv, stringsAsFactors=F)
pca$label <- ifelse((pca$tf %in% c('PITX3', 'HMX2')) & (pca$barcode %in% c('TGCATC20NGA', 'TAGTGG20NCG')), pca$tf, '')



p1 <- getNmrPlot()
p2 <- getLogos()
p3 <- plot_grid(p1, p2, labels = c('', 'C'), ncol=2, rel_widths = c(2,3.5))

p4 <- getPCAplot(pca)

p <- plot_grid(p3, p4, labels = c('', 'D'), ncol=1, rel_heights = c(3,2.5))

pdf('fig_biological_validation_both.pdf', width=10, height = 10)
print(p)

pdf('fig_biological_validation.pdf', width=11, height = 6)
print(p3)


pdf('fig_pca.pdf', width=8, height = 5)
p5 <- getPCAplot(pca,FALSE)
print(p5)

quit()
