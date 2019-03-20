source("seqlogo_scripts/common_funcs.R")
library(ggthemes)

promiscuous = c('HHHMHH', 'HHHHMH', 'HHHHHH')

pwms = read.csv('../seqlogos/d0/new.seqlogo.MGW.allcycles.1.1.csv', stringsAsFactors=F)
pwms = pwms[pwms$context == 'bg',]
pwms = pwms[pwms$shapemer %in% promiscuous,]
print(head(pwms))


all_cols = names(pwms)
mat_cols = c('base', all_cols[grepl("X\\d+", all_cols)])
group_cols = setdiff(all_cols, mat_cols)
seq_length = length(mat_cols) - 1

df <- pwms %>%
          group_by_(.dots = group_cols) %>%
          nest()

df$data <- lapply(df$data, getPwm)

normalizePwm <- function(m) {
  print(m)
  nseq = max(colSums(m))
  print(nseq)
  m = m / nseq
  print(m)
  m
}

df$data <- lapply(df$data, normalizePwm)

print(df)

df <- df %>% select(shapemer, cycle, tf, motif, family, data)

print(df)

xx = df %>% group_by(shapemer,cycle, family) %>% summarize(num_tf=length(tf))
min_tf = min(xx$num_tf)
print(min_tf)
df = df %>% group_by(shapemer,cycle, family) %>% sample_n(min_tf)

print(df)

df = df %>% group_by(shapemer,cycle) %>% do(data = Reduce("+", .$data))
print(df)

#df %>% unnest() 
shapemers <- unique(df$shapemer)
cycles <- unique(df$cycle)
contexts <- c('bg') #unique(df$context)

df$Round = df$cycle

df$shapemer = factor(df$shapemer, levels=promiscuous, ordered=T)

p <- ggplot() + 
     geom_logo(data = df, col_scheme = cs) +
     theme_logo() +
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
     scale_y_continuous(breaks=seq(0,2,0.5), lim=c(-0.45, 2))

 p <- p + facet_grid(shapemer ~ Round,
                     labeller = labeller(Round = label_both))
p <- p + theme(strip.text.y = element_blank(), strip.text.x = element_text(size=12))

print(names(df))

shapemers = unique(df[, c('shapemer', 'Round')])
print(shapemers)

shapemers = shapemers %>% mutate(label = strsplit(as.character(shapemer), split=''), x=strsplit(as.character('3,4,5,6,7,8'), split=',')) %>% unnest()
print(shapemers)

shapemers$x = as.integer(shapemers$x)

p <- p + geom_label(data=shapemers, aes(label=label, x=x), y = -0.3, size=2.5) + coord_cartesian(clip = 'off')

 #p <- p + ggtitle(opt$title)

height = 1.0 + 0.75*length(shapemers)
width = 1.0 + 2.0*length(contexts)*length(cycles)

print(height)
print(width)

 pdf('fig_seqlogo_promiscuous.pdf', height=height, width=width)
 print(p)

 quit()

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

