#!/usr/bin/env Rscript
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(grid)
require(gtable)
library(ggrepel)
library(png)
source("results_scripts/tf_utils.R")
source("results_scripts/color_utils.R")


mytheme <- theme(legend.position="top",
                 legend.margin=margin(1,1,1,1),
                 legend.justification="right",
                 legend.box.margin=margin(1,1,-07,1),
                 legend.background = element_rect(colour = 'gray', fill = '#FEFEF1'),
                 plot.background = element_rect(fill = "transparent",colour = NA),
                 #plot.margin = margin(2,1,1,1),
                 #panel.spacing.y=unit(2, "lines"),
                 strip.background=element_blank())


getHisto <- function(qf) {
  for_map <- category_abs_colors[unlist(unique(qf$ctx))]
  rev_map <- setNames(names(for_map), unname(for_map))
  histo <- ggplot(qf, aes(x=pvalue, fill=ctx)) + 
             geom_histogram(aes(y=..density..), binwidth=.05, position='dodge')+
             scale_fill_manual(values=for_map)
  df <- ggplot_build(histo)$data[[1]]
  df <- mutate(df, ctx=rev_map[fill])
  df[,c('x', 'ymax', 'ctx')]
}

getOnlyHistoInfo <- function(qf) {
  qf %>%
    group_by(shape,levels,en_th,family,tag) %>%
    do(getHisto(.))
}

getECR <- function(v1, v2) {
  if (v2 < 1e-3) {
    return(1000.0)
  }
  return(round(v1/v2, 3))

}

getRatio <- function(df) {
  df1 <- df[df$ctx=='experiment',]
  df2 <- df[df$ctx=='control',]
  ratio <- getECR(df1$ymax[1], df2$ymax[1])
  ratio2 <- getECR(sum(df1$ymax), sum(df2$ymax))
  label <- paste0('ECR[', deparse("0.05"), ']==',deparse(sprintf("%6.1f", ratio)))
  label2 <- paste0('ECR[', deparse("0.10"), ']==',deparse(sprintf("%6.1f", ratio2)))
  quality <- ifelse(ratio>1.1, 'excellent', ifelse(ratio>1, 'good', 'bad'))
  support <- ifelse(ratio>1, 'yes', 'no')
  specific <- paste('Specific: ', ifelse(ratio>1, '\u2714', '\u2718'))
  #print(specific)
  data.frame(ratio=ratio,ratio2=ratio2,label=label,label2=label2,quality=quality,support=support, specific=specific)
}

interpolateQvalue <- function(df, p) {
  left <- rbind(df[df$pvalue<=p, 'qvalue'],data.frame(qvalue=-1.0))
  right <- rbind(df[df$pvalue>=p, 'qvalue'],data.frame(qvalue=1.0))
  return(min(abs(max(left)),min(right)))
}


getQvalueFDR <- function(df) {
  df1 <- df[df$ctx=='experiment',]
  print('df1')
  print(df1)
  if (nrow(df1) < 1) {
    return(data.frame(ratio=NA,ratio=NA,label='NA',label2='NA'))
  }
  print(df1[df1$pvalue<=0.05, 'qvalue'])
  #print(max(df1[df1$pvalue<=0.05, 'qvalue']))
  #ratio <- abs(max(rbind(df1[df1$pvalue<=0.05, 'qvalue'],data.frame(qvalue=-1.0))))
  #print(ratio)
  #ratio2 <- abs(max(rbind(df1[df1$pvalue<=0.10, 'qvalue'],data.frame(qvalue=-1.0))))
  #print(ratio2)
  ratio <- interpolateQvalue(df1, 0.05)
  print(ratio)
  ratio2 <- interpolateQvalue(df1, 0.10)
  print(ratio2)
  #ratio <- round(min(ratio, min(df1[df1$pvalue>=0.05, 'qvalue'])),3)
  #ratio2 <- round(min(ratio2, min(df1[df1$pvalue>=0.10, 'qvalue'])),3)
  ratio <- round(ratio,3)
  ratio2 <- round(ratio2,3)
  #print(paste("Ratio", ratio))
  #print(paste("Ratio2", ratio2))
  label <- paste0('qvalue[', deparse("0.05"), ']==',deparse(sprintf("%6.3f", ratio)))
  print(label)
  label2 <- paste0('qvalue[', deparse("0.10"), ']==',deparse(sprintf("%6.3f", ratio2)))
  print(label2)
  selection <- paste('Selection: ', ifelse(ratio<=0.1, '\u2714', '\u2718'))
  print(selection)
  return(data.frame(ratio3=ratio,ratio4=ratio2,label3=label,label4=label2, selection=selection))
}

my_labeller <- label_bquote(
  rows = paste("Enrichment threshold:", .(en_th)), #.(am) / alpha,
  cols = paste("Shape", .(shape), "Levels=", .(levels)) #.(vs) ^ .(cyl)
)

enrichThreshold <- function(strings) {
    paste("Enrichment threshold:", strings)
}
labelLevels <- function(strings) {
    strings <- strsplit(strings, "[|]")
    vapply(strings, function(x) {
        paste(x, collapse = sprintf('\u2264'))
    }, FUN.VALUE = character(1))
}

plotQvalueHisto <- function(qf, title, facets, more_details=TRUE) {
  pi00 <- qf %>% group_by(shape,family,ctx) %>%
                 mutate(label=round(max(nullratio), 3))
  df <- getOnlyHistoInfo(qf)
  left_bins <- df[df$x < 0.075,]
  print(left_bins)
  left_bins$label <- format(round(left_bins$ymax, 1))
  print(left_bins)
  left_bins$label.x <- 0.25+6*left_bins$x
  left_bins$label.y <- 0.65*max(left_bins$ymax)
  tf <- left_bins %>%
    group_by(shape,family, en_th, levels) %>%
    do(getRatio(.))
  print(qf)
  print(title)
  xf <- qf %>%
    group_by(shape,family, en_th, levels) %>%
    do(getQvalueFDR(.))
  qf <- merge(qf, tf, all.x=TRUE)
  tf <- merge(tf, xf)
  tf$specific[tf$ratio3 > 0.1] <- NA
  histo <- ggplot(qf, aes(x=pvalue)) 
  #if (more_details) {
  #  histo <- histo + 
  #        geom_segment(data=left_bins, aes(x = x, y = ymax, xend = label.x, yend = label.y, colour = ctx)) +
  #        geom_label(data=left_bins, aes(x=label.x, y=label.y, label=label, color=ctx))
  #}
  only_experiment <- FALSE #unique(qf$ctx)[1] == 'experiment'
  print(only_experiment)
  histo <- histo +
          geom_histogram(aes(y=..density.., fill=ctx, support=support), binwidth=.05, position='dodge') +
          facet_grid(facets, labeller = labeller(en_th=enrichThreshold, levels=labelLevels, .multi_line = FALSE)) + 
          geom_text(data=xf, aes(label=label3), x=-Inf, y=Inf,hjust=-0.25,vjust=1.5, parse=TRUE)
  if (!only_experiment) {
    histo <- histo +
          geom_text(data=tf, aes(label=label), x=-Inf, y=Inf,hjust=-0.5,vjust=2.5, parse=TRUE) 
#          geom_label(data=tf, inherit.aes=FALSE, aes(label=specific),  x= Inf, y=Inf, hjust=1.00, vjust=-0.12, parse=FALSE) 
##  } else {
##  histo <- histo +
##          geom_line(data=pi00, aes(x = pvalue, y = nullratio), color = 'blue', linetype = "solid", size = 0.5) +
##          geom_text(data=pi00, aes(label = paste("hat(pi)[0] ==", label), x = 0.65, y = nullratio+2.5), parse = TRUE, colour = "black" ) 
  }
  if (more_details) {
    histo <- histo + 
#         geom_text(data=tf, aes(label=label2), x=Inf, y=Inf,hjust=1.3,vjust=2.75, parse=TRUE) +
#          geom_text(data=xf, aes(label=label4), x=-Inf, y=Inf,hjust=-0.25,vjust=2.75, parse=TRUE) +
          labs(title = title)
  }
  histo <- histo +
          #annotate("label", data=xf, inherit.aes=FALSE, aes(label=selection), x=0, y=Inf,hjust=0.5,vjust=-0.12, parse=FALSE) +
#          geom_label(data=xf, inherit.aes=FALSE, aes(label=selection), x=-Inf, y=Inf, hjust=0, vjust=-0.12, parse=FALSE) +
          #annotation_custom(g, xmin=-Inf, xmax=0.15,  ymin=-Inf, ymax=Inf+20) +
          theme_few() +
          #theme(panel.border=element_blank(), axis.line=element_line()) +
          #scale_linetype_manual(values=c("excellent"='solid', "good"='dotted', 'bad'='blank')) +
          scale_fill_skp() + scale_color_skp() + mytheme + guides(fill = guide_legend(nrow = 1), color=FALSE) +
          labs(
            fill = '',
            linetype = '',
            x = 'P-value p',
            y = 'Density of the histogram-bin corresponding to P-value p'
          ) + 
          theme(plot.title = element_text(hjust = 0.5, vjust= 0.25, margin=margin(0,0,-15,0)),
                strip.text.x = element_text(face="bold", margin=margin(1,1,1,1)))
  if (only_experiment) {
    histo <- histo + guides(fill=FALSE)
  }
  pg <- ggplot_build(histo)
  #highlight <- pg$data[[2]] %>%
  #               filter(xmax<=0.025, support=='yes') %>%
  #               mutate(ymin=0.0, colour='yellow', size=0.5)
  #pg$data[[2]] <- rbind(pg$data[[2]], highlight)
  gt <- ggplot_gtable(pg)
  gt$layout$clip[substr(gt$layout$name, 1, 5) == "panel"] <- "off"
#  gt <- gtable_add_rows(gt, unit(gt$heights[[3]], 'cm'), 2)
#  gt <- gtable_add_grob(gt, 
#                     list(rectGrob(gp = gpar(col = NA, fill = gray(0.5))),
#                          textGrob("Variable 2", gp = gpar(col = gray(1)))),
#                     3, 4, 3, 6, name = paste(runif(2)))
# gt <- gtable_add_rows(gt, unit(1/8, "line"), 3)
  
  #print(gt$layout)
  return(gt)
}


args = commandArgs(trailingOnly=TRUE)

infile = args[1]
shape_levels_file = args[2]
selected_pdf = args[3]
combined_pdf = args[4]
separate_pdf = args[5]


shape_levels <- read.csv(shape_levels_file, stringsAsFactors=FALSE)

#shape_levels$levels <- gsub(':', sprintf('\u2264'), shape_levels$levels)
shape_levels$new_levels_type <- shape_levels$levels_type
shape_levels$new_levels_type[shape_levels$levels_type == 'other1'] <- 'alternative'
shape_levels$new_levels_type[shape_levels$levels_type == 'other2'] <- 'alternative2'
shape_levels$new_levels_type[shape_levels$levels_type == 'publish'] <- 'main-text'
shape_levels$levels <- paste(paste(paste('shape', shape_levels$shape, sep=':'), paste('levels', shape_levels$new_levels_type, sep=':'), sep=', '), shape_levels$levels, sep='\n')
shape_levels <- shape_levels[, c('shape', 'levels_type', 'levels')]
shape_levels$shape <- factor(shape_levels$shape, levels=c('MGW', 'HelT', 'ProT', 'Roll'), ordered=TRUE)
shape_levels$levels_type <- factor(shape_levels$levels_type, levels=c('publish', 'other1', 'other2'), ordered=TRUE)
shape_levels <- shape_levels[order(shape_levels$shape, shape_levels$levels_type), ]
shape_levels$levels <- factor(shape_levels$levels, levels=shape_levels$levels, ordered=TRUE)

print(shape_levels)

all <- read.csv(infile, stringsAsFactors=FALSE, colClasses=c(en_th='character'))
all$shape <- factor(all$shape, levels=c('MGW', 'HelT', 'ProT', 'Roll'), ordered=TRUE)
all$ctx <- factor(all$ctx, levels=all_categories, ordered=TRUE)
all$family <- factor(all$family, levels=c('bHLH', 'ETS', 'homeodomain'), ordered=TRUE)


all$tag <- paste(all$family, all$shape, sep=', ')

all <- merge(all, shape_levels, all.x = TRUE)


print(head(all))

ets <- plotQvalueHisto(all[all$family=='ETS',], 'Family: ETS', levels~en_th)
homeodomain <- plotQvalueHisto(all[all$family=='homeodomain',], 'Family: homeodomain', levels~en_th)
bhlh <- plotQvalueHisto(all[all$family=='bHLH',], 'Family: bHLH', levels~en_th)
cairo_pdf(file=separate_pdf, bg='white', width=15, height=28, onefile=TRUE)
grid.newpage()
grid.draw(bhlh)
grid.newpage()
grid.draw(ets)
grid.newpage()
grid.draw(homeodomain)
dev.off()


combined <- plotQvalueHisto(all[all$en_th=='1.20',], 'Enrichment Threshold: 1.20', levels~family)
cairo_pdf(file=combined_pdf, bg='white', width=15, height=28, onefile=TRUE)
grid.newpage()
grid.draw(combined)
dev.off()

selected <- all[all$en_th=='1.20' & all$levels_type=='publish', ]
p <- plotQvalueHisto(selected, '', shape~family, TRUE)
cairo_pdf(file=selected_pdf, bg='transparent', width=7.34, height=5)
grid.draw(p)
dev.off()
