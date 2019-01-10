#!/usr/bin/env Rscript
require(plyr)
require(ggplot2)
require(ggthemes)
require(gridExtra)
require(grid)
require(dplyr)

inset <- read.csv(text='en_th, shape, height, xmin, ymin
1.10,MGW,  0.50, 20, 0.15
1.10,HelT, 0.70, 40, 0.32
1.10,ProT, 0.65, 35, 0.35
1.10,Roll, 0.15, 1.6, 0.15
1.20,MGW,  0.33, 15, 0.20
1.20,HelT, 0.33, 20, 0.20
1.20,ProT, 0.33, 30, 0.25
1.20,Roll, 0.33, 35, 0.20
', stringsAsFactors = FALSE, colClasses=c("en_th"="character"))

shapes <- c("MGW", "HelT", "ProT", "Roll")
en_ths <- c('1.10', '1.20')


color_values = c(
'Low @ 1.20'  = '#808080',
'Low'         = '#808080',
'High @ 1.20' = '#f58231',
#'High'        = '#e6beff',
#'High'        = '#aaffc3',
'High'        = '#f58231',
'High @ 1.10' = '#808000',
'Low @ 1.10'  = '#aa6e28'
)
labels = c(
'High @ 1.20' = 'High @ 1.20',
'High'        = 'High'       ,
'Low @ 1.20'  = 'Low @ 1.20' ,
'Low'         = 'Low'        ,
'High @ 1.10' = 'High @ 1.10',
'Low @ 1.10'  = 'Low @ 1.10' 
)


annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data=NULL) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
      geom = ggplot2:::GeomCustomAnn,
      inherit.aes = TRUE, params = list(grob = grob, 
        xmin = xmin, xmax = xmax, 
        ymin = ymin, ymax = ymax))
}

mytheme <- theme(legend.position="top",
                 legend.margin=margin(1,1,1,1),
                 legend.justification="right",
                 legend.box.margin=margin(1,1,-07,1),
                 legend.background = element_rect(colour = 'gray', fill = '#FEFEF1'),
                 #plot.background = element_rect(fill = "transparent",colour = NA),
                 #plot.margin = margin(2,1,1,1),
                 #panel.spacing.y=unit(2, "lines"),
                 #strip.background=element_blank()
                 )


plotEnrichedShapes <- function(df) {
  #print(df)
  #print('in plotenrich')
  #print(unique(df$shape))
  p <- ggplot(df, aes(x=x, y=promiscuity, fill=promiscuous)) +
        #facet_wrap(~shape, scales='free',ncol=2) +
        facet_wrap(~shape, ncol=2) +
        geom_bar(stat='identity', position='dodge') +
        theme_few() +
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        #mytheme +
        #scale_fill_skp() +
        #theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.35))+
        #theme(legend.position="none") +
        #scale_x_discrete(expand = c(0, 2)) + 
        scale_fill_manual(values=color_values, labels=labels) +
        labs(x="Shapemer x", y="Mimimum fraction of TFs in any family for\nwhich x is enriched in aptamers with no-motif", fill='') +
        guides(fill=FALSE)
}


getPlotWithInset <- function(all) {

  tmp <- ddply(all, .(shape, kmer), summarize, maxpromiscuity = max(promiscuity))
  tmp <- ddply(tmp, .(shape), mutate, x = rank(-maxpromiscuity, ties.method="first"))
  all <- merge(all, tmp, all.x=TRUE, id.vars=c('shape'))

  all <- ddply(all, .(shape), transform, xmin=max(xmin), ymin=max(ymin))
  #print(length(unique(all$en_th)))
  #print(unique(all$en_th))
  if (length(unique(all$en_th)) > 1 ) {
  all$promiscuous <- paste(all$promiscuous, all$en_th, sep=' @ ')
  all$promiscuous <- factor(all$promiscuous, levels=c('High @ 1.10', 'High @ 1.20', 'Low @ 1.10', 'Low @ 1.20'), ordered=TRUE)
  } else {
  all$promiscuous <- factor(all$promiscuous, levels=c('High', 'Low'), ordered=TRUE)
  }
  #print(all)

  limits <- merge(all, inset, all.x=TRUE) %>% 
            group_by(shape) %>%
            filter(promiscuous %in% c('High @ 1.10', 'High @ 1.20', 'High')) %>%
            summarize(xmax=max(x)+1)
  
  #print('before calling plotenrich')
  #print(unique(all$shape))
  p <-  plotEnrichedShapes(all)

  with_limits <- merge(all,limits)
  with_limits <- with_limits[with_limits$x < with_limits$xmax, ]
  
  insets <- dlply(with_limits, c("shape"), function(d){
      annotation_custom2(grob = ggplotGrob(plotEnrichedShapes(d) +
                                   geom_text(aes(label=kmer), angle=90, y=-Inf, hjust=-0.2,family='Courier') +
                                   scale_y_continuous(position = "right")+
          theme(legend.position="none") +
                                   theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                                        strip.text.x = element_blank())), 
        data = d,
        xmin=max(d$xmin), xmax=Inf, ymin=max(d$ymin), ymax=Inf)
      })
  
  p <- p + insets
}

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outfile = args[2]
selected_pdf = args[3]
combined_pdf = args[4]

print(outfile)
print(selected_pdf)
print(combined_pdf)

print(infile)

all <- read.csv(infile, stringsAsFactors=F, colClasses=c("en_th"="character"))

#print(all)
#print(all[all$en_th == '1.20', ])


all <- merge(all, inset, all.x = TRUE)
all <- ddply(all, .(shape, kmer, en_th), transform, promiscuous = ifelse(promiscuity > height, 'High', 'Low'))
all$shape <- factor(all$shape, levels=c('MGW', 'HelT', 'ProT', 'Roll'))

#print(head(all))

write.csv(all[all$promiscuous == 'High', c('shape', 'kmer', 'en_th')], file=outfile, row.names=FALSE)


p <- getPlotWithInset(all)
pdf(file=combined_pdf, width=15, height=9)
print(p)

p <- getPlotWithInset(all[all$en_th == '1.20',]) 
pdf(file=selected_pdf, width=10, height=6)
print(p)



