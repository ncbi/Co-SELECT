#!/usr/bin/env Rscript

#####################
# Original code 
#####################

#loading libraries
library(ggplot2)
library(ggthemes)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)

pca <- read.csv(args[1], stringsAsFactors=F)
pca$label <- ifelse(pca$tf == 'ISX', paste(pca$tf, pca$barcode, sep='.'), pca$tf)

delta <- 0.1
ymax <- max(pca$pca.two) + delta
ymin <- min(pca$pca.two) - delta


# Create dummy variables 
bhlh <- unique(pca[pca$family=='bHLH', c('subcategory', 'family')])
ets <- unique(pca[pca$family=='ETS', c('subcategory', 'family')])
homeodomain <- unique(pca[pca$family=='homeodomain', c('subcategory', 'family')])


mycolors = read.csv(text="serial,subcategory,color,family
8,TAAA,#00BFC4,homeodomain
7,TAAT,#00B0F6,homeodomain
9,Other,#9590FF,homeodomain
4,CGGAA,#A3A500,ETS
5,GGAAG,#39B600,ETS
;6,Other,#00BF7D,ETS
1,CACGTG,#F8766D,bHLH
2,CATATG,#D89000,bHLH
3,Other,#E76BF3,bHLH", comment.char=';', stringsAsFactors=F)

mycolors = mycolors[order(mycolors$serial),]

homeodomain$subcategory[homeodomain$subcategory=='TAAA'] <- 'TAAA (mostly Hox)'
pca$subcategory[pca$subcategory=='TAAA'] <- 'TAAA (mostly Hox)'
mycolors$subcategory[mycolors$subcategory=='TAAA'] <- 'TAAA (mostly Hox)'

bhlh$subcategory = factor(bhlh$subcategory, levels=mycolors$subcategory[mycolors$family=='bHLH'], ordered=T)
homeodomain$subcategory = factor(homeodomain$subcategory, levels=mycolors$subcategory[mycolors$family=='homeodomain'], ordered=T)
ets$subcategory = factor(ets$subcategory, levels=mycolors$subcategory[mycolors$family=='ETS'], ordered=T)


print(bhlh)
print(ets)
print(homeodomain)

pdf(args[2], width=9, height=6)

pca$category = paste(pca$subcategory, pca$family, sep='.')
mycolors$category = paste(mycolors$subcategory, mycolors$family, sep='.')
print(pca)

## Create plot
ggplot() + 
## Plot the four dummy layer outside of the intended plotting area  
  geom_point(data = ets, aes(x = 0, y = ymax+1, shape = subcategory)) + 
  geom_point(data = homeodomain, aes(x = 0, y = ymax+1,fill = subcategory)) +
  geom_line(data = bhlh, aes(x = 0, y = ymax+1, linetype = subcategory)) +
## Add in your real plot goal  
  #geom_bar(data = services, aes(category, fill=subcategory)) +
  geom_point(data = pca, aes(x=pca.one, y=pca.two, color=category), size = 5, alpha = 0.8) +
  geom_text_repel(data = pca, aes(x=pca.one, y=pca.two, label=label), size = 1.5) +
## Remove the Fill legend  
  #scale_fill_hue(guide="none") +
  scale_alpha(guide="none") +
  scale_size(guide="none") +
  scale_color_manual(guide="none", values=setNames(mycolors$color, mycolors$category)) +
## Override the guide aesthetics to make them look like fill colors  
  guides(shape = guide_legend(order=2, override.aes = list(colour = mycolors$color[mycolors$family=='ETS'], fill = NA, shape = 15, size = 8),title="ETS", label.theme=element_text(family='mono')),
         linetype = guide_legend(order=1, override.aes = list(colour = mycolors$color[mycolors$family=='bHLH'], shape = 15, size = 8),title = "bHLH", label.theme=element_text(family='mono')),
         fill = guide_legend(order=3, override.aes = list(colour = mycolors$color[mycolors$family=='homeodomain'], shape = 15,size = 8), title = "homeodomain", label.theme=element_text(family='mono'))) +
## Adjust the plot range to hide all the extra layers  
  ylim(ymin,ymax) +
  theme_few() +
  labs(x="PCA dimension 1",
       y="PCA dimension 2",
       title="Clustering of TF Experiments by PCA of Shapemers Enriched in Motif-free Aptamers") 
#> Warning: Removed 3 rows containing missing values (geom_point).
#> Warning: Removed 2 rows containing missing values (geom_point).
#> Warning: Removed 3 rows containing missing values (geom_point).
#> geom_path: Each group consist of only one observation. Do you need to adjust the group aesthetic?
