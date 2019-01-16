#!/usr/bin/env Rscript

library(ggplot2)
library(mixtools)
library(grid)
library(gtable)
library(ggthemes)
library(plyr)
library(tidyr)
library(dplyr)

top_dir = '/panfs/pan1/dnashape/ShapeDistributions'
shapes = c('MGW', 'HelT', 'ProT', 'Roll')
comp = c('MGW' = 2, 'HelT' = 2, 'ProT' = 2, 'Roll' = 2)

top = data.frame()
df = data.frame()
limits = data.frame()
gaussians = tibble()
functions = tibble()
for (shape in shapes) {
  load(paste(top_dir, paste0('unmixed_gaussians_', shape, '_k', comp[shape], '.dist'), sep = '/'))
  top = rbind(top, data.frame(x = mixmdl$x, shape = shape, plot = 'A'))
  gaussians = bind_rows(gaussians, tibble(mu = mixmdl$mu, sigma = mixmdl$sigma, lambda=mixmdl$lambda, comp=paste0('C',1:comp[shape]), shape = shape, plot='A'))
  den = density(mixmdl$x)
  df = rbind(df, data.frame(x = den$x, shape = shape, den = den$y, y=rep(c(0.00, 0.05), each=length(den$x)), plot = 'B'))
  limits = rbind(limits, data.frame(shape=shape, max=max(den$x), min=min(den$x)))
  functions = bind_rows(functions, tibble(shape=shape, f=list(approxfun(den$x, den$y))))
}


print(gaussians)
print(functions)

dis <- read.csv('shape_levels.csv', stringsAsFactors=FALSE)
dis <- dis[dis$levels_type == 'publish', c('shape', 'levels')]
dis = merge(dis, limits)
dis$levels = paste(dis$min, dis$levels, dis$max, sep='|')
print(dis)

getCutoffLabelPos <- function(x) {
  print(x)
  x = unlist(strsplit(x, '[|]'))
  vals = as.numeric(x[c(T,F)])
  labs = x[c(F,T)]
  n = length(vals)
  pos = (vals[2:n] + vals[1:(n-1)])/2
  tibble(cutoff = list(vals[2:(n-1)]), label = list(labs), pos = list(pos))
}

dis <- dis %>% 
    group_by(shape) %>%
    do(getCutoffLabelPos(.$levels))

print(dis)



applyFunc <- function(f, x) {
  lapply(seq_along(f), function(n) {
    do.call(f[[n]], list('v' = x[[n]]$cutoff))
  })
}

cutoffs = dis %>%
    select(shape, cutoff) %>%
    unnest() %>%
    mutate(selected = 'yes') %>%
    group_by(shape) %>%
    nest() %>%
    inner_join(functions, by='shape') %>%
    mutate(yend = applyFunc(f, data)) %>%
    select(-f) %>%
    unnest() %>%
    mutate(plot = 'B')

print(cutoffs)


labels <- dis %>% 
    select(shape, label, pos) %>%
    unnest() %>%
    mutate(plot = 'B')
print(labels)

print('computing vlines')

print(gaussians)


deleteKnown <- function(x) {
  lapply(x, function(x) {
    for (i in 1:nrow(x)) {
      for (j in 1:nrow(x)) {
        if ((x$plot[i] == 'B') && (x$plot[j] == 'A')) {
          if (abs(x$cutoff[i] - x$cutoff[j]) < 0.1) {
            x$comp[i] = 'CX'
            x$selected[j] = x$selected[i]
          }
        }
      }
    }
    x = x[x$comp != 'CX',]
  })
}

vlines = gaussians %>%
    group_by(shape, comp) %>%
    mutate(cutoff = list(c(mu+sigma, mu-sigma))) %>%
    mutate(yend = list(lambda*dnorm(c(mu+sigma, mu-sigma), mean=mu, sd=sigma))) %>%
    unnest(cutoff, yend) %>%
    select(shape, comp, plot, cutoff, yend) %>%
    mutate(selected = 'no') %>%
    bind_rows(cutoffs %>% mutate(comp='C0')) %>%
    group_by(shape) %>%
    nest() %>%
    mutate(data = deleteKnown(data)) %>%
    unnest() %>%
    mutate(plot = 'A')


df = ddply(df, .(shape), mutate, den = den/max(den, na.rm=T))


pdf("Rplots.pdf", width=12, height=4)

df$shape = factor(df$shape, levels=shapes)
top$shape = factor(top$shape, levels=shapes)
labels$shape = factor(labels$shape, levels=shapes)
cutoffs$shape = factor(cutoffs$shape, levels=shapes)
vlines$shape = factor(vlines$shape, levels=shapes)

normaldens <- ddply(gaussians, c("shape", "plot", "comp"), function(df) {
  shape = unique(df$shape)
  plot = unique(df$plot)
  comp = unique(df$comp)
  lim = limits[limits$shape == shape,]
  x = with(lim, seq(min, max, length = 100))
  data.frame(
    x = x, plot = plot, comp = comp,
    density = df$lambda * dnorm(x, df$mu, df$sigma)
  )
})

p <- ggplot(top, aes(x)) +
       geom_histogram(aes(y = ..density..), color= 'grey', alpha = 0.4, linetype='solid', size = 0.25, fill=NA) +                        
       geom_line(data = normaldens, aes(y = density, color= comp), size=0.5, alpha = 0.9) +
       geom_segment(data = vlines, aes(x = cutoff, y=0, xend = cutoff, yend=yend, color=comp, linetype=selected), size = 0.5, alpha = 0.9) +
       geom_line(aes(y = ..density..), color = 'red', size = 1, alpha = 0.9, stat = 'density') +  
       geom_raster(data = df, aes(y = y, fill = den)) +
       geom_vline(data = cutoffs, aes(xintercept = cutoff), alpha = 0.9) +
       geom_text(data = cutoffs, aes(x = cutoff, label = cutoff), y=+Inf, vjust = -0.25, size=3) +
       geom_text(data = labels, aes(x = pos, label = label), y=0.025) +
       facet_grid(plot~shape, scale = "free", space='free_y') +
       scale_linetype_manual(values = c('yes'='solid', 'no'='dashed')) +
       scale_fill_gradient(low = 'gray', high = 'red') +
       scale_color_manual(values = c('C0' = 'brown', 'C1' = 'darkgreen', 'C2' = 'blue', 'C3' = 'brown')) +
       theme_bw() + guides(fill = FALSE, linetype = FALSE, color = FALSE) +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
       theme(strip.text.y = element_blank()) +
       theme(panel.spacing.y = unit(0.8, "lines"))


gt = ggplot_gtable(ggplot_build(p))
gt <- gtable_filter(gt, 
                     "(background|panel|strip|axis-[tbr]|xlab|ylab|guide-box|title|axis-l-1)",
                     trim=FALSE)
gt$layout$clip[grepl("^panel", gt$layout$name)] <- "off"
grid.draw(gt)

