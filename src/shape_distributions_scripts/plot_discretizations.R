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
adjusts = c(
 'MGW' = 1,
 'HelT' = 2,
 'ProT' = 1,
 'Roll' = 1
)

top = data.frame()
df = data.frame()
limits = data.frame()
gaussians = tibble()
functions = tibble()
for (shape in shapes) {
  load(paste(top_dir, paste0('unmixed_gaussians_', shape, '_k', comp[shape], '.dist'), sep = '/'))
  top = rbind(top, data.frame(x = mixmdl$x, shape = shape, plot = 'A'))
  gaussians = bind_rows(gaussians, tibble(mu = mixmdl$mu, sigma = mixmdl$sigma, lambda=mixmdl$lambda, comp=paste0('C',1:comp[shape]), shape = shape, plot='A'))
  den = density(mixmdl$x, adjust=adjusts[shape])
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


df = ddply(df, .(shape), mutate, relden = den/max(den, na.rm=T))


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

colors = c(
'C0' = '#1b9e77',
'C1' = '#d95f02',
'C2' = '#7570b3',
'C3' = '#e7298a'
)

p <- ggplot(top, aes(x)) +
       geom_histogram(aes(y = ..density..), color= 'grey', alpha = 0.4, linetype='solid', size = 0.25, fill=NA) +                        
       geom_line(data = normaldens, aes(y = density, color= comp), size=0.5, alpha = 0.9) +
       geom_segment(data = vlines, aes(x = cutoff, y=0, xend = cutoff, yend=yend, color=comp, linetype=selected), size = 0.5, alpha = 0.9) +
       geom_line(data = df %>% mutate(plot='A'), aes(y = den, color = 'C0'), size = 0.7, alpha = 0.9) +  
       #geom_line(aes(y = ..density..), color = 'black', size = 0.5, alpha = 0.9, linetype = 'dotted', stat = 'density') +  
       geom_raster(data = df, aes(y = y, fill = relden)) +
       geom_vline(data = cutoffs, aes(xintercept = cutoff), alpha = 0.9) +
       geom_text(data = cutoffs, aes(x = cutoff, label = cutoff), y=-Inf, vjust = +1.5, family='serif', size=3, color='grey30') +
       geom_text(data = labels, aes(x = pos, label = label), y=0.025) +
       facet_grid(plot~shape, scale = "free", space='free_y') +
       scale_x_continuous(sec.axis = dup_axis(name=NULL, labels=NULL)) +
       scale_linetype_manual(values = c('yes'='solid', 'no'='dashed'),
                             labels = c('yes'='SD and/or chosen as a cutoff', 'no' = 'Standard deviation (SD)')) +
       scale_fill_gradient(low = 'gray', high = 'red') +
       scale_color_manual(values = c('C0' = 'red', 'C1' = 'darkgreen', 'C2' = 'blue', 'C3' = 'brown'),
                          labels = c('C0' = 'Original',
                                     'C1' = 'One Gaussian component', 
                                     'C2' = 'Other Gaussian component')) +
       labs(x = "Value of a shape feature", y='Density') +
       #scale_color_manual(values = colors) +
       theme_tufte() + guides(fill = FALSE, 
                              linetype = guide_legend(title='Vertical line:', label.hjust=1, order=2),
                              color=guide_legend(title='Density curve:', label.hjust=1, order = 1)) +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
       theme(strip.text.y = element_blank()) +
       theme(panel.spacing.y = unit(-0.05, "lines")) +
       theme(panel.border = element_blank()) +
       theme(axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")) +
       # Remove only the legend title
       theme(legend.position='bottom') +
       theme(legend.text = element_text(margin = margin(r = 10, unit = "pt"))) +
       theme(plot.margin = margin(1, 1, 0, 1, "mm")) +
       theme(axis.title.x = element_text(vjust=-4.5)) 


gt = ggplot_gtable(ggplot_build(p))

panels <- grep("panel", gt$layout$name)
top <- unique(gt$layout$t[panels])
## min(top) gives the row of top row of panels
## min(top) - 2 is the extra top axis, min(top) - 1 is the strip
## max(top) gives the row of bottom row of panels
## max(top) +1 gives the row of bottom axes
## intersperse a copy of the bottom axes
new_gt <- gtable:::rbind_gtable(gt[seq.int(min(top)-3), ], gt[seq.int(min(top)-1, min(top)),], "first")
new_gt <- gtable:::rbind_gtable(new_gt, gt[max(top)+1, ], "first")
new_gt <- gtable:::rbind_gtable(new_gt, gt[min(top)-2, ], "first")
new_gt <- gtable:::rbind_gtable(new_gt, gt[seq(min(top)+1, max(top)),], "first")
new_gt <- gtable:::rbind_gtable(new_gt, gt[seq(max(top)+2, nrow(gt)),], "first")

# Get the ylabel grob
ylabl = gtable_filter(gt, "ylab-l", trim=F)

# Insert the strip grobs into the new rows
new_gt = gtable_add_grob(new_gt, ylabl$grobs[[1]],  t=min(top)-1, l=ylabl$layout$l[[1]], b=max(top), r=ylabl$layout$r[[1]], name='ylab-l')

gt <- new_gt
gt <- gtable_filter(gt, 
                     "(background|panel|strip|axis-[tbr]|xlab|ylab|guide-box|title|subtitle|caption|tag|axis-l-1)",
                     trim=FALSE)
gt$layout$clip[grepl("^panel", gt$layout$name)] <- "off"
grid.draw(gt)

