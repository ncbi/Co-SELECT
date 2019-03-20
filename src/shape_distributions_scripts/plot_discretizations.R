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

dis <- read.csv('shape_levels.csv')
dis <- dis[dis$levels_type == 'publish', c('shape', 'levels')]
dis = merge(dis, limits)
dis$levels = paste(dis$min, dis$levels, dis$max, sep='|')
print(dis)

getCutoffLabelPos <- function(x) {
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

# now dis is nested, cannot be unnested because the
# numbers of cutoffs and labels are different
print(dis)


applyFunc <- function(f, x) {
  lapply(seq_along(f), function(n) {
    do.call(f[[n]], list('v' = x[[n]]$cutoff))
  })
}

# unnest cutoffs from dis and also compute the y values
# of the original density plot at the cutoff levels 
cutoffs = dis %>%
    select(shape, cutoff) %>%
    unnest() %>%
    mutate(selected = 'yes') %>%
    group_by(shape) %>%
    nest() %>%
    inner_join(functions, by='shape') %>%
    mutate(yend = applyFunc(f, data)) %>%
    select(-f) %>%
    unnest()

print(cutoffs)


labels <- dis %>% 
    select(shape, label, pos) %>%
    unnest() %>%
    mutate(plot = 'B')

print(labels)

setCutoffSource <- function(x) {
  lapply(x, function(x) {
    for (i in 1:nrow(x)) {
      for (j in 1:nrow(x)) {
        if ((x$plot[i] == 'B') && (x$plot[j] == 'A')) {
          if (abs(x$cutoff[i] - x$cutoff[j]) < 0.1) {
            x$comp[i] = x$comp[j]
            x$plot[i] = 'AB'
          }
        }
      }
    }
    x
  })
}

# For a vertical line, it could be either a standard deviate of
# a component or it is an adhoc cutoff based on intersection / local minima
# In the first case it has comp = C1, C2 etc and in the second case 
# it has comp = CX

# Initially we start with all cutoffs but later eliminate those which 
# have a counterpart in the standard deviates

vlines = gaussians %>%
    group_by(shape, comp) %>%
    mutate(cutoff = list(c(mu+sigma, mu-sigma))) %>%
    mutate(yend = list(lambda*dnorm(c(mu+sigma, mu-sigma), mean=mu, sd=sigma))) %>%
    unnest(cutoff, yend) %>%
    select(shape, comp, plot, cutoff, yend) %>%
    mutate(selected = 'no') %>%
    bind_rows(cutoffs %>% mutate(comp='CX', plot='B')) %>%
    group_by(shape) %>%
    nest() %>%
    mutate(data = setCutoffSource(data)) %>%
    unnest()

vlines_A = vlines %>% filter(plot != 'AB') %>% mutate(plot = 'A')
vlines_B = vlines %>% filter(plot != 'A') %>% mutate(plot = 'B')

df = ddply(df, .(shape), mutate, relden = den/max(den, na.rm=T))


pdf("fig_discretization.pdf", width=10, height=4)

df$shape = factor(df$shape, levels=shapes)
top$shape = factor(top$shape, levels=shapes)
labels$shape = factor(labels$shape, levels=shapes)
cutoffs$shape = factor(cutoffs$shape, levels=shapes)
vlines_A$shape = factor(vlines_A$shape, levels=shapes)
vlines_B$shape = factor(vlines_B$shape, levels=shapes)
gaussians$shape = factor(gaussians$shape, levels=shapes)

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

dots = data.frame(x=c(4.65, 32.04, -9.32, -4.99), y=c(0.22, 0.03, 0.088, 0.058), shape=c('MGW', 'HelT', 'ProT', 'ProT'), plot = 'A')


p <- ggplot(top, aes(x)) +
       geom_histogram(aes(y = ..density..), color= 'grey', alpha = 0.4, linetype='solid', size = 0.25, fill=NA) +                        
       geom_line(data = normaldens, aes(y = density, color= comp), size=0.5, alpha = 0.9) +
       geom_segment(data = vlines_A, aes(x = cutoff, y=0, xend = cutoff, yend=yend, color=comp, linetype=selected), size = 0.5, alpha = 0.9) +
       geom_line(data = df %>% mutate(plot='A'), aes(y = den, color = 'C0'), size = 0.7, alpha = 0.9) +  
       geom_raster(data = df, aes(y = y, fill = relden)) +
       geom_segment(data = vlines_B, aes(x = cutoff, y = -Inf, xend = cutoff, yend = +Inf, color = comp), alpha = 0.9) +
       geom_text(data = vlines_B, aes(x = cutoff, label = cutoff), y=-Inf, vjust = +1.5, family='serif', size=3, color='grey30') +
       geom_text(data = labels, aes(x = pos, label = label), y=0.025) +
       geom_point(data = dots, aes(x=x, y=y), size=3, color='darkviolet') +
       facet_grid(plot~shape, scale = "free", space='free_y') +
       scale_x_continuous(sec.axis = dup_axis(name=NULL, labels=NULL)) +
       scale_linetype_manual(values = c('yes'='solid', 'no'='31'),
                             labels = c('yes'='Cutoff', 'no' = 'Standard deviation')) +
       scale_fill_gradient(low = 'gray', high = 'red') +
       scale_color_manual(values = c('C0' = 'red', 'CX' = 'darkviolet', 'C1' = 'darkgreen', 'C2' = 'blue', 'C3' = 'brown'),
                          labels = c('C0' = 'Original density', 'CX' = 'Intersection / local minima',
                                     'C1' = 'One Gaussian component', 
                                     'C2' = 'Other Gaussian component')) +
       labs(x = "Value of a shape feature", y='Density') +
       theme_tufte() + guides(fill = FALSE, 
                              linetype = guide_legend(title='', label.hjust=1, order=2),
                              color=guide_legend(title='', label.hjust=1, order = 1)) +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
       theme(strip.text.y = element_blank()) +
       theme(panel.spacing.y = unit(-0.05, "lines")) +
       theme(panel.border = element_blank()) +
       theme(axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")) +
       # Remove only the legend title
       theme(legend.position='bottom') +
       theme(legend.text = element_text(margin = margin(r = 5, unit = "pt"))) +
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

