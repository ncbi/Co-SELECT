# Load the required packages
require(ggplot2)
require(ggseqlogo)
require(plyr)
require(dplyr)
require(tidyr)
require(rlist)
require(gridExtra)

pwms = read.csv('SEQLOGO/seqlogo_info_PITX3_TGCATC20NGA_NEW.csv', stringsAsFactors=F)
pwms$cycle <- 4
print(pwms)

cycles <- unique(pwms$cycle)
contexes <- unique(pwms$context)

cs = make_col_scheme(
       chars = c('A', 'C', 'G', 'T', 'U', 'X', 'B', 'D', 'E'),
       groups = c('A', 'C', 'G', 'T', 'U', 'full', 'fgfull', 'motif-containing', 'motif-free'),
       cols = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839', 'white', 'green', 'blue', 'red'))


getPwm <- function(df) {
  mat <- as.matrix(df[, !(names(df) %in% c('base'))])
  if (all(mat == 0)) {
    mat <- cbind(mat, c(1,0,0,0))
  }
  row.names(mat) <- df$base
  mat
}

pwms <- pwms %>%
          group_by(shapemer, cycle, context) %>%
          nest() %>%
          mutate(name = paste(shapemer, cycle, context))

mylist <- lapply(pwms$data, getPwm)
names(mylist) <- pwms$name
pwms$context <- ifelse(pwms$context == 'fg', 'motif-containing', 'motif-free')

pwms$seq_group = factor(pwms$name, ordered=TRUE)


d <- ggseqlogo(mylist, col_scheme=cs, ncol=length(contexes)*length(cycles)) + theme(panel.background = element_rect(fill = NA, color = "black")) +
  ylim(0,2) +
  geom_rect(data = pwms, aes(fill = context), xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1) +
  #ggtitle(paste(shape, 'shapemers for', tf, primer))+
  scale_x_continuous(breaks=1:10, lim=c(0.5,10.5))

print(d)

