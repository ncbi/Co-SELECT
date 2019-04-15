# Load the required packages
suppressMessages(require(ggplot2))
suppressMessages(require(spseqlogo))
suppressMessages(require(plyr))
suppressMessages(require(dplyr))
suppressMessages(require(tidyr))
suppressMessages(require(rlist))
suppressMessages(require(grid))
suppressMessages(require(gridExtra))

cs = make_col_scheme(
       chars = c('A', 'C', 'G', 'T', 'U', 'X', 'B', 'D', 'E', 'F', 'H', 'I', 'J'),
       groups = c('A', 'C', 'G', 'T', 'U', 'full', 'fgfull', 'motif-containing', 'motif-free', 'fg', 'bg', 'motif-containing-total', 'motif-containing-enriched'),
       cols = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839', 'white', 'green', 'blue', 'red', 'blue', 'red', 'blue', 'blue'))

getPwm <- function(df) {
  mat <- as.matrix(df[, !(names(df) %in% c('base'))])
  if (all(mat == 0)) {
    mat <- cbind(mat, c(1,0,0,0))
  }
  row.names(mat) <- df$base
  mat
}

getLogoPlotOld <- function(pwms) {
  all_cols = names(pwms)
  mat_cols = c('base', all_cols[grepl("X\\d+", all_cols)])
  group_cols = setdiff(all_cols, mat_cols)
  seq_length = length(mat_cols) - 1
  
  pwms <- pwms %>%
            group_by_(.dots = group_cols) %>%
            nest()
  
  pwms$data <- lapply(pwms$data, getPwm)
  pwms$nseq <- lapply(pwms$data, function(x) max(colSums(x[, 1:seq_length])))
  pwms <- pwms %>% unnest(nseq)

  pwms$context <- ifelse(pwms$context == 'fg', 'motif-containing', 'motif-free')

  print(pwms)
  
  p <- ggplot() + 
         geom_logo(data = pwms, col_scheme = cs) +
         theme_logo() +
         theme_bw() +
         theme(panel.background = element_rect(fill = NA, color = "black")) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         theme(legend.position='bottom', legend.box = "horizontal") +
         guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
         geom_rect(data = pwms, aes(fill = context), xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1) +
         scale_x_continuous(breaks=1:seq_length, lim=c(0.5, seq_length + 0.5)) +
         geom_label(data = pwms, aes(label = nseq), x=-Inf, y=Inf, hjust = -0.1, vjust = +1.25) +
         ylim(0,2)
}


getLogoPlot <- function(pwms) {
  all_cols = names(pwms)
  mat_cols = c('base', all_cols[grepl("X\\d+", all_cols)])
  group_cols = setdiff(all_cols, mat_cols)
  seq_length = length(mat_cols) - 1
  
  pwms <- pwms %>%
            group_by_(.dots = group_cols) %>%
            nest()
  
  pwms$data <- lapply(pwms$data, getPwm)
  pwms$nseq <- lapply(pwms$data, function(x) max(colSums(x[, 1:seq_length])))
  pwms <- pwms %>% unnest(nseq)

  pwms$context <- ifelse(pwms$context == 'fg', 'motif-containing', 'motif-free')

  print(pwms)
  
  p <- ggplot() + 
         geom_logo(data = pwms, col_scheme = cs) +
         theme_logo() +
         theme_bw() +
         theme(panel.background = element_rect(fill = NA, color = "black")) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         theme(legend.position='bottom', legend.box = "horizontal") +
         guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
         #geom_rect(data = pwms, aes(fill = context), xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.1) +
         scale_x_continuous(breaks=1:seq_length, lim=c(0.5, seq_length + 0.5)) +
         #geom_label(data = pwms, aes(label = nseq), x=-Inf, y=Inf, hjust = -0.1, vjust = +1.25) +
         theme(axis.line = element_blank()) +
         ylim(0,2)
}

