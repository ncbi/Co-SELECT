#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

infile = '../data/homeodomain/PITX3/TGCATC20NGA/4.non.TAAT.d2.bg.MGW.S4.080M5.130H5.750X.mer.pwm'
outfile = 'hoho.csv'

print(infile)
print(outfile)

df = read.csv(infile)

print(df)

getPWM <- function(df) {
    mat = as.matrix(df[, paste0('X', 1:10)])
    rownames(mat) = df$base
    return(mat)
}

all_cols = names(df)
mat_cols = c('base', all_cols[grepl("X\\d+", all_cols)])
group_cols = setdiff(all_cols, mat_cols)

df <- df %>%
          group_by_(.dots = group_cols) %>%
          nest()

df$data <- lapply(df$data, getPWM)


getMedianIC <- function(mat) {
  print(mat)
  return(1.0)
}

print(lapply(df$data, getMedianIC))

print(df)
quit()

df = df.groupby(['shapemer', 'seqmer']).agg({'count': sum}).reset_index('seqmer')

def makeMatrix(x):
  cnt = x['count']
  seq = x['seqmer']
  m = np.zeros((4, len(seq)), dtype=np.int64)
  for j, b in enumerate(seq):
    m['ACGT'.find(b), j] += cnt
  return(m)

df['mat'] = df.apply(makeMatrix, axis=1)

df['mat'].apply(lambda x: pd.DataFrame(x, index=pd.Index(list('ACGT'), name='base'), \
                                       columns=['X{}'.format(i+1) for i in range(x.shape[1])]).unstack()) \
         .groupby('shapemer').agg(sum) \
         .stack() \
         .reset_index() \
         .to_csv(outfile, index=False)

