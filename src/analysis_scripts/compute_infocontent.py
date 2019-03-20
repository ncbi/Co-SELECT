#!/usr/bin/env python

from itertools import groupby
import subprocess as sp
import pandas as pd
import numpy as np
import sys

infile = '../data/homeodomain/PITX3/TGCATC20NGA/4.non.TAAT.d2.bg.MGW.S4.080M5.130H5.750X.mer.pwm'
outfile = 'haha.csv'

if len(sys.argv) > 1:
  infile = sys.argv[1]
  outfile = sys.argv[2]

print infile
print outfile

df = pd.read_csv(infile)
print(df)
df = df.set_index('shapemer').unstack()
print(df)
exit(0)

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

