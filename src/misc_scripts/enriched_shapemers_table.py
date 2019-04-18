import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__))+'/..')
from results_scripts.more_results import *
import pandas as pd

top_results_dir = '../results/d0'
ltype = 'publish'
cycle = 4
lflank = rflank = 1

infile = '%s/%s/denriched_same.%d.l%d.r%d.csv' % (top_results_dir, ltype, cycle, lflank, rflank)
tf_motif_file = './tf_coremotif.csv'
outfile = './table_enriched_shapemers.tex'

print(infile)

df = pd.read_csv(infile)
df = df[(df['en_th'] == 1.2) & (df['label'] == 'both')]

print(df.head())



df = get_enriched_shapemers_in_dataframe(df, tf_motif_file)
df = df.dropna(axis=1, how='all')

df = df.fillna('')
print(df)
print(df.columns)
print(len(df.columns))

df.columns = df.columns.droplevel()


cols = pd.DataFrame(index = df.columns).reset_index('shapemer').groupby('shape').agg({'shapemer':len}).reset_index()
csum = [0] + cols['shapemer'].cumsum().tolist()
counts = cols['shapemer'].tolist()
print cols

import re

lines = df.to_latex(multicolumn=True, multirow=True, multicolumn_format='c', bold_rows=True, column_format='ll>{\\ttfamily}c'+'c'*sum(counts))
lines = lines.split('\n')
print(lines[3])
lines[3] = re.sub('([SMHX]+)', r'\\rotatebox{90}{\\texttt{\1}}', lines[3])
print(lines[3])
lines[3] = re.sub(r'(\\textbf{shapemer})', r'\\rotatebox{30}{\\textbf{\1}}', lines[3])
print(lines[3])
lines.insert(3, ' '.join(['\\cmidrule(%s){%d-%d}' % ('r'*counts[i], 4+csum[i], 4+csum[i+1]-1) for i in range(len(counts))]))
#del lines[4]
for i in range(len(lines)):
  print "**********"
  print lines[i]
  lines[i] = re.sub(r"\\cline{2-76}", '', lines[i])
  lines[i] = re.sub(r'(\\textbf{homeodomain|ETS|bHLH})', r'\\rotatebox{90}{\\textbf{\1}}', lines[i])
  lines[i] = re.sub(r'(\\textbf{family})', r'\\rotatebox{45}{\\textbf{\1}}', lines[i])
  #lines[i] = re.sub(r'(\bX\b)', r'\\ding{51}', lines[i]) # tick
  lines[i] = re.sub(r'(\bX\b)', r'\\ding{54}', lines[i])
  print lines[i]

with open(outfile, 'w') as f:
    f.write('\n'.join(lines))

import pylatex as pl
from pylatex.utils import italic, NoEscape
from pylatex.package import Package
from pylatex.base_classes import Command
import os

doc = pl.Document('basic')
doc.documentclass = Command(
        'documentclass',
        options=['12pt', 'landscape'],
        arguments=['standalone'],
    )
doc.packages.append(Package('array'))
doc.packages.append(Package('graphicx'))
doc.packages.append(Package('booktabs'))
doc.packages.append(Package('multirow'))
doc.packages.append(Package('pifont'))
doc.append(NoEscape('\n'.join(lines)))

doc.generate_pdf(os.path.splitext(outfile)[0]+'_standalone', clean_tex=True)
