#!/usr/bin/env python
import os

lflank = rflank = 1
top_results_dir = '/panfs/pan1/aptax/FinalResults'
infile = "{}/dqvalue_l{}.r{}.csv".format(top_results_dir, lflank, rflank)
out_tex = 'out.tex'
out_pdf = 'out.pdf'

print infile

import argparse
parser = argparse.ArgumentParser(description='Generate main results of Co-SELECT in a table.')
parser.add_argument("qvalue_file", default=infile, help="input file containing all qvalues computed by Co-SELECT")
parser.add_argument("out_latex_file", default=out_tex, help="output latex file containing the tabular section only")
parser.add_argument("out_pdf_file", default=out_pdf, help="output standalone pdf file containing the table only")
opt = parser.parse_args()

print opt

import pandas as pd

df = pd.read_csv(opt.qvalue_file, dtype={'en_th':object})

df = df[df['ctx'] == 'experiment']
df = df[df['levels_type'] == 'publish']
df = df[df['en_th'].isin(['1.10', '1.20'])]

df = df[['tf.x', 'primer.x', 'family', 'shape', 'rvalue', 'en_th']]
df = df.rename(index=str, columns={"family": "Family", "shape": "Shape", 'en_th':'Threshold'})

percent_significant_one = lambda x: sum(x <= 0.05)*100.0/len(x)
percent_significant_one.__name__ = '0.05'

percent_significant_two = lambda x: sum(x <= 0.10)*100.0/len(x)
percent_significant_two.__name__ = '0.10'


df = df.groupby(['Family', 'Shape', 'Threshold']).agg({'rvalue': [percent_significant_one, percent_significant_two]}).round(2)
df.columns = df.columns.droplevel().set_names('FDR')
df = df.stack()

df = df.reset_index()

shapes = ['MGW', 'HelT', 'ProT', 'Roll']

df['Family'] = pd.Categorical(df['Family'], categories=['bHLH', 'ETS', 'homeodomain'], ordered=True)
df['Shape'] = pd.Categorical(df['Shape'], categories=shapes, ordered=True)

df = df.set_index(['Family', 'Threshold', 'FDR', 'Shape'])

df = df.unstack().unstack()
df.columns = df.columns.droplevel()
#df.index.set_names('', inplace=True)

lines = df.to_latex(multicolumn=True, multirow=True, column_format='lc'+('' + 'r'*2)*len(shapes), multicolumn_format='c', bold_rows=True)
lines = lines.split('\n')
lines.insert(3, ' '.join(['\\cmidrule(lr){%d-%d}' % (2*i+df.index.nlevels+1, 2*i+df.index.nlevels+2) for i in range(len(shapes))]))
#del lines[4]

with open(opt.out_latex_file, 'w') as f:
    f.write('\n'.join(lines))

import pylatex as pl
from pylatex.utils import italic, NoEscape
from pylatex.package import Package

doc = pl.Document('basic')
doc.packages.append(Package('booktabs'))
doc.packages.append(Package('multirow'))
doc.append(NoEscape('\n'.join(lines)))

doc.generate_pdf(os.path.splitext(opt.out_pdf_file)[0], clean_tex=True)

