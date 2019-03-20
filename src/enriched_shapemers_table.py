#from results_scripts.more_results import *
import pandas as pd

top_results_dir = '../results/d0'
ltype = 'publish'
cycle = 4
lflank = rflank = 1
fg_type = 'd0'

infile = '%s/%s/denriched_same.%d.l%d.r%d.csv' % (top_results_dir, ltype, cycle, lflank, rflank)
print(infile)

df = pd.read_csv(infile)
print(df)

df = df[(df['en_th'] == 1.2) & (df['label'] == 'both')]
print(df)
print(df.columns)



#outfile = '%s/enriched_shapemers_%s_cycle%d_%s.xlsx' % (top_results_dir, ltype, cycle, fg_type)
#csvfile = '%s/table_enriched_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, ltype, cycle, lflank, rflank)
outfile = 'tab_haha.xlsx'
csvfile = 'table_haha.csv'

tf_motif_file = './tf_coremotif.csv'


def save_dataframe(writer, sheet, df):
  df.dropna(axis='index', how='all', inplace=True)
  df.dropna(axis='columns', how='all', inplace=True)
  df.to_excel(writer, sheet_name = sheet)
  #for j, width in enumerate(get_col_widths(df)):
  #  writer.sheets[sheet].set_column(j, j, width)
  writer.sheets[sheet].freeze_panes(df.columns.nlevels+1, df.index.nlevels)


def get_enriched_shapemers_in_excel(combined, motif_file, xlsfile, csvfile):
  writer = pd.ExcelWriter(xlsfile, engine='xlsxwriter')
  tf_motif = pd.read_csv(motif_file)
  #combined = pd.read_csv(infile)
  shapes = combined['shape'].drop_duplicates()
  combined = combined[['family.x', 'tf.x', 'primer.x', 'kmer', 'label', 'shape']]
  combined['label'] = combined['label'].apply(lambda x: 'X' if x == 'both' else 'o')
  combined = combined.rename(index=str, columns={"family.x": "family", "primer.x": "primer", "tf.x": "tf", "kmer": "shapemer"})
  combined = combined.merge(tf_motif)
  combined = combined.rename(index=str, columns={"primer": "barcode"})
  tmp = combined[['tf', 'barcode']].drop_duplicates()
  tmp['exp'] = tmp.groupby('tf')['barcode'].transform(lambda x: [y+1 for y in range(len(x))] if (len(x) > 1) else [0])
  print tmp
  combined = combined.merge(tmp)
  print combined
  combined['tf'] = combined.apply(lambda x: "{} ({})".format(x['tf'], x['exp']) if x['exp'] else x['tf'], axis='columns')
  print combined
  print combined[combined['exp'] == 2]
  combined = combined.drop(['barcode', 'exp'], axis='columns')
  combined = combined.set_index(['family', 'tf', 'motif', 'shapemer','shape'])
  combined = combined.unstack().unstack()
  for i, shape in enumerate(shapes):
    df = combined.iloc[:, combined.columns.get_level_values('shape')==shape]
    save_dataframe(writer, shape, df)
  save_dataframe(writer, "all", combined)
  writer.save()
  return(combined)

df = get_enriched_shapemers_in_excel(df, tf_motif_file, outfile, csvfile)

df = df.fillna('')

df.columns = df.columns.droplevel()
#df.index.set_names('', inplace=True)


print df.index

#df.index.set_levels(['\\texttt{%s}' % x for x in df.index.levels[2].tolist()], level=2, inplace=True)

print df.index

cols = pd.DataFrame(index = df.columns).reset_index('shapemer').groupby('shape').agg({'shapemer':len}).reset_index()
csum = [0] + cols['shapemer'].cumsum().tolist()
counts = cols['shapemer'].tolist()
print cols
print cols

import re

#lines = df.to_latex(multicolumn=True, multirow=True, multicolumn_format='c', bold_rows=True)
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

#print lines

with open('hehehehe.tex', 'w') as f:
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

doc.generate_pdf(os.path.splitext('hehehehehehehe.pdf')[0], clean_tex=True)
