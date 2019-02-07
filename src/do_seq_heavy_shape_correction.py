import pandas as pd

df = pd.read_csv('../results/d0/publish/mic_correction_same.4.l1.r1.csv')

print(df)

fisher = pd.read_csv('../results/d0/publish/combined_fisher.4.l1.r1.csv')

print(fisher)
fisher = fisher[['en_th', 'shape', 'family.x', 'tf.x', 'primer.x', 'FSBN','FSBS', 'FNBN', 'FNBS', 'pvalue']]
fisher = fisher.rename(index=str, columns={'family.x': 'family', 'tf.x':'tf', 'primer.x':'primer', 'pvalue':'oldpvalue'})

print(fisher)

df = df.merge(fisher)
print(df)

import scipy.stats as stats

df['pvalue'] = df.apply(lambda x: stats.fisher_exact([[x.FSBS-x.FSBSIgnore, x.FSBN], [x.FNBS-x.FNBSIgnore, x.FNBN]], alternative='greater')[1], axis='columns')


import statsmodels.stats.multitest as smm

df['rvalue'] = df.groupby(['en_th', 'family', 'shape'])['pvalue'].transform(lambda x: smm.multipletests(x, method='fdr_bh')[1])

print(df)


def get_col_widths(dataframe):
  # First we find the maximum length of the index columns   
  idx_max = [max([len(str(s)) for s in dataframe.index.get_level_values(idx)] + [len(str(idx))]) for idx in dataframe.index.names]
  # Then, we concatenate this to the max of the lengths of column name and its values for each column, left to right
  return idx_max + [max([len(str(s)) for s in dataframe[col].values] + \
                        [len(str(x)) for x in col] if dataframe.columns.nlevels > 1 else [len(str(col))]) for col in dataframe.columns]

def save_dataframe(writer, sheet, df):
  df.to_excel(writer, sheet_name = sheet)
  for j, width in enumerate(get_col_widths(df)):
    writer.sheets[sheet].set_column(j, j, width)
  writer.sheets[sheet].freeze_panes(df.columns.nlevels+1, df.index.nlevels)


writer = pd.ExcelWriter('detailed.xlsx', engine='xlsxwriter')
combined = df.copy()
shapes = combined['shape'].drop_duplicates()
combined = combined[['family', 'tf', 'primer', 'oldpvalue', 'pvalue', 'FSBS',  'FSBN', 'FNBS', 'FNBN', 'FSBSIgnore', 'FNBSIgnore', 'en_th', 'shape']]
combined = combined.rename(index=str, columns={'en_th':'threshold'})
combined = combined.rename(index=str, columns={"primer": "barcode"})
combined = combined.round({'pvalue': 3, 'oldpvalue': 3})
combined = combined.set_index(['family', 'tf', 'barcode', 'threshold','shape'])
combined = combined.unstack().unstack()
for i, shape in enumerate(shapes):
  tdf = combined.iloc[:, combined.columns.get_level_values('shape')==shape]
  save_dataframe(writer, shape, tdf)
save_dataframe(writer, 'AllShapes', combined)
writer.save()




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

with open('haha.tex', 'w') as f:
    f.write('\n'.join(lines))

import pylatex as pl
from pylatex.utils import italic, NoEscape
from pylatex.package import Package
import os

doc = pl.Document('basic')
doc.packages.append(Package('booktabs'))
doc.packages.append(Package('multirow'))
doc.append(NoEscape('\n'.join(lines)))

doc.generate_pdf(os.path.splitext('haha.tex')[0], clean_tex=True)

