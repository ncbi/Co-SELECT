import pandas as pd

def get_col_widths(dataframe):
  # First we find the maximum length of the index columns   
  idx_max = [max([len(str(s)) for s in dataframe.index.get_level_values(idx)] + [len(str(idx))]) for idx in dataframe.index.names]
  # Then, we concatenate this to the max of the lengths of column name and its values for each column, left to right
  return idx_max + [max([len(str(s)) for s in dataframe[col].values] + \
                        [len(str(x)) for x in col] if dataframe.columns.nlevels > 1 else [len(str(col))]) for col in dataframe.columns]

def save_dataframe(writer, sheet, df):
  df.dropna(axis='index', how='all', inplace=True)
  df.dropna(axis='columns', how='all', inplace=True)
  df.to_excel(writer, sheet_name = sheet)
  for j, width in enumerate(get_col_widths(df)):
    writer.sheets[sheet].set_column(j, j, width)
  writer.sheets[sheet].freeze_panes(df.columns.nlevels+1, df.index.nlevels)


def get_enriched_shapemers_in_excel(infile, motif_file, xlsfile, csvfile):
  writer = pd.ExcelWriter(xlsfile, engine='xlsxwriter')
  tf_motif = pd.read_csv(motif_file)
  combined = pd.read_csv(infile)
  shapes = combined['shape'].drop_duplicates()
  combined = combined[['family.x', 'tf.x', 'primer.x', 'kmer', 'label', 'en_th', 'shape']]
  combined = combined[combined['label'].isin(['both','bg'])]
  combined['label'] = combined['label'].apply(lambda x: 'X' if x == 'both' else 'o')
  combined = combined.rename(index=str, columns={"family.x": "family", "primer.x": "primer", "tf.x": "tf", "kmer": "shapemer", 'en_th':'threshold'})
  combined = combined.merge(tf_motif)
  combined = combined.rename(index=str, columns={"primer": "barcode"})
  combined.drop('label', axis=1).to_csv(csvfile, index=False)
  combined = combined.set_index(['family', 'tf', 'barcode', 'motif', 'threshold', 'shapemer','shape'])
  combined = combined.unstack().unstack().unstack()
  for i, shape in enumerate(shapes):
    df = combined.iloc[:, combined.columns.get_level_values('shape')==shape]
    save_dataframe(writer, shape, df)
    res_pitx3 =  df[df.index.get_level_values('tf').isin(['PITX3', 'GBX1'])]
    if not res_pitx3.empty:
      sheet = "%s-PITX3" % shape
      save_dataframe(writer, sheet, res_pitx3)
  writer.save()


def get_enriched_shapemers_in_dataframe(combined, motif_file):
  tf_motif = pd.read_csv(motif_file)
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
  return(combined)


def get_detailed_results_in_excel(infile, motif_file, outfile, rename_dict):
  writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
  tf_motif = pd.read_csv(motif_file)
  combined = pd.read_csv(infile)
  shapes = combined['shape'].drop_duplicates()
  combined = combined[['family.x', 'tf.x', 'primer.x', 'pvalue', 'FSBS',  'FSBN', 'FNBS', 'FNBN', 'discretization', 'en_th', 'shape']]
  combined = combined.rename(index=str, columns={"family.x": "family", "primer.x": "primer", "tf.x": "tf", 'en_th':'threshold'})
  combined = combined.merge(tf_motif)
  combined = combined.rename(index=str, columns={"primer": "barcode"})
  combined = combined.round({'pvalue': 3})
  combined['discretization'] = combined['discretization'].map(rename_dict)
  combined = combined.set_index(['family', 'tf', 'barcode', 'motif', 'discretization', 'threshold','shape'])
  combined = combined.unstack().unstack()
  for i, shape in enumerate(shapes):
    df = combined.iloc[:, combined.columns.get_level_values('shape')==shape]
    save_dataframe(writer, shape, df)
  save_dataframe(writer, 'AllShapes', combined)
  writer.save()

