import pandas as pd

def get_col_widths(dataframe):
    # Then, we concatenate this to the max of the lengths of column name and its values for each column, left to right
    return [max([len(str(s))+3 for s in dataframe[col].values] + [len(col)]) for col in dataframe.columns]

def set_width(writer, sheet, df):
  for j, width in enumerate(get_col_widths(df)):
    writer.sheets[sheet].set_column(j, j, width)


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
  combined = combined.unstack(fill_value='').unstack(fill_value='').unstack(fill_value='')
  for i, shape in enumerate(shapes):
    df = combined.iloc[:, combined.columns.get_level_values('shape')==shape]
    df.to_excel(writer, sheet_name=shape)
    res_pitx3 =  df[df.index.get_level_values('tf').isin(['PITX3', 'GBX1'])]
    if not res_pitx3.empty:
      res_pitx3.to_excel(writer, sheet_name="%s-PITX3" % shape)
  writer.save()



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
    df.to_excel(writer, sheet_name=shape)
  combined.to_excel(writer, sheet_name='AllShapes')
  writer.save()

