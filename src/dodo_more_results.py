from dodo_common import *

level_names = {'main-text': 'publish', 'alternative': 'other1'}

def get_col_widths(dataframe):
    # Then, we concatenate this to the max of the lengths of column name and its values for each column, left to right
    return [max([len(str(s))+3 for s in dataframe[col].values] + [len(col)]) for col in dataframe.columns]

def set_width(writer, sheet, df):
  for j, width in enumerate(get_col_widths(df)):
    writer.sheets[sheet].set_column(j, j, width)

def save_data(writer, sheet, df):
  df = df.set_index(['family', 'tf', 'barcode', 'motif', 'threshold', 'shapemer'])
  df = df.unstack(fill_value='').unstack()
  df.to_excel(writer, sheet_name=sheet)
  #set_width(shape,res)


def save_all(writer, sheet, df):
  df = df.set_index(['family', 'tf', 'barcode', 'motif', 'threshold', 'shapemer', 'shape'])
  df = df.unstack().unstack(fill_value='').unstack()
  df.to_excel(writer, sheet_name=sheet)
  #set_width(shape,res)

def save_detailed(writer, sheet, df):
  df = df.set_index(['family', 'tf', 'barcode', 'motif', 'discheme', 'threshold'])
  df = df.unstack()
  df.to_excel(writer, sheet_name=sheet)


def save_detailed_all(writer, sheet, df):
  df = df.set_index(['family', 'tf', 'barcode', 'motif', 'discheme', 'threshold', 'shape'])
  df = df.unstack().unstack()
  df.to_excel(writer, sheet_name=sheet)


def get_enriched_shapemers(outfile, ltype):
  writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
  tf_motif = pd.read_csv('tf_coremotif.csv')
  #combined = pd.DataFrame()
  for i, shape in enumerate(shapes):
    res = pd.DataFrame()
    #for ltype in level_names.keys():
    for th in en_thresholds:
        infile = "%s/%s/%s/denriched_same.%s.l1.r1.csv" % (top_results_dir, ltype, th, shape)
        df = pd.read_csv(infile)
        df = df[['family.x', 'tf.x', 'primer.x', 'kmer', 'label']]
        df = df[df['label'].isin(['both','bg'])]
        df = df.rename(index=str, columns={"family.x": "family", "primer.x": "primer", "tf.x": "tf", "kmer": "shapemer"})
        df = df.merge(tf_motif)
        df = df.rename(index=str, columns={"primer": "barcode"})
        df = df[['family', 'tf', 'barcode', 'motif', 'shapemer', 'label']]
        #df['discheme'] = ltype
        df['threshold'] = th
        df['label'] = df['label'].apply(lambda x: 'X' if x == 'both' else 'o')
        res = res.append(df)
    save_data(writer, shape, res)
    res_pitx3 =  res[res.tf.isin(['PITX3', 'GBX1'])]
    if not res_pitx3.empty:
      save_data(writer, "%s-PITX3" % shape, res_pitx3)
    res = res[['family', 'tf', 'barcode', 'motif', 'threshold', 'shapemer', 'label']]
    res['shape'] = shape
    #combined = combined.append(res)
  
  #save_all(writer, 'AllShapes',combined)
  #writer.save()



def get_detailed_results(outfile):
  writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
  tf_motif = pd.read_csv('tf_coremotif.csv')
  combined = pd.DataFrame()
  for i, shape in enumerate(shapes):
    res = pd.DataFrame()
    for ltype in level_names.keys():
      for th in en_thresholds:
        infile = "%s/%s/%s/dfisher_same.%s.l1.r1.csv" % (top_results_dir, level_names[ltype], th, shape)
        df = pd.read_csv(infile)
        df = df[['family.x', 'tf.x', 'primer.x', 'pvalue', 'FSBS',  'FSBN', 'FNBS', 'FNBN']]
        df = df.rename(index=str, columns={"family.x": "family", "primer.x": "primer", "tf.x": "tf"})
        df = df.merge(tf_motif)
        df = df.rename(index=str, columns={"primer": "barcode"})
        df = df[['family', 'tf', 'barcode', 'motif', 'pvalue', 'FSBS',  'FSBN', 'FNBS', 'FNBN']]
        df['discheme'] = ltype
        df['threshold'] = th
        res = res.append(df)
    save_detailed(writer, shape, res)
    res = res[['family', 'tf', 'barcode', 'motif', 'pvalue', 'discheme', 'threshold']]
    res = res.round({'pvalue': 3})
    res['shape'] = shape
    combined = combined.append(res)
  save_detailed_all(writer, 'AllShapes',combined)


def task_get_enriched_shapemers():
  for lflank, rflank in flank_configs:
    for ltype in level_names.keys():
      infiles = ["%s/%s/%s/denriched_same.%s.l1.r1.csv" % (top_results_dir, level_names[ltype], en_th, shape_type) for en_th in en_thresholds for shape_type in shapes]
      print(infiles)
      outfile = '%s/enriched_shapemers_%s_l%d.r%d.xlsx' % (top_results_dir, ltype, lflank, rflank)
      yield {
        'name'      : outfile,
        'actions'   : [(get_enriched_shapemers, [outfile, level_names[ltype]])],
        'file_dep'  : infiles,
        'targets'   : [outfile],
        'clean'     : True,
      }


def task_get_detailed_results():
  for lflank, rflank in flank_configs:
    infiles = ["%s/%s/%s/dfisher_same.%s.l1.r1.csv" % (top_results_dir, levels_type, en_th, shape_type) for en_th in en_thresholds for shape_type in shapes for levels_type in discrete_levels_type]
    outfile = '%s/detailed_results_l%d.r%d.xlsx' % (top_results_dir, lflank, rflank)
    yield {
      'name'      : outfile,
      'actions'   : [(get_detailed_results, [outfile])],
      'file_dep'  : infiles,
      'targets'   : [outfile],
      'clean'     : True,
    }
