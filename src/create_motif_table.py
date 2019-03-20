import pandas as pd
inventory_dir = '.'

tfs1 = pd.read_csv(inventory_dir + '/tf_run_coselect.csv')
tfs2 = pd.read_csv(inventory_dir + '/tf_coremotif.csv')
tfs3 = pd.read_csv(inventory_dir + '/tf_inventory_jolma_ronshamir.csv')

tfs = tfs1.merge(tfs2).merge(tfs3)

tfs = tfs[['family', 'tf', 'primer', 'motif']]

#tfs = tfs.set_index(['family', 'tf', 'primer'])
tfs = tfs.rename(index=str, columns={"primer": "barcode"})
tmp = tfs[['tf', 'barcode']].drop_duplicates()
tmp['exp'] = tmp.groupby('tf')['barcode'].transform(lambda x: [y+1 for y in range(len(x))] if (len(x) > 1) else [0])
print tmp
tfs = tfs.merge(tmp)
print tfs
tfs['tf'] = tfs.apply(lambda x: "{} ({})".format(x['tf'], x['exp']) if x['exp'] else x['tf'], axis='columns')
print tfs
print tfs[tfs['exp'] == 2]

tfs["tf.lower"] = tfs["tf"].str.lower()
tfs.sort_values(by=['family', 'tf.lower'], inplace=True)

tfs = tfs.reset_index()
print tfs
tfs.index += 1
print tfs
tfs = tfs.drop(['barcode', 'exp', 'tf.lower', 'index'], axis='columns')

tfs.to_latex('coremotifs_to_tomodify.tex', longtable=True)


