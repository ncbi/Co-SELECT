#!/usr/bin/env python


import argparse
parser = argparse.ArgumentParser(description='Cluster TF experiments by the enriched shapemers.')
parser.add_argument("infile", default='', help="input file containing all enriched shapemers")
parser.add_argument("promiscuous_file", default='', help="input file containing all promiscuous shapemers")
parser.add_argument("heatmap_pdf", default='', help="output pdf heatmap")
parser.add_argument("pca_csv", default='', help="output csv containing PCA information")
opt = parser.parse_args()

print opt


import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

sns.set(color_codes=True)


promiscuous = pd.read_csv(opt.promiscuous_file)
promiscuous = promiscuous[promiscuous['en_th'] == 1.2]
promiscuous = promiscuous[promiscuous['promiscuity'] >= 0.33]
promiscuous = promiscuous[['shape', 'kmer']]
to_exclude = promiscuous.apply(lambda x: tuple(x), axis='columns')

df = pd.read_csv(opt.infile) #, header = [1,2,3], index_col=[0,1,2,3], low_memory=False)
df['value'] = 1
df = df.set_index(['family', 'tf', 'barcode', 'motif', 'threshold', 'shapemer','shape'])
df = df.unstack(fill_value=0).unstack(fill_value=0).unstack(fill_value=0)
df.columns = df.columns.droplevel()
print df.columns.get_level_values('threshold')
df = df.iloc[:, df.columns.get_level_values('threshold') == 1.2]
df.columns = df.columns.droplevel('threshold')

df2 = df.iloc[:, ~df.columns.isin(to_exclude)]

def get_cluster(df, fname, metric='hamming'):
  families = df.index.get_level_values('family')
  row_pal = sns.husl_palette(len(families), s=.75)
  row_lut = dict(zip(families, row_pal))
  row_colors = pd.DataFrame(df.index.map(lambda x: row_lut[x[0]]), index=df.index, columns=['family'])
  
  shapes = df.columns.get_level_values('shape')
  col_pal = sns.husl_palette(len(shapes), s=.75)
  col_lut = dict(zip(shapes, col_pal))
  col_colors = pd.DataFrame(df.columns.map(lambda x: col_lut[x[0]]), index=df.columns, columns=['shape'])
  
  g = sns.clustermap(df,
                     col_cluster = False,
                     metric = metric,
                     method = 'average',
                     figsize=(30, 20),
                     row_colors = row_colors,
                     col_colors = col_colors)
  
  g.savefig(fname)
  
get_cluster(df, opt.heatmap_pdf)
#get_cluster(df2, 'cluster_wo_promiscuous_actual_enrichment.pdf', 'euclidean')


from sklearn.decomposition import PCA

pca = PCA(n_components=3)
pca_result = pca.fit_transform(df2.values)

df3 = pd.DataFrame()

df3['pca.one'] = pca_result[:,0]
df3['pca.two'] = pca_result[:,1] 
df3['pca.three'] = pca_result[:,2]
df3['family'] = df2.index.map(lambda x: x[0])
df3['motif'] = df2.index.map(lambda x: x[3])
df3['subcategory'] = df2.index.map(lambda x: x[3] if x[3] in ['TAAA', 'TAAT', 'CGGAA', 'GGAAG', 'CACGTG', 'CATATG'] else 'Other')

df3.to_csv(opt.pca_csv, index=False)

print 'Explained variation per principal component: {}'.format(pca.explained_variance_ratio_)

