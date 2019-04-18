import sys, os, re
import pandas as pd
import subprocess as sp
from PyPDF2 import PdfFileMerger

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from params_config import *
from analysis_scripts.tf_process import *
from analysis_scripts.shape_info import *
from analysis_scripts.get_coverage import getCoverage
from round0_scripts.simple_predict import *

DOIT_CONFIG = {'check_file_uptodate': 'timestamp', 'verbosity':2}

class TaskInfo:
  def __init__(self, tf, primer, family, accessions, motifs, cycles, distance_thresholds):
    self.tf = tf
    self.primer = primer
    self.family = family
    self.accessions = accessions
    self.motifs = motifs
    self.cycles = cycles
    self.shape_length = 6  # len(motif) + 2
    self.datadir = '%s/%s/%s' % (top_data_dir, tf, primer)
    self.tf_info = TFInfo(tf, primer, family)
    self.distances = distance_thresholds


tf_info   = pd.read_csv(tf_info_file)
tf_motif  = pd.read_csv(tf_motif_file)
tf_run    = pd.read_csv(tf_run_coselect_file)
accession = pd.read_csv(nonzero_accession_file)

tfs = tf_info.merge(tf_motif).merge(tf_run)
tfs['distance'] = tfs['motif'].str.len() - 2
tfs['accessions'] = ''
for cycle in cycles:
  tfs['cycle'] = cycle
  tfs['accessions'] = tfs['accessions'] + ',' + tfs['cycle'].map(lambda x:str(x)) + ':' + tfs.merge(accession)['accession']
tfs['accessions'] = tfs['accessions'].map(lambda x: {int(y.split(':')[0]):y.split(':')[1] for y in x.split(',')[1:]})


motif_dist = tfs[['motif', 'distance']].drop_duplicates()
motifs = motif_dist['motif'].tolist()
distances = motif_dist['distance'].tolist()


def combine_csvs(configs, outfile):
  index_cols = list(set(configs.columns) - set(['infile']))
  df = pd.DataFrame()
  for indx, r in configs.iterrows():
    tmp = pd.read_csv(r['infile'])
    for col in index_cols:
      tmp[col] = r[col]
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def merge_pdfs(infiles, outfile):
  merger = PdfFileMerger()
  for pdf in infiles:
    merger.append(open(pdf, 'rb'))
  with open(outfile, 'wb') as fout:
    merger.write(fout)


def merge_pdfs_exclude_empty(infiles, outfile):
  keepfiles = filter(lambda x: sp.check_output('pdfinfo %s | grep Title:' % x, shell=True).find("no logo") < 0, infiles)
  merge_pdfs(keepfiles, outfile)


def merge_pdfs_from_list(inlist, outfile):
  infiles = pd.read_csv(inlist)['pdf']
  merge_pdfs(infiles, outfile)


def combine_data_frames(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (tf1, primer1, family1, tf2, primer2, family2) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['tf.x'] = tf1
    tmp['primer.x'] = primer1
    tmp['family.x'] = family1
    tmp['tf.y'] = tf2
    tmp['primer.y'] = primer2
    tmp['family.y'] = family2
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def combine_qvalues(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (en_th, shape_type, levels_type, family, ctx) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['en_th'] = en_th
    tmp['shape'] = shape_type
    tmp['levels_type'] = levels_type
    tmp['family'] = family
    tmp['ctx'] = ctx
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def combine_enriched(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (en_th, shape_type) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['en_th'] = en_th
    tmp['shape'] = shape_type
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def combine_detailed_results(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (en_th, shape_type) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['en_th'] = en_th
    tmp['shape'] = shape_type
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def combine_detailed_results2(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, ltype in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['discretization'] = ltype
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def combine_data_frames_th_shape(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (en_th, shape_type) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['en_th'] = en_th
    tmp['shape'] = shape_type
    df = df.append(tmp, ignore_index=True)
  print df
  df.to_csv(outfile, index=False)

