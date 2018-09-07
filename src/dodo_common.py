import sys
import os, re
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from params_config import *
from analysis_scripts.tf_process import *
from analysis_scripts.shape_info import *
from analysis_scripts.get_coverage import getCoverage
from round0_scripts.simple_predict import *

DOIT_CONFIG = {'check_file_uptodate': 'timestamp', 'verbosity':2}

class TaskInfo:
  def __init__(self, tf, primer, family, accession, motifs, cycles, distance_thresholds):
    self.tf = tf
    self.primer = primer
    self.family = family
    self.accession = accession
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
tfs['final'] = 4
tfs = tfs.merge(accession, left_on=['tf', 'primer', 'final'], right_on=['tf', 'primer', 'cycle'])

print tfs

motif_dist = tfs[['motif', 'distance']].drop_duplicates()
motifs = motif_dist['motif'].tolist()
distances = motif_dist['distance'].tolist()


