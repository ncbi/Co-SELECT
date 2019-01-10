#!/usr/bin/env python

import re, os, sys
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.realpath(__file__))+'/..'))
print sys.path

from analysis_scripts.tf_process import *
from analysis_scripts.shape_info import *
from analysis_scripts.tf_utils import *


orig_data_dir = '/panfs/pan1/aptax/FinalData'
top_data_dir = '/panfs/pan1/aptax/FinalData'
top_results_dir = '/panfs/pan1/dnashape/FinalResults/d0'
top_pwms_dir = 'SEQLOGO'

shape_type = 'MGW'
levels_type = 'publish'
lflank = 1
rflank = 1
threshold = '1.20'
MAX_SHAPEMERS = 10

tf = 'MAX'
family = 'bHLH'
primer = 'TGACCT20NGA'
motif = 'CACGTG'
cycle = 4


tf = 'PITX3'
family = 'homeodomain'
primer = 'TGCATC20NGA'
motif = 'TAAT'


import argparse
parser = argparse.ArgumentParser(description='Get the pwm of the sequence logo for a set of shapemers within the given contexts of a TF experiment')
parser.add_argument("tf",       help = "transcription factor name")
parser.add_argument("primer",   help = "barcode")
parser.add_argument("family",   help = "tf family")
parser.add_argument("motif",    help = "motif used by Co-SELECT")
parser.add_argument("cycle",    type = int, help = "SELEX round")
parser.add_argument("orig_dir", help = "top directory for original sequence files")
parser.add_argument("data_dir", help = "top directory for derived files")
parser.add_argument("res_dir",  help = "top results directory")
parser.add_argument("out_file", help = "file in which the output pwm would be saved")
opt = parser.parse_args()

tf = opt.tf
primer = opt.primer
family = opt.family
motif = opt.motif
cycle = opt.cycle
orig_data_dir = opt.orig_dir
top_data_dir = opt.data_dir
top_results_dir = opt.res_dir
out_file = opt.out_file


def debug(variable):
    print variable, '=', repr(eval(variable))


debug('tf')
debug('primer')
debug('family')
debug('motif')
debug('cycle')
debug('orig_data_dir')
debug('top_data_dir')
debug('top_results_dir')
debug('out_file')


info = TFInfo(tf, primer, family)
dist = len(motif) - 2

shape_levels_str = shape_info[levels_type][shape_type].getLevelsStr()
skip = shape_info[levels_type][shape_type].skip

seq_file = "%s/%s" % (orig_data_dir, info.getSequenceFile(cycle))
shape_file = "%s/%s" % (orig_data_dir, info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
count_file = "%s.cnt" % (seq_file)
fg_file = "%s/%s" % (top_data_dir, info.getContextFile(cycle, motif, dist, 'fg'))
bg_file = "%s/%s" % (top_data_dir, info.getContextFile(cycle, motif, dist, 'bg'))
parts_file = "%s/%s" % (top_data_dir, info.getFgPartsFile(cycle, motif, dist, lflank, rflank))

debug('seq_file')
debug('shape_file')
debug('count_file')
debug('fg_file')
debug('bg_file')
debug('parts_file')


enrich_file = '/'.join([top_results_dir, levels_type, threshold, info.family, info.tf, info.primer, '.'.join(['enriched', info.tf, info.primer, shape_type, str(cycle), motif, str(lflank), str(rflank), 'csv'])])

enr = pd.read_csv(enrich_file).sort_values(by='bg4.en.rank')
enr = enr[enr.label == 'both'].append(enr[enr.label == 'bg']).head(MAX_SHAPEMERS)
enr = enr[['kmer', 'label']].rename(index=str, columns={'label':'enrichment'})
print(enr)

promiscuous_file = '%s/highly_promiscuous_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
promis = pd.read_csv(promiscuous_file, dtype={'en_th':object})
promis = promis[(promis['shape'] == shape_type) & (promis['en_th'] == threshold)]
promis['promiscuity'] = 'high'
promis = promis[['kmer', 'promiscuity']]
print(promis)

combined = pd.merge(enr, promis, how='outer').rename(index=str, columns={'kmer':'shapemer'}) #['HHHHMH', 'HMHHXH', 'HMHHXX']
combined = combined.set_index(['shapemer'])
print(combined)

shapemers = combined.index.tolist()
debug('shapemers')



from itertools import izip
import pandas as pd
import numpy as np


shape_length = len(shapemers[0])
seq_length = shape_length + 4


def add_pwm(pwm, kmer, cnt):
  for i, c in enumerate(kmer):
    pwm['ACGT'.find(c), i] += cnt


def getNextSID(f):
  try:
    return int(f.readline().rstrip())
  except:
    return float('Inf')

def check_update_pwm_bg(shp, seq, cnt):
  if shp in shapemers:
      add_pwm(pwms.loc[shp, 'bg']['mat'], seq, cnt)

def check_update_pwm_fg(shp, seq, cnt):
  rev_seq = rev_comp(seq)
  rev_shp = shp[::-1]
  for s,p in [(seq, shp), (rev_seq, rev_shp)]:
    #print s, p
    if p in shapemers:
      #print 'adding', s, p
      add_pwm(pwms.loc[p, 'fg']['mat'], s, cnt)

#def print_bg_window(seq, shape, kmers, s, e, count, pwm, candidates):
#  wl = max(0, s-2)
#  wr = min(len(shape), e-2) 
#  window = shape[wl:wr]
#  #print s, e, window
#  for k in range(len(window)-shape_length+1):
#    check_update_pwm(window[k:k+shape_length], seq[wl+k:wl+k+shape_length+4], count, kmers, pwm, candidates, 'bg')
#
#
def print_shape_window(seq, shape, pos, motif_length, lflank, rflank, skip, count):
    # SEQ START POS = pos - lflank
    # SHAPE START POS = pos - lflank - 2 irrespective of shape feature
    # SEQ END POS = pos + rflank + k  where k is motif length
    # SHAPE END POS = pos + rflank + k - 2      for MGW,ProT and
    #               = pos + rflank + k - 2 + 1  for Roll, HelT
    wl = max(0, pos-lflank-2)
    wr = min(len(shape), pos+motif_length+rflank-skip) 
    window = shape[wl:wr]
    for k in range(len(window)-shape_length+1):
      #print window[k:k+shape_length], seq[wl+k:wl+k+shape_length+4], count
      check_update_pwm_fg(window[k:k+shape_length], seq[wl+k:wl+k+shape_length+4], count)

def ProcessFull(seq, shp, cnt):
  rev_seq = rev_comp(seq)
  rev_shp = shp[::-1]
  for s,p in [(seq, shp), (rev_seq, rev_shp)]:
    for i in range(len(p)-shape_length+1):
      check_update_pwm_bg(p[i:i+shape_length], s[i:i+shape_length+4], cnt)


#def ProcessBg(pwm, seq, shp, cnt, kmers):
#  rev_seq = complement(seq[::-1])
#  rev_shp = shp[::-1]
#  for s,p in [(seq,shp), (rev_seq, rev_shp)]:
#    for start, end in split_pos('[C]{'+str(shape_length+4)+',}', s):
#      print_bg_window(s, p, kmers, start, end, cnt, pwm, ['bg-longC'])
#
def ProcessFg(seq, shp, cnt, f, i, motif):
  words = f.readline().rstrip().split()
  assert(int(words[0]) == i)

  num_forward = int(words[1])
  num_reverse = int(words[2])
  positions = words[3].split(',')
  #print num_forward, num_reverse, positions
  for j in range(num_forward):
    i = int(positions[j])
    #print seq
    #print shp
    #print i, seq[i:i+len(motif)], motif
    assert(seq[i:i+len(motif)] == motif)
    print_shape_window(seq, shp, i, len(motif), lflank, rflank, skip, cnt)
  
  rev_seq = rev_comp(seq)
  rev_shp = shp[::-1]
  for j in range(num_forward, len(positions)):
    i = int(positions[j])
    #print rev_seq
    #print i, rev_seq[i:i+len(motif)], motif
    assert(rev_seq[i:i+len(motif)] == motif)
    print_shape_window(rev_seq, rev_shp, i, len(motif), lflank, rflank, skip, cnt)


#contexts = ['full', 'fgfull', 'fg', 'bg', 'bg-longC']
contexts = ['fg', 'bg']

pwms = pd.DataFrame(index = pd.MultiIndex.from_product([shapemers, contexts], names = ["shapemer", "context"]))
pwms['mat'] = pwms.index.map(lambda x: np.zeros((4, seq_length), dtype=np.int64))


 
with open(seq_file) as f, open(count_file) as g, open(shape_file) as h, open(fg_file) as fg, open(bg_file) as bg, open(parts_file) as parts:
  next_fg = getNextSID(fg)
  next_bg = getNextSID(bg)
  for i, (s,c,p) in enumerate(izip(f,g,h)):
    cnt = int(c)
    seq = s.rstrip()
    shp = p.rstrip()
    candidates = [] 
    if (i == next_fg):
      #candidates.append('fgfull')
      #print seq
      #print "  %s" % shp
      ProcessFg(seq, shp, cnt, parts, i, motif)
      next_fg = getNextSID(fg)
    if (i == next_bg):
      seq = seq[info.lbc_len : -info.rbc_len]
      shp = shp[info.lbc_len : -info.rbc_len]
      ProcessFull(seq, shp, cnt)
      next_bg = getNextSID(bg)
 
pwms['mat'].apply(lambda x: pd.DataFrame(x, index=pd.Index(list('ACGT'), name='base'), \
                                            columns=['X{}'.format(x+1) for x in range(seq_length)]).unstack()) \
                 .stack().join(combined).reset_index().to_csv(out_file, index=False)
