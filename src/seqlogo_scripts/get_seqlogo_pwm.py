#!/usr/bin/env python

import re, os, sys
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.realpath(__file__))+'/..'))
print sys.path

from analysis_scripts.tf_process import *
from analysis_scripts.shape_info import *
from analysis_scripts.tf_utils import *


orig_data_dir = '/panfs/pan1/aptax/FinalData'
top_data_dir = '/panfs/pan1/aptax/FinalData'
top_pwms_dir = 'SEQLOGO'

shape_type = 'MGW'
levels_type = 'publish'
lflank = 1
rflank = 1
skip = 2

tf = 'MAX'
family = 'bHLH'
primer = 'TGACCT20NGA'
motif = 'CACGTG'
cycle = 4


tf = 'PITX3'
family = 'homeodomain'
primer = 'TGCATC20NGA'
motif = 'TAAT'

info = TFInfo(tf, primer, family)
dist = len(motif) - 2

shape_levels_str = shape_info[levels_type][shape_type].getLevelsStr()

seq_file = "%s/%s" % (orig_data_dir, info.getSequenceFile(cycle))
shape_file = "%s/%s" % (orig_data_dir, info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
count_file = "%s.cnt" % (seq_file)
fg_file = "%s/%s" % (top_data_dir, info.getContextFile(cycle, motif, dist, 'fg'))
bg_file = "%s/%s" % (top_data_dir, info.getContextFile(cycle, motif, dist, 'bg'))
parts_file = "%s/%s" % (top_data_dir, info.getFgPartsFile(cycle, motif, dist, lflank, rflank))

shapemers = ['HHHHMH', 'HMHHXH', 'HMHHXX']

out_file = top_pwms_dir + '/seqlogo_info_%s_%s_NEW.csv' % (tf, primer)

def debug(variable):
    print variable, '=', repr(eval(variable))

#import argparse
#parser = argparse.ArgumentParser(description='Get the pwm of the sequence logo for a set of shapemers within the given contexts of a TF experiment')
#parser.add_argument("seq_file",   default = seq_file,   help = "input sequence file")
#parser.add_argument("shape_file", default = shape_file, help = "input shape file")
#parser.add_argument("count_file", default = count_file, help = "sequence count file")
#parser.add_argument("fg_file",    default = fg_file,    help = "file containing ids of sequences containing motif or partial motifs")
#parser.add_argument("bg_file",    default = bg_file,    help = "file containing ids of motif-free sequences")
#parser.add_argument("parts_file", default = parts_file, help = "file containing segments of motif or partial motifs")
#parser.add_argument("motif",      default = motif,      help = "sequence motif")
#parser.add_argument("primer",     default = primer,     help = "primer used in the SELEX experiment")
#parser.add_argument("shapemers",  default = ','.join(shapemers), help="comma separated list of shapemers")
#parser.add_argument("out_file",   default = out_file,    help = "output file containing pwm of all shapemers")
#opt = parser.parse_args()
#
#seq_file = opt.seq_file
#shape_file = opt.shape_file
#count_file = opt.count_file
#fg_file = opt.fg_file
#bg_file = opt.bg_file
#parts_file = opt.parts_file
#shapemers = opt.shapemers.split(',')
#out_file = opt.out_file

debug('seq_file')
debug('shape_file')
debug('count_file')
debug('fg_file')
debug('bg_file')
debug('parts_file')
debug('shapemers')
debug('out_file')


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
      assert(len(seq)== 20)
      shp = shp[info.lbc_len : -info.rbc_len]
      ProcessFull(seq, shp, cnt)
      next_bg = getNextSID(bg)
 
pwms['mat'].apply(lambda x: pd.DataFrame(x, index=pd.Index(list('ACGT'), name='base'), \
                                            columns=['X{}'.format(x+1) for x in range(seq_length)]).unstack()) \
                 .stack().reset_index().to_csv(out_file, index=False)
