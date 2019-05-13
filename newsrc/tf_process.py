import re
import sys, os
sys.path.insert(0, os.path.dirname("../src/"))
from analysis_scripts.tf_utils import *
from analysis_scripts.kmer_utils import *

def get_bc(barcode):  
  res = re.match('([ACGT]+)\d+N([ACGT]+)', barcode)
  return res.group(1), res.group(2)

def get_bc_len(barcode):
  left, right = get_bc(barcode)
  return len(left), len(right)


def partition_aptamers(fg_type, seq_file, motif, nbr_file, fg_file, bg_file, barcode):
  distance_threshold = len(motif) - 2
  lbc_len, rbc_len = get_bc_len(barcode)
  if fg_type == 'd1enriched':
    d1_nbrs = pd.read_csv(nbr_file)['seqmer'].tolist()
  else:
    d1_nbrs = []
  with open(seq_file) as f, open(fg_file, "w") as g0, open(bg_file, "w") as g1:
    for i, l in enumerate(f):
       seq = l.rstrip()
       #seq = seq[2:-2]
       seq = seq[lbc_len : -rbc_len]
       rev = rev_comp(seq)
       d = hamming_all(seq, motif)
       rd = hamming_all(rev, motif)
       if (d == 0) or (rd == 0):
         print(i, file=g0)
       elif (d >= distance_threshold) and (rd >= distance_threshold):
         print(i, file=g1)
       else:
         if fg_type == 'd1all':
           if (d == 1) or (rd == 1):
             print(i, file=g0)
         elif fg_type == 'd1enriched':
           d = min([hamming_all(seq, nbr) for nbr in d1_nbrs])
           rd = min([hamming_all(rev, nbr) for nbr in d1_nbrs])
           if (d == 0) or (rd == 0):
             print(i, file=g0)

 
def gen_bg_seqmers(seq_file, count_file, sid_file, seqmer_file, seq_length, barcode):
  lbc_len, rbc_len = get_bc_len(barcode)
  with open(seqmer_file, 'w') as g:
    for (sid1, seq), (sid2, cnt) in zip(filterSeq(seq_file, sid_file), filterSeq(count_file, sid_file)):
      assert(sid1 == sid2)
      count = int(cnt)
      seq = seq[lbc_len : -rbc_len]
      rev = rev_comp(seq)
      for k in range(len(seq)-seq_length+1):
        rep = seq[k:k+seq_length]
        print(rep, count, file=g)
      for k in range(len(rev)-seq_length+1):
        rep = rev[k:k+seq_length]
        print(rep, count, file=g)

