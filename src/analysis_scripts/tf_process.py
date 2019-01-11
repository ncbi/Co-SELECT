from itertools import product, izip
from tf_utils import *
from kmer_utils import *
import pandas as pd

def getRepresent(shape):
  rev = shape[::-1]
  return shape if (shape < rev) else rev

class TFInfo:
  def __init__(self, tf, primer, family):
    self.tf = tf
    self.primer = primer
    self.family = family
    res = re.match('([ACGT]+)\d+N([ACGT]+)', primer)
    self.lbc_len = len(res.group(1))
    self.rbc_len = len(res.group(2))
    
  def getSequenceFile(self, cycle):
    return "%s/%s/%s/%d.non" % (self.family if cycle else 'NoFamily', self.tf if cycle else 'ZeroCycle', self.primer, cycle)

  def getDiscreteShapeFile(self, cycle, shape_type, shape_levels_str):
    return "%s.%s.%s" % (self.getSequenceFile(cycle), shape_type, shape_levels_str)
  
  
  def getContextFile(self, cycle, motif, distance_threshold, ctx):
    fname = "%s.%s" % (self.getSequenceFile(cycle), motif)
    fname += ".d%d" % (distance_threshold)
    fname += ".%s" % (ctx)
    return fname
  
  def getContextPrefix(self, cycle, motif, distance_threshold, lflank, rflank, ctx):
    fname = self.getContextFile(cycle, motif, distance_threshold, ctx)
    if ctx == 'fg': 
      fname += ".l%d.r%d" % (lflank, rflank)
    return fname

  def getFgPartsFile(self, cycle, motif, distance_threshold, lflank, rflank):
    return self.getContextPrefix(cycle, motif, distance_threshold, lflank, rflank, 'fg') + ".part"
  
  def getContextedShapemerFile(self, cycle, motif, dist_th, lflank, rflank, ctx, shape, shape_levels_str):
    return self.getContextPrefix(cycle, motif, dist_th, lflank, rflank, ctx) + '.' + shape + '.' + shape_levels_str + '.mer'
  
  def getContextedShapemerCountFile(self, cycle, motif, dist_th, lflank, rflank, ctx, shape, shape_levels_str):
    return self.getContextedShapemerFile(cycle, motif, dist_th, lflank, rflank, ctx, shape, shape_levels_str) + '.cnt'
 

  def partition_aptamers_all_d1(self, seq_file, motif, distance_threshold, fg_file, bg_file):
    with open(seq_file) as f, open(fg_file, "w") as g0, open(bg_file, "w") as g1:
      for i, l in enumerate(f):
         seq = l.rstrip()
         #seq = seq[2:-2]
         seq = seq[self.lbc_len : -self.rbc_len]
         rev = rev_comp(seq)
         d = hamming_all(seq, motif)
         rd = hamming_all(rev, motif)
         if (d <= 1) or (rd <= 1):
           print >>g0, i
         if (d >= distance_threshold) and (rd >= distance_threshold):
           print >>g1, i
 
  def partition_aptamers(self, fg_type, seq_file, motif, distance_threshold, nbr_file, fg_file, bg_file):
    if fg_type == 'd1enriched':
      d1_nbrs = pd.read_csv(nbr_file)['seqmer'].tolist()
    else:
      d1_nbrs = []
    with open(seq_file) as f, open(fg_file, "w") as g0, open(bg_file, "w") as g1:
      for i, l in enumerate(f):
         seq = l.rstrip()
         #seq = seq[2:-2]
         seq = seq[self.lbc_len : -self.rbc_len]
         rev = rev_comp(seq)
         d = hamming_all(seq, motif)
         rd = hamming_all(rev, motif)
         if (d == 0) or (rd == 0):
           print >>g0, i
         elif (d >= distance_threshold) and (rd >= distance_threshold):
           print >>g1, i
         else:
           if fg_type == 'd1all':
             if (d == 1) or (rd == 1):
               print >>g0, i
           elif fg_type == 'd1enriched':
             d = min([hamming_all(seq, nbr) for nbr in d1_nbrs])
             rd = min([hamming_all(rev, nbr) for nbr in d1_nbrs])
             if (d == 0) or (rd == 0):
               print >>g0, i


  def partition_aptamers_orig(self, seq_file, motif, distance_threshold, fg_file, bg_file):
    with open(seq_file) as f, open(fg_file, "w") as g0, open(bg_file, "w") as g1:
      for i, l in enumerate(f):
         seq = l.rstrip()
         #seq = seq[2:-2]
         seq = seq[self.lbc_len : -self.rbc_len]
         rev = rev_comp(seq)
         d = hamming_all(seq, motif)
         rd = hamming_all(rev, motif)
         if (d == 0) or (rd == 0):
           print >>g0, i
         if (d >= distance_threshold) and (rd >= distance_threshold):
           print >>g1, i

  def partition_aptamers_matching_pair(self, seq_file, motif, distance_threshold, fg_file, bg_file):
    with open(seq_file) as f, open(fg_file, "w") as g0, open(bg_file, "w") as g1:
      for i, l in enumerate(f):
         seq = l.rstrip()
         #seq = seq[2:-2]       # WE NO LONGER DISCARD THE TWO BASES AT THE ENDS 
         seq = seq[self.lbc_len : -self.rbc_len]
         rev = rev_comp(seq)
         d = hamming_all(seq, motif)
         rd = hamming_all(rev, motif)
         if (d == 0) or (rd == 0):
           print >>g0, i
         else:
           #if (d >= distance_threshold) and (rd >= distance_threshold):
           #  print >>g1, i
           #d = hamming_pair_match_with_barcode(seq, motif, self.lbc_len, self.rbc_len)
           #rd = hamming_pair_match_with_barcode(rev, motif, self.rbc_len, self.lbc_len)
           d = hamming_pair_match(seq, motif)
           rd = hamming_pair_match(rev, motif)
           if (d <= 1) and (rd <= 1):
             print >>g1, i
           
  def gen_fg_parts_orig(self, seq_file, sid_file, motif, parts_file):
    with open(parts_file, "w") as g:
      for i, seq in filterSeq(seq_file, sid_file):
         #seq = seq[2:-2]
         seq = seq[self.lbc_len : -self.rbc_len]
         rev = rev_comp(seq)
         pos = find_all_occur(seq, motif)
         rpos = find_all_occur(rev, motif)
         pos = [x + self.lbc_len for x in pos]
         rpos = [x + self.rbc_len for x in rpos]
         print >>g, i, len(pos), len(rpos), ",".join([str(x) for x in pos+rpos])
 
  def gen_fg_parts_all_d1_hamming_nbr(self, seq_file, sid_file, motif, parts_file):
    with open(parts_file, "w") as g:
      for i, seq in filterSeq(seq_file, sid_file):
         #seq = seq[2:-2]
         seq = seq[self.lbc_len : -self.rbc_len]
         rev = rev_comp(seq)
         pos = find_all_occur_hamming_nbr(seq, motif, 1)
         rpos = find_all_occur_hamming_nbr(rev, motif, 1)
         pos = [x + self.lbc_len for x in pos]
         rpos = [x + self.rbc_len for x in rpos]
         print >>g, i, len(pos), len(rpos), ",".join([str(x) for x in pos+rpos])

 
  def gen_fg_parts(self, fg_type, seq_file, sid_file, nbr_file, motif, parts_file):
    print fg_type
    if fg_type == 'd1enriched':
      d1_nbrs = pd.read_csv(nbr_file)['seqmer'].tolist()
      d1_nbrs = list(set(d1_nbrs) - set([motif]))
    else:
      d1_nbrs = []
    with open(parts_file, "w") as g:
      for i, seq in filterSeq(seq_file, sid_file):
         #seq = seq[2:-2]
         seq = seq[self.lbc_len : -self.rbc_len]
         rev = rev_comp(seq)
         if fg_type == 'd1all':
           pos = find_all_occur_hamming_nbr(seq, motif, 1)
           rpos = find_all_occur_hamming_nbr(rev, motif, 1)
         else: # both for d0 and d1enriched
           pos = find_all_occur(seq, motif)
           rpos = find_all_occur(rev, motif)
         if fg_type == 'd1enriched':
           for x in d1_nbrs:
             pos += find_all_occur(seq, x)
             rpos += find_all_occur(rev, x)
         pos = [x + self.lbc_len for x in sorted(pos)]
         rpos = [x + self.rbc_len for x in sorted(rpos)]
         print >>g, i, len(pos), len(rpos), ",".join([str(x) for x in pos+rpos])


  def discretize_shape(self, shape_info, shape_file, discrete_file):
    with open(discrete_file, "w") as g:
      for shape_str in shape_iter(shape_info, shape_file):
        # shape_str is already chomped for 'skip' bases at both end
        g.write(shape_str + '\n')

  def print_shape_window(self, shape, seq_id, pos, motif_length, shape_length, lflank, rflank, skip, count, f):
    # SEQ START POS = pos - lflank
    # SHAPE START POS = pos - lflank - 2 irrespective of shape feature
    # SEQ END POS = pos + rflank + k  where k is motif length
    # SHAPE END POS = pos + rflank + k - 2      for MGW,ProT and
    #               = pos + rflank + k - 2 + 1  for Roll, HelT
    wl = max(0, pos-lflank-2)
    wr = min(len(shape), pos+motif_length+rflank-skip) 
    window = shape[wl:wr]
    #print window
    for k in range(len(window)-shape_length+1):
      #print window[k:k+shape_length], getRepresent(window[k:k+shape_length])
      rep = getRepresent(window[k:k+shape_length])
      for j in range(count):
        print >>f, rep, seq_id
      #print rep, seq_id

 
  def gen_fg_shapemers(self, shape_info, shape_length, seq_file, shape_file, count_file, sid_file, motif, lflank, rflank, parts_file, shapemer_file):
    with open(parts_file) as f, open(shapemer_file, 'w') as g:
      for l, (sid0, seq), (sid1, shape), (sid2, cnt) in izip(f, filterSeq(seq_file, sid_file), filterSeq(shape_file, sid_file), filterSeq(count_file, sid_file)):
        #print seq
        #print "%s%s" % (''.join([' ']*shape_info.skip), shape)
        words = l.rstrip().split(' ')
        seq_id = int(words[0])
        assert(seq_id == sid0)
        assert(seq_id == sid1)
        assert(seq_id == sid2)
        count = int(cnt)
        num_forward = int(words[1])
        num_reverse = int(words[2])
        positions = words[3].split(',')
        for i in range(num_forward):
          pos = int(positions[i])
          #print '->', pos, seq[pos:pos+len(motif)]
          assert(hamming_single(seq[pos:pos+len(motif)], motif) <= 1)
          self.print_shape_window(shape, seq_id, pos, len(motif), shape_length, lflank, rflank, shape_info.skip, count, g)
  
        shape = shape[::-1]
        rev = rev_comp(seq)
        for i in range(num_forward, len(positions)):
          pos = int(positions[i])
          #print '<-', pos, rev[pos:pos+len(motif)]
          assert(hamming_single(rev[pos:pos+len(motif)], motif) <= 1)
          self.print_shape_window(shape, seq_id, pos, len(motif), shape_length, lflank, rflank, shape_info.skip, count, g)
 
  def gen_bg_shapemers(self, shape_info, shape_length, shape_file, count_file, sid_file, shapemer_file):
    with open(shapemer_file, 'w') as g:
      for (sid1, shape), (sid2, cnt) in izip(filterSeq(shape_file, sid_file), filterSeq(count_file, sid_file)):
        assert(sid1 == sid2)
        count = int(cnt)
        shape = shape[self.lbc_len : -self.rbc_len]
        for k in range(len(shape)-shape_length+1):    # TODO: we need rep of either strand
          rep = getRepresent(shape[k:k+shape_length])
          rev = getRepresent(shape[k:k+shape_length][::-1])
          for j in range(count):
            print >>g, rep, sid1
            print >>g, rev, sid1

