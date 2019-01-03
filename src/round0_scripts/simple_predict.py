import pandas as pd
import os

#def getSequenceFile(tf, primer, cycle):
#  return "%s/%s/%d.non" % (tf, primer, cycle)
#
#def getContextFile(tf, primer, cycle, motif, lflank, rflank, ctx):
#  fname = "%s.%s" % (getSequenceFile(tf, primer, cycle), motif)
#  if ctx == 'fg':
#    fname += ".l%d.r%d" % (lflank, rflank)
#  fname += ".%s" % (ctx)
#  return fname
#
#def getContextedShapemerCountFile(tf, primer, cycle, motif, lflank, rflank, ctx, shape):
#  return getContextFile(tf, primer, cycle, motif, lflank, rflank, ctx) + '.' + shape  + '.mer.cnt'

def getCycle0ProbFile(probability_dir, ctx, shape, motif, dist, lflank, rflank, shape_levels_str):
  fname = "%s/prob.%s.%s.%s.%s" % (probability_dir, ctx, motif, shape, shape_levels_str)
  fname += ".d%s" % (dist)
  if ctx == 'fg':
    fname += ".l%d.r%d" % (lflank, rflank)
  fname += ".txt"
  return fname



def train_all(infiles, outfile):
  print "training", outfile
  count = {}
  total = 0
  for fname in infiles:
      if (not os.path.exists(fname)):
        print fname, "does not exists. skipping ..."
        continue
      with open(fname) as g:
        for l in g:
          vals = l.strip().split()
          kmer = vals[1]
          cnt = int(vals[0])
          try:
            count[kmer] += cnt
          except:
            count[kmer] = cnt
          total += cnt

  with open(outfile, 'w') as f:
    for kmer in count:
      f.write("%s\t%d\t%.20f\n" % (kmer, count[kmer], 1.0*count[kmer]/total))

if __name__ == "__main__":
  topdir = '/panfs/pan1/aptax/PRJEB14744-by-original-name'
  inventory_dir = '../inventory'
  zdf = pd.read_csv(inventory_dir + '/our_data_zero_cycle.csv')
  zdf = zdf[zdf['primer'].str.contains('14N')==False]
  primers = zdf['primer'].tolist()
  for shape in ['MGW', 'ProT', 'HelT', 'Roll']:
    for context in ['fg', 'bg']:
      for motif in ['TAAT']:
        for l,r in [(1,1)]:
          train_all(topdir, primers, context, shape, motif, l, r, 'tmp_haha.txt')
