from itertools import product
import subprocess as sp
from analysis_scripts.kmer_utils import *
from analysis_scripts.shape_info import *

levels_type = 'publish'
shape_types = ['MGW', 'HelT', 'ProT', 'Roll']
shape_types = ['HelT', 'Roll']
shape_types = ['MGW']

dnashape_exe = '/home/pals2/Work/aptax/DNAshapeR/src/dnashape'
topdir = '/panfs/pan1/dnashape/seqmer2shapemer'


for seq_len in [9,10]:
  seq_file = topdir + "/all_{}N_seq.txt".format(seq_len) 
  with open(seq_file, "w") as f:
    for x in product('ACGT', repeat=seq_len):
        seq = ''.join(x)
        f.write(">{0}\n{0}\n".format(seq))

for shape_type in shape_types:
  seq_len = 10
  if shape_type in ['HelT', 'Roll']:
    seq_len = 9
  seq_file = topdir + "/all_{}N_seq.txt".format(seq_len) 

  sp.call("%s %s %s" % (dnashape_exe, seq_file, shape_type), shell=True)

  shape = shape_info[levels_type][shape_type]
  shape_file = seq_file + '.' + shape_type
  print seq_file
  print shape_file

  sep=','
  shape_mers = {}
  map_file = topdir + "/map_seq2shape_{}.csv".format(shape_type)
  with open(map_file, "w") as f:
    f.write("%s,%s\n" % ("seq", "shapemer"))
    
    for header, seq in fasta_iter(shape_file, sep):
      shape_mer = shape.encodeList(seq.split(sep))
      assert(len(shape_mer) == 6)
      rev = shape_mer[::-1]
      if (rev < shape_mer): shape_mer = rev
      f.write("%s,%s\n" % (header, shape_mer))
      try:
        shape_mers[shape_mer].append(header)
      except:
        shape_mers[shape_mer] = [header]
    
    
    print "%s,%s,%s,%s" % ('SHAPE', 'SHAPEMER', 'NUM_SEQS', 'SEQS')
    for mer in sorted(shape_mers.keys()):
      print "%s,%s,%s,%s" % (shape_type, mer, len(shape_mers[mer]), ';'.join(shape_mers[mer]))

