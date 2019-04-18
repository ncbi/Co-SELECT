from dodo_common import *
from params_config import *
from doit.tools import run_once
from promiscuous_scripts.gen_maps_seq_shape import *

en_thresholds = ["1.20"]
discrete_levels_type = ["publish"]
cycles = [4]
shapes = ['MGW']
seq_lengths = [9, 10]
num_iter = 2

topdir = '/panfs/pan1/dnashape/seqmer2shapemer'

def task_gen_rand_seq():
  """Generate random sequences of length 9 or 10"""
  for seq_len in seq_lengths:
    seq_file = topdir + "/all_{}N_seq.txt".format(seq_len) 
    yield {
      'name'      : seq_file,
      'actions'   : [(gen_rand_seq, [seq_len, seq_file])],
      'uptodate'  : [run_once],
      'targets'   : [seq_file],
      'clean'     : True,
    }
  

def task_gen_shape2seq_map():
  """Generate mapping of sequence mers to shapemes"""
  for shape_type in shapes:
    for levels_type in discrete_levels_type:
      seq_len = 9 if shape_type in ['HelT', 'Roll'] else 10
      seq_file = topdir + "/all_{}N_seq.txt".format(seq_len) 
      map_file = topdir + "/map_seq2shape_{}.csv".format(shape_type)
      yield {
        'name'      : map_file,
        'actions'   : [(gen_shape2seq_map, [shape_type, levels_type, seq_file, map_file])],
        'file_dep'  : [seq_file],
        'targets'   : [map_file],
        'clean'     : True,
      }


def task_get_bg_summary():
  """Get information about enrichment of shapemers in motif-free pool (bg)
     of all TF experiments"""
  for shape_type in shapes:
    for levels_type in discrete_levels_type:
      for en_th in en_thresholds:
        outfile = "%s/%s/%s/bg_summary_%s.csv" % (top_results_dir, levels_type, en_th, shape_type)
        yield {
          'name'      : outfile,
          'actions'   : ["promiscuous_scripts/get_bg_summary.R -s %s -o %s" % (shape_type, outfile)],
          'file_dep'  : [],
          'targets'   : [outfile],
          'clean'     : True,
        }



def task_permulation_test():
  """Compute statistical significance of promiscuous shapemers using permutation test"""
  for shape_type in shapes:
    for levels_type in discrete_levels_type:
      for en_th in en_thresholds:
        for cycle in cycles:
          bgfile = "%s/%s/%s/bg_summary_%s.csv" % (top_results_dir, levels_type, en_th, shape_type)
          mapfile = topdir + "/map_seq2shape_{}.csv".format(shape_type)
          promisfile = '%s/%s/%s/promiscuous_%s.%s.%d.l1.r1.csv' % (top_results_dir,
              levels_type, en_th, 'same', shape_type, cycle)
          outfile = '%s/%s/%s/promiscuous_significance_%s.%s.%d.l1.r1.csv' % (
              top_results_dir, levels_type, en_th, 'same', shape_type, cycle)
          yield {
            'name'      : outfile,
            'actions'   : ["promiscuous_scripts/get_enr_random_seq.R "
                "-n %d -m %s -b %s -p %s -o %s" % (num_iter, mapfile, bgfile, promisfile, outfile)],
            'file_dep'  : [mapfile, bgfile, promisfile],
            'targets'   : [outfile],
            'clean'     : True,
          }

