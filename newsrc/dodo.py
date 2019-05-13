from judi import add_param, show_param_db, Task, File
import pandas as pd

import sys, os
#sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.dirname("../src/"))

from analysis_scripts.tf_utils import *
from tf_process import *

dnashape_exe              = '../DNAshapeR/src/dnashape'
ena_project_file          = '../src/PRJEB14744.txt'
nonzero_accession_file    = '../src/PRJEB14744_nonzero_cycle.csv'
zero_accession_file       = '../src/PRJEB14744_zero_cycle.csv'
tf_info_file              = '../src/tf_inventory_jolma_ronshamir.csv'
tf_motif_file             = '../src/tf_coremotif.csv'
tf_run_coselect_file      = '../src/tf_run_coselect.csv'

orig_data_dir             = '../data'
seqmer_data_dir           = '../seqmerdata'
download_dir              = '../downloads'
top_shape_dist_dir        = '../shapedists'
store_dir                 = '../store'

contexts                  = ['fg', 'bg']
fg_types                  = ['d0', 'd1all', 'd1enriched']
fg_type                   = 'd1all'
fg_type                   = 'd1enriched'
fg_type                   = 'd0'

cycles                    = [1, 2, 3, 4] if fg_type == 'd0' else [4]
cycles                    = [3, 4] if fg_type == 'd0' else [4]

cycles                    = [4]

shapes                    = ["MGW", "ProT", "HelT", "Roll"]
shapes                    = ["MGW"]

flank_configs             = [(1,1)]

discrete_levels_type      = ['publish', 'other1', 'other2']

en_thresholds             = ['1.05', '1.10', '1.20', '1.50']

discrete_levels_type      = ['publish', 'other1']

seqmer_len                = 10


tf_info   = pd.read_csv(tf_info_file)
tf_motif  = pd.read_csv(tf_motif_file)
tf_run    = pd.read_csv(tf_run_coselect_file)

tfs = tf_info.merge(tf_motif).merge(tf_run)
tfs = tfs[['tf', 'primer', 'family', 'motif']]

#to remove
#tfs = tfs.head(2)


tfs_nz = tfs.assign(key=1).merge(pd.DataFrame({'cycle': cycles}).assign(key=1)).drop('key', 1)

accession = pd.read_csv(nonzero_accession_file).drop(columns='batch')
tfs_nz = tfs_nz.merge(accession)
print(tfs_nz)


zero_acc  = pd.read_csv(zero_accession_file)
barcodes = tfs[['primer', 'motif']].drop_duplicates()
tfs_z  = zero_acc.merge(barcodes)
tfs_z['cycle']  = 0
tfs_z['tf']     = 'ZeroCycle'
tfs_z['family'] = 'NoFamily'
print(zero_acc)

params = tfs_nz.append(tfs_z, sort=False)
params = params.rename(columns={'primer':'barcode'})

add_param(params)
show_param_db()


class Preprocess(Task):
  """ Unzip fastq files, keep only sequence info of those containing only ACGT """
  mask = ['family', 'motif']
  inputs = {'fastq': File('original', path=lambda x: "{}/{}.fastq.gz".format(download_dir, x['accession']), mask=mask)}
  targets = {'seq': File('NoN.fastq', root=store_dir, mask=mask),
             'cnt': File('NoN.fastq.cnt', root=store_dir, mask=mask)}
  actions = [(unzip_seq_filter_N, ['#barcode', '$fastq', '$seq', '$cnt'])]


class Partition(Task):
  """ Partition aptamer sequences into motif-containing (fg) and motif-free (bg)
  based on distance from MOTIF """
  mask = ['family']
  inputs = {'seq': Preprocess.targets['seq'],
            #'nbr': File('nbr.txt', root=store_dir)
            }
  targets = {'fg': File('fg.txt', mask=mask, root=store_dir),
             'bg': File('bg.txt', mask=mask, root=store_dir)}
  actions = [(partition_aptamers, [fg_type, '$seq', '#motif', 'hahah', '$fg', '$bg', '#barcode'])]


class GetBgSeqmers(Task):
  """ Get seqmers in the motif-free pool"""
  mask = ['family']
  inputs = {'bg': Partition.targets['bg'],
            'seq': Preprocess.targets['seq'],
            'cnt': Preprocess.targets['cnt'],
            }
  targets = {'out': File('seqmer.txt', mask=mask, root=store_dir)}
  actions = [(gen_bg_seqmers, ['$seq', '$cnt', '$bg', '$out', seqmer_len, '#barcode'])]


class CountBgSeqmers(Task):
  """ Count seqmers in the motif-free pool"""
  mask = ['family']
  inputs = {'inp': GetBgSeqmers.targets['out']}
  targets = {'out': File('seqmer.txt.cnt', mask=mask, root=store_dir)}
  actions = [('cat {} | sort | awk \'BEGIN {{OFS="\t"}} ($1 == last || last == "") {{sum += $2}} ($1 != last && last != "") {{print last, sum; sum = $2}} {{last = $1}} END {{print last, sum}}\' > {}', ['$inp', '$out'])]

