from itertools import groupby
import subprocess as sp
import pandas as pd

def filterSeq(full, part):
  i=-1
  with open(full) as f, open(part) as g:
    for gl in g:
      keep = int(gl.rstrip().split()[0])
      #print keep
      for i, l in enumerate(f, i+1):
        #print "    ", i, l.rstrip()
        if (i==keep):
          yield keep, l.rstrip()
          break

def getCoverage(seq_id_file, cnt_file, shapemer_seqid_file, shapemer_coverage_file):
  print seq_id_file
  print cnt_file
  num_seq = int(sp.check_output("cat %s | wc -l" % seq_id_file, shell=True))
  print num_seq
  count = {}
  for i, cnt in filterSeq(cnt_file, seq_id_file):
    count[i] = int(cnt)
  full_count = sum([count[x] for x in count.keys()])

  print shapemer_seqid_file
  print shapemer_coverage_file
  #print sorted(count.keys())
  #print count[1435383]

  with sp.Popen("cat %s | sort | uniq -c" % (shapemer_seqid_file), stdout=sp.PIPE, shell=True).stdout as f, \
       open(shapemer_coverage_file, 'w') as g, \
       open(shapemer_coverage_file+'.long', 'w') as h:
    for key, grp in groupby(f, lambda l: l.rstrip().split()[1]):
      tmp = [[int(l.rstrip().split()[2]), int(l.rstrip().split()[0])] for l in grp]
      #print key, tmp
      seqs = [t[0] for t in tmp]
      occur = [t[1] for t in tmp]
      seq_sum = sum([count[x] for x in seqs])
      occur_sum = sum([occur[i]*count[seqs[i]] for i in range(len(seqs))])
      print >>g, "%s,%d,%0.20f,%d,%0.20f,%d,%0.20f" % (key, len(seqs), 1.0*len(seqs)/num_seq, seq_sum, 1.0*seq_sum/full_count, occur_sum, 1.0*occur_sum/seq_sum)
      #print >>h, "%s,%s" % (key, ';'.join([str(x) for x in sorted(seqs)]))


def processTf(topdir, tf, bc):
  print(tf, bc)
  #for c in [0,1, 2, 3, 4]:
  #for c in [0,4]:
  for c in [0,4]:
    for ctx in ['fg', 'nmd2']:
      for shape in ['MGW']: #['ProT']: # , MGW
      #for shape in ['ProT', 'HelT']: #['ProT']: # , MGW
        tfdir = '{0}/{1}/{2}'.format(topdir, tf, bc)
        cntfile = '{1}/{0}.non.cnt'.format(c, tfdir)
        sidfile = '{2}/{0}.non.{1}'.format(c, ctx, tfdir)
        merfile = '{2}/{0}.non.{1}.{3}.mer'.format(c, ctx, tfdir, shape)
        covfile = '{2}/{0}.non.{1}.{3}.mer.cov'.format(c, ctx, tfdir, shape)
        getCoverage(sidfile, cntfile, merfile, covfile)

if __name__ == "__main__":
  topdir = '/panfs/pan1/aptax/NewTmpExcludeSingleton'
  topdir = '/panfs/pan1/aptax/NewTmp'
  topdir = '/panfs/pan1/aptax/PRJEB14744-by-original-name'

  #tf_info = pd.read_csv("tf_list.csv")
  #print(tf_info.head())
  #for i in range(10):
  #  processTf(topdir, tf_info.iloc[i]['TF'], tf_info.iloc[i]['BARCODE'])
  
  #processTf('GSC', 'TGGGAC20NGA')
  #processTf('ALX3', 'TGTAAA20NAAG')
  #processTf('DLX3', 'TATGTT20NCG')
  #processTf('SHOX2', 'TACGTC20NTGC')
  processTf(topdir, 'BSX', 'TATGAA20NCG')
