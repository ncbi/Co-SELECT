import pandas as pd
import re
from analysis_scripts.tf_utils import rev_comp

d1 = pd.read_csv('tf_run_coselect_131.csv')
d2 = pd.read_csv('tf_coremotif.csv')

d = d1.merge(d2)

print d

def identify_left_barcode(x):
  res = re.match('([ACGT]+)\d+N([ACGT]+)', x['primer'])
  return pd.Series({'left': res.group(1)})

def identify_right_barcode(x):
  res = re.match('([ACGT]+)\d+N([ACGT]+)', x['primer'])
  return pd.Series({'right': res.group(2)})



def hamming_pair(s, t, r):
  mat = sum([s[i:i+2] == t[i:i+2] for i in range(len(s)-1)])
  has_random = sum([sum(r[i:i+2]) if s[i:i+2] == t[i:i+2] else 0 for i in range(len(s)-1)])
  #print [[s[i:i+2], t[i:i+2]] for i in range(len(s)-1)]
  #print has_random
  return mat if has_random else 0

def hamming_pair_match(seq, pat, lbc_len, rbc_len):
  n = len(seq)
  k = len(pat)
  rand_part = [0]*lbc_len + [1]*(n-lbc_len-rbc_len) + [0]*rbc_len
  return max([hamming_pair(seq[i:i+k], pat, rand_part[i:i+k]) for i in range(n-k+1)])


def check_valid(barcode, motif, how):
  m = len(motif)
  b = len(barcode)
  #print b-m+1
  for i in range(max(0, b-m+1, b-2), b-1):
    #print "i=", i, "m-b+i-1", m-b+i-1
    assert(m-b+i-1 < m-1)
    for j in range(m-b+i):
      #print "b=", b, "i=", i, "m=", m, "j=", j , j+b-i-1
      assert(j+b-i-1 < m)
      #print i, j, barcode[i:i+2], motif[j:j+2]
      if barcode[i:i+2] == motif[j:j+2]:
        if barcode[b-1] == motif[j+b-i-1]:
          if j+b-i <= m-1:
            print barcode, motif, i, j, barcode[i:i+2], how, "!"
            return pd.Series({how+'m': barcode[i:i+2], how: 1})
        else:          
          if j+b-i <= m-2:
            print barcode, motif, i, j, barcode[i:i+2], how, "@"
            #print "j=", j, "m-4=", m-4
            return pd.Series({how+'m': barcode[i:i+2], how: 1})
  return pd.Series({how+'m': '', how: 0})


def check_valid_right(barcode, motif, how):
  m = len(motif)
  b = len(barcode)
  for i in range(min(1,m-2,b-1)):
    for j in range(i+1,m-1):
      #print how, barcode, i, barcode[i:i+2], motif, j, motif[j:j+2]
      if barcode[i:i+2] == motif[j:j+2]:
        assert(j-i >= 0)
        assert(j-i < m)
        if barcode[0] == motif[j-i]:
          if j-i >= 1:
            #print barcode, motif, i, j, barcode[i:i+2], how, "*"
            return pd.Series({how+'m': barcode[i:i+2], how: 1})
        else:          
          if j-i >= 2:
            #print barcode, motif, i, j, barcode[i:i+2], how, "+"
            return pd.Series({how+'m': barcode[i:i+2], how: 1})
  return pd.Series({how+'m': '', how: 0})


def check_left_motif(x):
  return check_valid(x['left'], x['motif'], 'lb')


def check_motif_right(x):
  return check_valid_right(x['right'], x['motif'], 'rb')


def check_motif_revcomp_right(x):
  return check_valid(rev_comp(x['right']), x['motif'], 'rcr')

def check_revcomp_left_motif(x):
  return check_valid_right(rev_comp(x['left']), x['motif'], 'rcl')


d = pd.concat([d, d.apply(identify_left_barcode, axis='columns')], axis='columns')
d = pd.concat([d, d.apply(check_left_motif, axis='columns')], axis='columns')


d = pd.concat([d, d.apply(identify_right_barcode, axis='columns')], axis='columns')
d = pd.concat([d, d.apply(check_motif_right, axis='columns')], axis='columns')

d['rrev'] = d['right'].apply(rev_comp)
d = pd.concat([d, d.apply(check_motif_revcomp_right, axis='columns')], axis='columns')

d['lrev'] = d['left'].apply(rev_comp)
d = pd.concat([d, d.apply(check_revcomp_left_motif, axis='columns')], axis='columns')

d['remove'] = d['lb']
d['remove'] = d['remove'] | d['rb']
d['remove'] = d['remove'] | d['rcr']
d['remove'] = d['remove'] | d['rcl']

print d
d.to_csv('trouble.csv')

d[d.remove == 0][['tf','primer']].to_csv("tf_run_coselect_x.csv", index=False)
