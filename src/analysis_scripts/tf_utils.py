from subprocess import Popen, PIPE
from itertools import izip
from re import finditer
import os, re


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def rev_comp(s):
  n = len(s)
  t = [None]*n
  for i in range(n):
    z = s[n-i-1]
    if (z == 'A'): t[i] = 'T'
    elif (z == 'C'): t[i] = 'G'
    elif (z == 'G'): t[i] = 'C'
    elif (z == 'T'): t[i] = 'A'
    else: print "Unknown base", z; exit(0)
  return ''.join(t)


def hamming_single(s, t):
  return sum([s[i] != t[i] for i in range(len(s))])

def hamming_all(seq, pat):
  n = len(seq)
  k = len(pat)
  return min([hamming_single(seq[i:i+k], pat) for i in range(n-k+1)])


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

def find_all_occur(seq, pat):
  #return [m.start() for m in finditer(pat, seq)]
  t = []
  pos = seq.find(pat)
  while (pos>=0):
    t.append(pos)
    #pos = seq.find(pat, pos+len(pat))
    pos = seq.find(pat, pos+1)
  return t

def unzip_seq_filter_N_old(in_file, seq_file, count_file):
  p = Popen(["zcat", in_file], stdout=PIPE)
  counts = {}
  with  p.stdout as f:
    for l1, l2, l3, l4 in izip(f,f,f,f):
      s = l2.rstrip()
      if (s.find('N') < 0): # doesnot have N
        try:
          counts[s] += 1
        except:
          counts[s] = 1
  with open(seq_file, "w") as g, open(count_file, "w") as h:
    for s,c in sorted(counts.items(), key=lambda x: (x[1],x[0]), reverse=True):
        print >>g, s
        print >>h, c

def unzip_seq_filter_N(barcode, in_file, seq_file, count_file):
  res = re.match('([ACGT]+)\d+N([ACGT]+)', barcode)
  left = res.group(1)
  right = res.group(2)
  print left, right
  p = Popen(["zcat", in_file], stdout=PIPE)
  counts = {}
  with  p.stdout as f:
    for l1, l2, l3, l4 in izip(f,f,f,f):
      s = l2.rstrip()
      if (s.find('N') < 0): # doesnot have N
        try:
          counts[s] += 1
        except:
          counts[s] = 1
  with open(seq_file, "w") as g, open(count_file, "w") as h:
    for s,c in sorted(counts.items(), key=lambda x: (x[1],x[0]), reverse=True):
        print >>g, ''.join([left, s, right])
        print >>h, c

