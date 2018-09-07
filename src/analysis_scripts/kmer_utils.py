from itertools import groupby, izip, product

def fasta_iter(fasta_name, sep=''):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">" FOR PYTHON 2.x user next instead of __next__
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = sep.join(s.strip() for s in faiter.next()) #FOR PYTHON 2.x user next instead of __next__
        yield header, seq


def filterSeq(full, part):
  i=-1
  with open(full) as f, open(part) as g:
    for gl in g:
      keep = int(gl.rstrip().split()[0])
      #print keep
      for i, l in enumerate(f, i+1):
        #print "    ", i, l
        if (i==keep):
          yield keep, l.rstrip()
          break

def filterFasta(full, part, sep):
  i=-1
  with open(part) as g:
    fi = fasta_iter(full, sep)
    for gl in g:
      keep = int(gl.rstrip().split()[0])
      for i, (header, seq) in enumerate(fi, i+1):
        #print "i", i, "header", header
        assert(int(header) == i+1)
        if (i==keep):
          yield header, seq
          break
    
def shape_count_iter(shape, shape_file, count_file):
  sep = ','
  with open(count_file) as g:
    for (header, seq), cnt in izip(fasta_iter(shape_file, sep), g):
      yield shape.encodeList(seq.split(sep)), int(cnt)

 
def shape_iter(shape, shape_file):
  sep = ','
  for header, seq in fasta_iter(shape_file, sep):
    yield shape.encodeList(seq.split(sep))


def shape_continuous_count_iter(shape_file, count_file):
  sep = ','
  with open(count_file) as g:
    for (header, seq), cnt in izip(fasta_iter(shape_file, sep), g):
      yield header, seq, int(cnt)

def filter_shape_continuous_count_iter(shape_file, count_file, filter_file):
  sep = ','
  for (header, seq), (i, cnt) in izip(filterFasta(shape_file, filter_file, sep), filterSeq(count_file, filter_file)):
    assert(int(header) == i+1)
    yield header, seq, int(cnt)


def filter_shape_count_iter(shape, shape_file, count_file, filter_file):
  sep = ','
  for (header, seq), (i, cnt) in izip(filterFasta(shape_file, filter_file, sep), filterSeq(count_file, filter_file)):
    assert(int(header) == i+1)
    yield header, shape.encodeList(seq.split(sep)), int(cnt)


def countShapeMers(shape, shape_files, kmer_length):
  count = {''.join(x):0 for x in product(shape.getAlpha(), repeat=kmer_length)}
  total = 0
  for shape_file, cnt_file in shape_files:
    for seq, wt in shape_count_iter(shape, shape_file, cnt_file):
      for i in range(len(seq)-kmer_length+1):
        kmer = seq[i:i+kmer_length]
        try:
          count[kmer] += wt
        except:
          count[kmer] = wt
        total += wt
  return count, total


def countFilteredShapeMers(shape, shape_files, kmer_length):
  count = {''.join(x):0 for x in product(shape.getAlpha(), repeat=kmer_length)}
  total = 0
  for shape_file, cnt_file, filter_file in shape_files:
    for header, seq, wt in filter_shape_count_iter(shape, shape_file, cnt_file, filter_file):
      for i in range(len(seq)-kmer_length+1):
        kmer = seq[i:i+kmer_length]
        try:
          count[kmer] += wt
        except:
          count[kmer] = wt
        total += wt
  return count, total


def loadKmerCounts(count_file, kmer_length):
  count = {}
  total = 0
  with open(count_file) as f:
    for l in f:
      val = l.rstrip().split()
      kmer = val[0]
      cnt = int(val[1])
      count[kmer] = cnt
      total += cnt
  return count, total

def saveKmerCounts(count_file, count):
  with open(count_file, "w") as f:
    for kmer in sorted(count.keys()):
      f.write("%s\t%d\n" % (kmer, count[kmer]))

