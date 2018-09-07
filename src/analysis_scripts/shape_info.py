class Alphabet:
  a = [
      "SH",
      "SMH",
      "SMHX",
  ]
  aold = [
      "LH",
      "LMH",
      "LMHS",
  ]
  def __init__(self, cutoffs = [3.967, 5.083]):
    num_levels = len(cutoffs) + 1
    if num_levels not in range(len(Alphabet.a[0]),len(Alphabet.a[-1])+1): num_levels = 2
    self.alpha = Alphabet.a[num_levels-len(Alphabet.a[0])]
    self.levels = cutoffs + [float('inf')]

  def __call__(self):
    return self.alpha

  def encode(self,x):
    for i, c in enumerate(self.levels):
      if (x < c): return self.alpha[i]

class Shape:
  def __init__(self, cutoffs, skip):
    self.cutoffs = cutoffs
    self.skip = skip
    self.alphabet = Alphabet(self.cutoffs)
    levels_str = ['']*(2*len(self.cutoffs)+1)
    levels_str[::2] = list(self.alphabet())
    levels_str[1::2] = ["%.3f" % x for x in self.cutoffs]
    self.levels_str = ''.join(levels_str)

  def getLevelsStr(self):
    return self.levels_str

  def getLevels(self):
    levels_str = ['']*(2*len(self.cutoffs)+1)
    levels_str[::2] = list(self.alphabet())
    levels_str[1::2] = ["%g" % x for x in self.cutoffs]
    return '|'.join(levels_str)
    

  def encodeList(self, l):
    return ''.join([self.alphabet.encode(float(x)) for x in l[self.skip : -self.skip]])

#'MGW' : Shape([3.9667, 5.083], 2),

shape_info = {
  'publish' : {
    'Roll' : Shape([-3, 1.2, 5], 1),
    'MGW'  : Shape([4.25, 5.4, 5.7], 2),
    'ProT' : Shape([-8.0, -5.5], 2),
    'HelT' : Shape([31.9, 34, 36], 1),
  },
  'other1' : {
    'Roll' : Shape([1.2], 1),
    'MGW'  : Shape([4.0, 4.6, 5.6], 2),
    'ProT' : Shape([-10.5, -5.0], 2),
    'HelT' : Shape([32, 34.5, 35.5], 1),
  }
}

if __name__ == "__main__":
  with open('shape_levels.csv', 'w') as f:
    f.write("{},{},{},{}\n".format('shape','levels_type','levels_str','levels'))
    for levels_type in ['publish', 'other1']:
      for shape in ['MGW', 'HelT', 'ProT', 'Roll']:
        print "{},{},{},{}".format(shape,levels_type,shape_info[levels_type][shape].getLevelsStr(),shape_info[levels_type][shape].getLevels())
        f.write("{},{},{},{}\n".format(shape,levels_type,shape_info[levels_type][shape].getLevelsStr(),shape_info[levels_type][shape].getLevels()))
