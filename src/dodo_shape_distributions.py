from dodo_common import *

def task_flatten_dist_files():
  for shape_type in shapes:
    infile = "%s/%s.dist" % (top_shape_dist_dir, shape_type)
    outfile = "%s/%s_flatten.dist" % (top_shape_dist_dir, shape_type)
    yield {
      'name'      : outfile,
      'actions'   : ["shape_distributions_scripts/flatten_dist.R %s %s" % (infile, outfile)],
      'file_dep'  : [infile],
      'targets'   : [outfile],
      'clean'     : True,
    }

def task_unmix_gaussians():
  for shape_type in shapes:
    for num_comp in [2,3,4,5]:
      infile = "%s/%s_flatten.dist" % (top_shape_dist_dir, shape_type)
      outfile = "%s/unmixed_gaussians_%s_k%d.dist" % (top_shape_dist_dir, shape_type, num_comp)
      yield {
        'name'      : outfile,
        'actions'   : ["shape_distributions_scripts/unmix_gaussians.R %s %d %s" % (infile, num_comp, outfile)],
        'file_dep'  : [infile],
        'targets'   : [outfile],
        'clean'     : True,
      }

