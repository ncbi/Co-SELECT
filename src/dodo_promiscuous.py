from dodo_common import *
from doit.tools import run_once

topresdir = '/panfs/pan1/aptax/results'

en_thresholds = ["1.10", "1.20"]
discrete_levels_type = ["publish"]
shapes = ['MGW', 'ProT', 'HelT', 'Roll']
lflank, rflank = 1, 1
cycle = 4

print en_thresholds

def task_detect_promiscuous():
  """ Get summary information for enrichment """
  for en_th in en_thresholds:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        enriched = '%s/%s/%s/denriched_%s.%s.l%d.r%d.csv' % (topresdir, levels_type, en_th, 'same', shape_type, lflank, rflank)
        output = '%s/%s/%s/promiscuous_%s.%s.l%d.r%d.csv' % (topresdir, levels_type, en_th, 'same', shape_type, lflank, rflank)
        yield {
          'name'      : output,
          'actions'   : ["../scripts/detect_promiscuous.R -s %s -t %s" % (shape_type, en_th)],
          'file_dep'  : [enriched],
          'targets'   : [output],
          'clean'     : True,
        }

def task_get_bg_summary():
  """ Get summary information for enrichment """
  for shape_type in shapes:
      output = '%s/%s/%s/bg_summary_%s.csv' % (topresdir, 'publish', "1.20", shape_type)
      yield {
        'name'      : output,
        'actions'   : ["../scripts/get_bg_summary.R -s %s" % (shape_type)],
        #'file_dep'  : [enriched],
        'targets'   : [output],
        'uptodate'  : [run_once],
        'clean'     : True,
      }


def task_promiscuous_significance():
  """ Get summary information for enrichment """
  for shape_type in shapes:
      summary = '%s/%s/%s/bg_summary_%s.csv' % (topresdir, 'publish', "1.20", shape_type)
      output = '%s/%s/%s/promiscuous_significance_%s.csv' % (topresdir, 'publish', "1.20", shape_type)
      yield {
        'name'      : output,
        'actions'   : ["../scripts/get_enr_random_seq.R -s %s -n %d" % (shape_type, 1000)],
        'file_dep'  : [summary],
        'targets'   : [output],
        'clean'     : True,
      }

#def task_significant_promiscuous():
#  """ Get summary information for enrichment """
#  for en_th in en_thresholds:
#    for shape_type in shapes:
#      for levels_type in discrete_levels_type:
#        shinfo = shape_info[levels_type][shape_type]
#        shape_levels_str = shinfo.getLevelsStr()
#        promiscuous = '%s/%s/%s/promiscuous_%s.%d.%s.l%d.r%d.csv' % (topresdir, levels_type, en_th, 'same', cycle, shape_type, lflank, rflank)
#        yield {
#          'name'      : output,
#          'actions'   : ["../scripts/detect_promiscuous.R %s %s %d %d %d %s %s" % (shape_type, en_th, 'same', cycle, lflank, rflank, levels_type, output)],
#          'file_dep'  : [enriched],
#          'targets'   : [output],
#          'clean'     : True,
#        }
#
#


#def task_get_cross_summary():
#  """ Get summary information for enrichment """
#  for i, row in cross.iterrows():
#    tf1, primer1, motif1, family1, dist1 = row[['tf_x', 'primer_x', 'motif_x', 'family_x', 'distance_x']]
#    tf2, primer2, motif2, family2, dist2 = row[['tf_y', 'primer_y', 'motif_y', 'family_y', 'distance_y']]
#    for en_th in en_thresholds:
#      for lflank, rflank in flank_configs:
#        for shape_type in shapes:
#          for levels_type in discrete_levels_type:
#            shinfo = shape_info[levels_type][shape_type]
#            shape_levels_str = shinfo.getLevelsStr()
#            resdir = '%s/%s/%s/%s_%s' % (topresdir, levels_type, en_th, min(family1, family2), max(family1, family2))
#            for cycle in [4]:
#                info1 = TFInfo(tf1, primer1)
#                info2 = TFInfo(tf2, primer2)
#                input_files = [top_data_dir +'/'+ info1.getContextedShapemerCountFile(cycle, motif1, int(dist1), lflank, rflank, 'fg', shape_type, shape_levels_str)]
#                input_files += [top_data_dir +'/'+ info2.getContextedShapemerCountFile(cycle, motif2, int(dist2), lflank, rflank, 'bg', shape_type, shape_levels_str)]
#                file_fisher = '.'.join(['fisher', tf1, primer1, motif1, tf2, primer2, motif2, shape_type, str(cycle), str(lflank), str(rflank), 'csv'])
#                output_file = '/'.join([resdir, file_fisher])
#                prob_files = [getCycle0ProbFile('fg', shape_type, motif1, dist1, lflank, rflank, shape_levels_str)]
#                prob_files += [getCycle0ProbFile('bg', shape_type, motif2, dist2, lflank, rflank, shape_levels_str)]
#                yield {
#                  'name'      : output_file,
#                  'actions'   : ["../scripts/get_cross_summary.R %s %s %s %d %s %s %s %d %d %s %d %d %s %s %s" % (tf1, primer1, motif1, dist1, tf2, primer2, motif2, dist2, cycle, shape_type, lflank, rflank, en_th, resdir, shape_levels_str)],
#                  'file_dep'  : input_files + prob_files,
#                  'targets'   : [output_file],
#                  'clean'     : True,
#                }
#
#def task_same_fisher():
#  for en_th in en_thresholds:
#      for lflank, rflank in flank_configs:
#        for shape_type in shapes:
#          for levels_type in discrete_levels_type:
#            shinfo = shape_info[levels_type][shape_type]
#            shape_levels_str = shinfo.getLevelsStr()
#            for cycle in [4]:
#                inputs = ['/'.join([topresdir, levels_type, en_th, task.family, task.tf, task.primer, '.'.join(['fisher', task.tf, task.primer, shape_type, str(cycle), motif, str(lflank), str(rflank), 'csv'])]) for task in nonzero_task_infos for motif, dist in izip(task.motifs, task.distances)]
#                target = '%s/%s/%s/dfisher_%s.%d.%s.l%d.r%d.csv' % (topresdir, levels_type, en_th, 'same', cycle, shape_type, lflank, rflank)
#                yield {
#                  'name'      : target,
#                  'actions'   : ["../scripts/fisher.R %s %s %s %d %d %d %s" % (shape_type, en_th, 'same', cycle, lflank, rflank, levels_type)],
#                  'file_dep'  : inputs,
#                  'targets'   : [target],
#                  'clean'     : True,
#                }
#
#def task_cross_fisher():
#  combinations = []
#  for i, row in cross.iterrows():
#    tf1, primer1, motif1, family1, dist1 = row[['tf_x', 'primer_x', 'motif_x', 'family_x', 'distance_x']]
#    tf2, primer2, motif2, family2, dist2 = row[['tf_y', 'primer_y', 'motif_y', 'family_y', 'distance_y']]
#    resdir = '%s_%s' % (min(family1, family2), max(family1, family2))
#    combinations.append((resdir, tf1, primer1, motif1, family1, dist1, tf2, primer2, motif2, family2, dist2))
#  for en_th in en_thresholds:
#      for lflank, rflank in flank_configs:
#        for shape_type in shapes:
#          for levels_type in discrete_levels_type:
#            shinfo = shape_info[levels_type][shape_type]
#            shape_levels_str = shinfo.getLevelsStr()
#            for cycle in [4]:
#                inputs = ['/'.join([topresdir, levels_type, en_th, resdir, '.'.join(['fisher', tf1, primer1, motif1, tf2, primer2, motif2, shape_type, str(cycle), str(lflank), str(rflank), 'csv'])]) for resdir, tf1, primer1, motif1, family1, dist1, tf2, primer2, motif2, family2, dist2 in combinations]
#                target = '%s/%s/%s/dfisher_%s.%d.%s.l%d.r%d.csv' % (topresdir, levels_type, en_th, 'cross', cycle, shape_type, lflank, rflank)
#                yield {
#                  'name'      : target,
#                  'actions'   : ["../scripts/fisher.R %s %s %s %d %d %d %s" % (shape_type, en_th, 'cross', cycle, lflank, rflank, levels_type)],
#                  'file_dep'  : inputs,
#                  'targets'   : [target],
#                  'clean'     : True,
#                }
