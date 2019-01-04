from dodo_common import *
from doit.tools import run_once

en_thresholds = ["1.10", "1.20"]
discrete_levels_type = ["publish", "other1"]
levels_names = {'publish': 'main-text', 'other1': 'alternative'}

print en_thresholds

def write_tfs_per_family(tfs, outfile):
  tmp = tfs[['tf', 'primer', 'family']]
  tmp = tmp.groupby('family').count()
  tmp['num_tfs'] = tmp['tf']
  tmp[['num_tfs']].to_csv(outfile, index=True) # [[]] needed for saving the column names
  


def combine_data_frames(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (en_th, shape_type) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['en_th'] = en_th
    tmp['shape'] = shape_type
    df = df.append(tmp, ignore_index=True)
  print df
  df.to_csv(outfile, index=False)



def task_count_per_family():
  output = '%s/tfs_per_family.csv' % (top_results_dir) 
  yield {
    'name'      : output,
    'actions'   : [(write_tfs_per_family, [tfs, output])],
    'file_dep'  : [tf_info_file, tf_run_coselect_file],
    'targets'   : [output],
    'clean'     : True,
  }
  

def task_detect_promiscuous():
  """ Get summary information for enrichment """
  tfs_per_family_file = '%s/tfs_per_family.csv' % (top_results_dir) 
  for en_th in en_thresholds:
    for shape_type in shapes:
      for cycle in cycles:
        for levels_type in discrete_levels_type:
          for lflank, rflank in flank_configs:
            enriched = '%s/%s/%s/denriched_%s.%s.%d.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, 'same', shape_type, cycle, lflank, rflank)
            output = '%s/%s/%s/promiscuous_%s.%s.%d.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, 'same', shape_type, cycle, lflank, rflank)
            yield {
              'name'      : output,
              'actions'   : ["promiscuous_scripts/detect_promiscuous.R -s %s -t %s -i %s -o %s -n %s" % (shape_type, en_th, enriched, output, tfs_per_family_file)],
              'file_dep'  : [enriched],
              'targets'   : [output],
              'clean'     : True,
            }

def task_combine_promiscuous():
  """ Get summary information for enrichment """
  for levels_type in discrete_levels_type:
    for lflank, rflank in flank_configs:
      for cycle in cycles:
        infiles = ['%s/%s/%s/promiscuous_%s.%s.%d.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, 'same', shape_type, cycle, lflank, rflank) for en_th in en_thresholds for shape_type in shapes]
        infos = [(en_th, shape_type) for en_th in en_thresholds for shape_type in shapes]
        output = '%s/promiscuous_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
        yield {
          'name'      : output,
          'actions'   : [(combine_data_frames, [infiles, infos, output])],
          'file_dep'  : infiles,
          'targets'   : [output],
          'clean'     : True,
        }


def task_plot_promiscuous():
  for levels_type in discrete_levels_type:
    for lflank, rflank in flank_configs:
      for cycle in cycles:
        infile = '%s/promiscuous_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
        outfile = '%s/highly_promiscuous_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
        selected_pdf = '%s/fig_promiscuous_shapemers_%s_cycle%d_%s.pdf' % (top_results_dir, levels_type, cycle, fg_type)
        combined_pdf = '%s/fig_promiscuous_shapemers_compare_threshold_%s_cycle%d_%s.pdf' % (top_results_dir, levels_type, cycle, fg_type)
        yield {
          'name'      : selected_pdf,
          'actions'   : ["promiscuous_scripts/plot_promiscuous.R %s %s %s %s" % (infile, outfile, selected_pdf, combined_pdf)],
          'file_dep'  : [infile],
          'targets'   : [outfile, selected_pdf, combined_pdf],
          'clean'     : True,
        }



