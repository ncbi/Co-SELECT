from dodo_common import *
from results_scripts.more_results import *

level_names = {'publish': 'main-text', 'other1': 'alternative'}


def task_generate_shape_levels_info():
  shape_levels_file = '%s/shape_levels.csv' % (top_results_dir)
  infile = 'analysis_scripts/shape_info.py'
  yield {
    'name'      : shape_levels_file,
    'actions'   : ["%s %s" % (infile, shape_levels_file)],
    'file_dep'  : [infile],
    'targets'   : [shape_levels_file],
    'clean'     : True,
  }



def task_plot_qvalues():
  for lflank, rflank in flank_configs:
    for cycle in cycles:
      infile = '%s/dqvalue.%d.l%d.r%d.csv' % (top_results_dir, cycle, lflank, rflank)
      shape_levels_file = '%s/shape_levels.csv' % (top_results_dir)
      selected_pdf = '%s/fig_qvalue_selected_cycle%d.pdf' % (top_results_dir, cycle)
      separate_pdf = '%s/fig_qvalue_separate_family_cycle%d.pdf' % (top_results_dir, cycle)
      combined_pdf = '%s/fig_qvalue_combined_family_cycle%d.pdf' % (top_results_dir, cycle)
      yield {
        'name'      : selected_pdf,
        'actions'   : ["results_scripts/plot_qvalue.R %s %s %s %s %s" % (infile, shape_levels_file, selected_pdf, separate_pdf, combined_pdf)],
        'file_dep'  : [infile, shape_levels_file],
        'targets'   : [selected_pdf, separate_pdf, combined_pdf],
        'clean'     : True,
      }



def task_get_fdr_table():
  for lflank, rflank in flank_configs:
    for cycle in cycles:
      infile = '%s/dqvalue.%d.l%d.r%d.csv' % (top_results_dir, cycle, lflank, rflank)
      table_tex = '%s/table_significant_tfs_at_fdr_cycle%d.tex' % (top_results_dir,cycle)
      table_pdf = '%s/table_significant_tfs_at_fdr_temp_cycle%d.pdf' % (top_results_dir,cycle)
      yield {
        'name'      : table_tex,
        'actions'   : ["results_scripts/compute_fdr_table.py %s %s %s" % (infile, table_tex, table_pdf)],
        'file_dep'  : [infile],
        'targets'   : [table_tex, table_pdf],
        'clean'     : True,
      }


def combine_enriched(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (en_th, shape_type) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['en_th'] = en_th
    tmp['shape'] = shape_type
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def task_combine_all_enriched_shapemers():
  for lflank, rflank in flank_configs:
    for levels_type in discrete_levels_type:
      for cycle in cycles:
        infiles = ["%s/%s/%s/denriched_same.%s.%d.l%d.r%d.csv" % (top_results_dir, levels_type, en_th, shape_type, cycle, lflank, rflank) for en_th in en_thresholds for shape_type in shapes]
        infos = [(en_th, shape_type)  for en_th in en_thresholds for shape_type in shapes]
        outfile = '%s/%s/denriched_same.%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
        yield {
          'name'      : outfile,
          'actions'   : [(combine_enriched, [infiles, infos, outfile])],
          'file_dep'  : infiles,
          'targets'   : [outfile],
          'clean'     : True,
        }


def task_get_enriched_shapemers_in_excel():
  for lflank, rflank in flank_configs:
    for ltype in discrete_levels_type:
      for cycle in cycles:
        infile = '%s/%s/denriched_same.%d.l%d.r%d.csv' % (top_results_dir, ltype, cycle, lflank, rflank)
        outfile = '%s/enriched_shapemers_%s_cycle%d.l%d.r%d.xlsx' % (top_results_dir, ltype, cycle, lflank, rflank)
        csvfile = '%s/table_enriched_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, ltype, cycle, lflank, rflank)
        yield {
          'name'      : outfile,
          'actions'   : [(get_enriched_shapemers_in_excel, [infile, tf_motif_file, outfile, csvfile])],
          'file_dep'  : [tf_motif_file, infile],
          'targets'   : [outfile, csvfile],
          'clean'     : True,
        }


def combine_detailed_results(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (en_th, shape_type) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['en_th'] = en_th
    tmp['shape'] = shape_type
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)

def combine_detailed_results2(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, ltype in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['discretization'] = ltype
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def task_combine_detailed_results():
  for lflank, rflank in flank_configs:
    for cycle in cycles:
      to_combine = []
      for levels_type in discrete_levels_type:
        infiles = ["%s/%s/%s/dfisher_same.%s.%d.l%d.r%d.csv" % (top_results_dir, levels_type, en_th, shape_type, cycle, lflank, rflank) for en_th in en_thresholds for shape_type in shapes]
        infos = [(en_th, shape_type) for en_th in en_thresholds for shape_type in shapes]
        outfile = '%s/%s/dfisher_same.%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
        to_combine.append(outfile)
        yield {
          'name'      : outfile,
          'actions'   : [(combine_detailed_results, [infiles, infos, outfile])],
          'file_dep'  : infiles,
          'targets'   : [outfile],
          'clean'     : True,
        }
      resfile = '%s/dfisher_same.%d.l%d.r%d.csv' % (top_results_dir, cycle, lflank, rflank)
      yield {
        'name'      : resfile,
        'actions'   : [(combine_detailed_results2, [to_combine, discrete_levels_type, resfile])],
        'file_dep'  : to_combine,
        'targets'   : [resfile],
        'clean'     : True,
      }


def task_get_detailed_results_in_excel():
  for lflank, rflank in flank_configs:
    for cycle in cycles:
      infile = '%s/dfisher_same.%d.l%d.r%d.csv' % (top_results_dir, cycle, lflank, rflank)
      outfile = '%s/detailed_results_cycle%d.l%d.r%d.xlsx' % (top_results_dir, cycle, lflank, rflank)
      yield {
        'name'      : outfile,
        'actions'   : [(get_detailed_results_in_excel, [infile, tf_motif_file, outfile, level_names])],
        'file_dep'  : [tf_motif_file, infile],
        'targets'   : [outfile],
        'clean'     : True,
      }


def task_cluster_by_enriched_shapemers():
  for levels_type in discrete_levels_type:
    for lflank, rflank in flank_configs:
      for cycle in cycles:
        promiscuous = '%s/promiscuous_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
        infile = '%s/table_enriched_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
        heatmap = '%s/fig_heatmap_shapemers_%s.cycle%d.l%d.r%d.pdf' % (top_results_dir, levels_type, cycle, lflank, rflank)
        pca_csv = '%s/pca_shapemers_%s.cycle%d.l%d.r%d.csv' % (top_results_dir, levels_type, cycle, lflank, rflank)
        pca_pdf = '%s/fig_pca_shapemers_%s.cycle%d.l%d.r%d.pdf' % (top_results_dir, levels_type, cycle, lflank, rflank)
        yield {
          'name'      : pca_csv,
          'actions'   : ["results_scripts/cluster_by_shapemers.py %s %s %s %s" % (infile, promiscuous, heatmap, pca_csv)],
          'file_dep'  : [infile, promiscuous],
          'targets'   : [heatmap, pca_csv],
          'clean'     : True,
        }
        yield {
          'name'      : pca_pdf,
          'actions'   : ["results_scripts/plot_pca.R %s %s" % (pca_csv, pca_pdf)],
          'file_dep'  : [pca_csv],
          'targets'   : [pca_pdf],
          'clean'     : True,
        }


