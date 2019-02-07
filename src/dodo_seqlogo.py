from dodo_common import *
from PyPDF2 import PdfFileMerger

shapes = ['MGW']
discrete_levels_type = ['publish']
en_thresholds = ['1.20', '1.10']


task_infos = []

for i, row in tfs.iterrows():
  tf = row['tf']
  bc = row['primer']
  motif = row['motif']
  family = row['family']
  dist = row['distance']
  accessions = row['accessions']
  task_infos.append(TaskInfo(tf, bc, family, accessions, [motif], cycles, [dist]))


def merge_pdfs(infiles, outfile):
  merger = PdfFileMerger()
  for pdf in infiles:
    merger.append(open(pdf, 'rb'))
  with open(outfile, 'wb') as fout:
    merger.write(fout)


def merge_pdfs_from_list(inlist, outfile):
  merger = PdfFileMerger()
  infiles = pd.read_csv(inlist)['pdf']
  for pdf in infiles:
    merger.append(open(pdf, 'rb'))
  with open(outfile, 'wb') as fout:
    merger.write(fout)

def task_compute_seqlogo_pwms():
  """ Get summary information for enrichment """
  for en_th in en_thresholds:
    for task in task_infos:
      for lflank, rflank in flank_configs:
        for shape_type in shapes:
          for levels_type in discrete_levels_type:
            shinfo = shape_info[levels_type][shape_type]
            shape_levels_str = shinfo.getLevelsStr()
            for i, cycle in enumerate(task.cycles):
              for motif, dist in izip(task.motifs, task.distances):
                enrich_file = '/'.join([top_results_dir, levels_type, en_th, task.family, task.tf, task.primer, 
                                        '.'.join(['enriched', task.tf_info.tf, task.tf_info.primer, shape_type, str(4), motif, str(lflank), str(rflank), 'csv'])])
                promiscuous_file = '%s/highly_promiscuous_%s_cycle%d.l%d.r%d.csv' % (top_results_dir, levels_type, 4, lflank, rflank)
                input_files = [  "%s/%s" % (orig_data_dir, task.tf_info.getSequenceFile(cycle)),
                                 "%s/%s" % (orig_data_dir, task.tf_info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str)),
                                 "%s/%s.cnt" % (orig_data_dir, task.tf_info.getSequenceFile(cycle)),
                                 "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'fg')),
                                 "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'bg')),
                                 "%s/%s" % (top_data_dir, task.tf_info.getFgPartsFile(cycle, motif, dist, lflank, rflank)),
                                 enrich_file, promiscuous_file
                              ]
                outdir = '/'.join([top_seqlogo_dir, levels_type, en_th, task.family, task.tf, task.primer])
                file_logo = '.'.join(['seqlogo', task.tf, task.primer, shape_type, str(cycle), motif, str(lflank), str(rflank), 'csv'])
                output_file = '/'.join([outdir, file_logo])
                ensure_dir(output_file)
                yield {
                  'name'      : output_file,
                  'actions'   : ["seqlogo_scripts/get_seqlogo_pwm.py %s %s %s %s %d %s %s %s %s %s %s %s" % (task.tf, task.primer, task.family, motif, cycle, en_th, orig_data_dir, top_data_dir, top_results_dir, enrich_file, promiscuous_file, output_file)],
                  'file_dep'  : input_files,
                  'targets'   : [output_file],
                  'clean'     : True,
                }
                if (i == 0): # we need to do for cycle 0 only once for all task.cycles
                  input_files = [  "%s/%s" % (orig_data_dir, task.tf_info.getSequenceFile(0)),
                                   "%s/%s" % (orig_data_dir, task.tf_info.getDiscreteShapeFile(0, shape_type, shape_levels_str)),
                                   "%s/%s.cnt" % (orig_data_dir, task.tf_info.getSequenceFile(0)),
                                   "%s/%s" % (top_data_dir, task.tf_info.getContextFile(0, motif, dist, 'fg')),
                                   "%s/%s" % (top_data_dir, task.tf_info.getContextFile(0, motif, dist, 'bg')),
                                   "%s/%s" % (top_data_dir, task.tf_info.getFgPartsFile(0, motif, dist, lflank, rflank)),
                                   enrich_file, promiscuous_file
                                ]
                  outdir = '/'.join([top_seqlogo_dir, levels_type, en_th, task.family, task.tf, task.primer])
                  file_logo = '.'.join(['seqlogo', task.tf, task.primer, shape_type, str(0), motif, str(lflank), str(rflank), 'csv'])
                  output_file = '/'.join([outdir, file_logo])
                  ensure_dir(output_file)
                  yield {
                    'name'      : output_file,
                    'actions'   : ["seqlogo_scripts/get_seqlogo_pwm.py %s %s %s %s %d %s %s %s %s %s %s %s" % (task.tf, task.primer, task.family, motif, 0, en_th, orig_data_dir, top_data_dir, top_results_dir, enrich_file, promiscuous_file, output_file)],
                    'file_dep'  : input_files,
                    'targets'   : [output_file],
                    'clean'     : True,
                  }

def combine_csvs(configs, outfile):
  index_cols = list(set(configs.columns) - set(['infile']))
  df = pd.DataFrame()
  for indx, r in configs.iterrows():
    tmp = pd.read_csv(r['infile'])
    for col in index_cols:
      tmp[col] = r[col]
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)



def task_combine_seqlogo_pwms():
  """ Get summary information for enrichment """
  for en_th in en_thresholds:
    for task in task_infos:
      for lflank, rflank in flank_configs:
        for shape_type in shapes:
          for levels_type in discrete_levels_type:
            for motif, dist in izip(task.motifs, task.distances):
              outdir = '/'.join([top_seqlogo_dir, levels_type, en_th, task.family, task.tf, task.primer])
              df = pd.DataFrame({'cycle':task.cycles+[0]})
              df['infile'] = df['cycle'].apply(lambda x: '/'.join([outdir, '.'.join(['seqlogo', task.tf, task.primer, shape_type, str(x), motif, str(lflank), str(rflank), 'csv'])]))
              outfile = "%s/%s" % (outdir, '.'.join(['seqlogo', task.tf, task.primer, shape_type, "allcycles", motif, str(lflank), str(rflank), 'csv']))
              yield {
                'name'      : outfile,
                'actions'   : [(combine_csvs, [df, outfile])],
                'file_dep'  : df['infile'].tolist(),
                'targets'   : [outfile],
                'clean'     : True,
              }


def task_plot_seqlogo_pwms():
  """ Get summary information for enrichment """
  for en_th in en_thresholds:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          for task in task_infos:
            for motif, dist in izip(task.motifs, task.distances):
              outdir = '/'.join([top_seqlogo_dir, levels_type, en_th, task.family, task.tf, task.primer])
              infile = "%s/%s" % (outdir, '.'.join(['seqlogo', task.tf, task.primer, shape_type, "allcycles", motif, str(lflank), str(rflank), 'csv']))
              outfile = "%s/%s" % (outdir, '.'.join(['seqlogo', task.tf, task.primer, shape_type, "allcycles", motif, str(lflank), str(rflank), 'pdf']))
              yield {
                'name'      : outfile,
                'actions'   : ["seqlogo_scripts/plot_seqlogo.R -i %s -t %s.%s.%s.%s -o %s" % (infile, task.family, task.tf, task.primer, motif, outfile)],
                'file_dep'  : [infile],
                'targets'   : [outfile],
                'clean'     : True,
              }


def task_combine_seqlogos():
  for en_th in en_thresholds:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          infiles = []
          for task in task_infos:
            for motif, dist in izip(task.motifs, task.distances):
              outdir = '/'.join([top_seqlogo_dir, levels_type, en_th, task.family, task.tf, task.primer])
              infile = "%s/%s" % (outdir, '.'.join(['seqlogo', task.tf, task.primer, shape_type, "allcycles", motif, str(lflank), str(rflank), 'csv']))
              infiles.append("%s/%s" % (outdir, '.'.join(['seqlogo', task.tf, task.primer, shape_type, "allcycles", motif, str(lflank), str(rflank), 'pdf'])))
          outfile = '%s/fig_seqlogo_enriched_shapemers_%s_%s_th%s.pdf' % (top_results_dir, levels_type, fg_type, en_th)
          yield {
            'name'      : outfile,
            'actions'   : [(merge_pdfs, [infiles, outfile])],
            'targets'   : [outfile],
            'file_dep'  : infiles,
            'clean'     : True,
          }



#def task_combine_promiscuous_seqlogo_pwms():
#  """ Get summary information for enrichment """
#  for en_th in en_thresholds:
#    for lflank, rflank in flank_configs:
#      for shape_type in shapes:
#        for levels_type in discrete_levels_type:
#          df = pd.DataFrame()
#          for task in task_infos:
#            for motif, dist in izip(task.motifs, task.distances):
#              indir = '/'.join([top_seqlogo_dir, levels_type, en_th, task.family, task.tf, task.primer])
#              infile = "%s/%s" % (indir, '.'.join(['seqlogo', task.tf, task.primer, shape_type, "allcycles", motif, str(lflank), str(rflank), 'csv']))
#              df = df.append(pd.DataFrame({'infile': [infile], 'tf':[task.tf], 'family':[task.family], 'primer':[task.primer], 'motif':[motif]}), ignore_index=True)
#          print(df)
#          outfile = "%s/%s" % (top_seqlogo_dir, '.'.join(['seqlogo', shape_type, "allcycles", str(lflank), str(rflank), 'csv']))
#          yield {
#            'name'      : outfile,
#            'actions'   : [(combine_csvs, [df, outfile])],
#            'file_dep'  : df['infile'].tolist(),
#            'targets'   : [outfile],
#            'clean'     : True,
#          }
#
#
#def task_plot_promiscuous_seqlogo():
#  """ Get summary information for enrichment """
#  for en_th in en_thresholds:
#    for lflank, rflank in flank_configs:
#      for shape_type in shapes:
#        for levels_type in discrete_levels_type:
#          infile = "%s/%s" % (top_seqlogo_dir, '.'.join(['seqlogo', shape_type, "allcycles", str(lflank), str(rflank), 'csv']))
#          pdflist = '%s/seqlogo_list_promiscuous_shapemers_%s_%s.csv' % (top_results_dir, levels_type, fg_type)
#          yield {
#            'name'      : pdflist,
#            'actions'   : ["seqlogo_scripts/plot_seqlogo_promiscuous.R -i %s -o %s" % (infile, pdflist)],
#            'file_dep'  : [infile],
#            'targets'   : [pdflist],
#            'clean'     : True,
#          }
#          outfile = '%s/fig_seqlogo_promiscuous_shapemers_%s_%s.pdf' % (top_results_dir, levels_type, fg_type)
#          yield {
#            'name'      : outfile,
#            'actions'   : [(merge_pdfs_from_list, [pdflist, outfile])],
#            'file_dep'  : [pdflist],
#            'targets'   : [outfile],
#            'clean'     : True,
#          }


