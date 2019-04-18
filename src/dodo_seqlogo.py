from dodo_common import *

task_infos = []

for i, row in tfs.iterrows():
  tf = row['tf']
  bc = row['primer']
  motif = row['motif']
  family = row['family']
  dist = row['distance']
  accessions = row['accessions']
  task_infos.append(TaskInfo(tf, bc, family, accessions, [motif], cycles, [dist]))


zfs = tfs.drop_duplicates(['primer','motif'])
zero_task_infos = []

for i, row in zfs.iterrows():
  tf = 'ZeroCycle'
  bc = row['primer']
  motif = row['motif']
  family = 'NoFamily'
  dist = row['distance']
  accessions = row['accessions']
  zero_task_infos.append(TaskInfo(tf, bc, family, accessions, [motif], [0], [dist]))


shapes = ['MGW']
discrete_levels_type = ['publish']
en_thresholds = ['1.20']
cycles = [1,2,3,4]


def task_get_fg_shapemers_seqmers():
  """ Generate shapemers from motif-containing oligos in selection rounds"""
  for task in task_infos:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in cycles:
            seq_file = "%s/%s" % (orig_data_dir, task.tf_info.getSequenceFile(cycle))
            shape_file = "%s/%s" % (orig_data_dir,
                task.tf_info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
            count_file = "%s.cnt" % (seq_file)
            for motif, dist in izip(task.motifs, task.distances):
              context_file = "%s/%s" % (top_data_dir,
                  task.tf_info.getContextFile(cycle, motif, dist, 'fg'))
              parts_file = "%s/%s" % (top_data_dir,
                  task.tf_info.getFgPartsFile(cycle, motif, dist, lflank, rflank))
              shapemer_file = "%s/%s.new" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      lflank, rflank, 'fg',shape_type, shape_levels_str))
              yield {
                'name'      : shapemer_file,
                'actions'   : [(task.tf_info.gen_fg_shapemers_seqmers, [shinfo,
                    task.shape_length, seq_file, shape_file, count_file,
                    context_file, motif, lflank, rflank, parts_file, shapemer_file])],
                'file_dep'  : [shape_file, count_file, context_file],
                'targets'   : [shapemer_file],
                'clean'     : True,
              }


def task_get_fg_shapemers_seqmers_cycle_zero():
  """ Generate shapemers from motif-containing oligos in initial pool"""
  for task in zero_task_infos:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in ['publish']:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in [0]:
            seq_file = "%s/%s" % (orig_data_dir,
                task.tf_info.getSequenceFile(cycle))
            shape_file = "%s/%s" % (orig_data_dir,
                task.tf_info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
            count_file = "%s.cnt" % (seq_file)
            for motif, dist in izip(task.motifs, task.distances):
              context_file = "%s/%s" % (top_data_dir,
                  task.tf_info.getContextFile(cycle, motif, dist, 'fg'))
              parts_file = "%s/%s" % (top_data_dir,
                  task.tf_info.getFgPartsFile(cycle, motif, dist, lflank, rflank))
              shapemer_file = "%s/%s.new" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      lflank, rflank, 'fg',shape_type, shape_levels_str))
              yield {
                'name'      : shapemer_file,
                'actions'   : [(task.tf_info.gen_fg_shapemers_seqmers, [shinfo,
                    task.shape_length, seq_file, shape_file, count_file,
                    context_file, motif, lflank, rflank, parts_file, shapemer_file])],
                'file_dep'  : [shape_file, count_file, context_file],
                'targets'   : [shapemer_file],
                'clean'     : True,
              }


def task_get_fg_shapemers_logo():
  """ Generate seqlogo of shapemers from motif-containing oligos in selection rounds"""
  for task in task_infos:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in cycles:
            for motif, dist in izip(task.motifs, task.distances):
              shapemer_file = "%s/%s.new" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      lflank, rflank, 'fg',shape_type, shape_levels_str))
              seqlogo_file = "%s/%s.pwm" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      lflank, rflank, 'fg',shape_type, shape_levels_str))
              yield {
                'name'      : seqlogo_file,
                'actions'   : ["analysis_scripts/get_logo_pwms.py %s %s" % (
                    shapemer_file, seqlogo_file)],
                'file_dep'  : [shapemer_file],
                'targets'   : [seqlogo_file],
                'clean'     : True,
              }

def task_get_fg_shapemers_logo_zero():
  """ Generate seqlogo of shapemers from motif-containing oligos in initial pool"""
  for task in zero_task_infos:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in [0]:
            for motif, dist in izip(task.motifs, task.distances):
              shapemer_file = "%s/%s.new" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      lflank, rflank, 'fg',shape_type, shape_levels_str))
              seqlogo_file = "%s/%s.pwm" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      lflank, rflank, 'fg',shape_type, shape_levels_str))
              yield {
                'name'      : seqlogo_file,
                'actions'   : ["analysis_scripts/get_logo_pwms.py %s %s" % (
                    shapemer_file, seqlogo_file)],
                'file_dep'  : [shapemer_file],
                'targets'   : [seqlogo_file],
                'clean'     : True,
              }



def task_get_bg_shapemers_seqmers():
  """ Generate shapemers from motif-free oligos in selection rounds"""
  for task in task_infos:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        for cycle in cycles:
          seq_file = "%s/%s" % (orig_data_dir, task.tf_info.getSequenceFile(cycle))
          shape_file = "%s/%s" % (orig_data_dir,
              task.tf_info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
          count_file = "%s.cnt" % (seq_file)
          for motif, dist in izip(task.motifs, task.distances):
              context_file = "%s/%s" % (top_data_dir,
                  task.tf_info.getContextFile(cycle, motif, dist, 'bg'))
              shapemer_file = "%s/%s.new" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      0, 0, 'bg',shape_type, shape_levels_str))
              yield {
                'name'      : shapemer_file,
                'actions'   : [(task.tf_info.gen_bg_shapemers_with_seqmers, [shinfo,
                    task.shape_length, seq_file, shape_file, count_file,
                    context_file, shapemer_file])],
                'file_dep'  : [shape_file, count_file, context_file],
                'targets'   : [shapemer_file],
                'clean'     : True,
              }


def task_get_bg_shapemers_seqmers_zero():
  """ Generate shapemers from motif-free oligos in initial pool"""
  for task in zero_task_infos:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        for cycle in [0]:
          seq_file = "%s/%s" % (orig_data_dir, task.tf_info.getSequenceFile(cycle))
          shape_file = "%s/%s" % (orig_data_dir,
              task.tf_info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
          count_file = "%s.cnt" % (seq_file)
          for motif, dist in izip(task.motifs, task.distances):
              context_file = "%s/%s" % (top_data_dir,
                  task.tf_info.getContextFile(cycle, motif, dist, 'bg'))
              shapemer_file = "%s/%s.new" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      0, 0, 'bg',shape_type, shape_levels_str))
              yield {
                'name'      : shapemer_file,
                'actions'   : [(task.tf_info.gen_bg_shapemers_with_seqmers, [shinfo,
                    task.shape_length, seq_file, shape_file, count_file,
                    context_file, shapemer_file])],
                'file_dep'  : [shape_file, count_file, context_file],
                'targets'   : [shapemer_file],
                'clean'     : True,
              }

def task_get_bg_shapemers_logo():
  """ Generate seqlogo of shapemers from motif-free oligos in selection rounds"""
  for task in task_infos:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        for cycle in cycles:
          for motif, dist in izip(task.motifs, task.distances):
              shapemer_file = "%s/%s.new" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      0, 0, 'bg',shape_type, shape_levels_str))
              seqlogo_file = "%s/%s.pwm" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      0, 0, 'bg',shape_type, shape_levels_str))
              yield {
                'name'      : seqlogo_file,
                'actions'   : ["analysis_scripts/get_logo_pwms.py %s %s" %(
                    shapemer_file, seqlogo_file)],
                'file_dep'  : [shapemer_file],
                'targets'   : [seqlogo_file],
                'clean'     : True,
              }

def task_get_bg_shapemers_logo_zero():
  """ Generate seqlogo of shapemers from motif-free oligos in initial pool"""
  for task in zero_task_infos:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        for cycle in [0]:
          for motif, dist in izip(task.motifs, task.distances):
              shapemer_file = "%s/%s.new" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      0, 0, 'bg',shape_type, shape_levels_str))
              seqlogo_file = "%s/%s.pwm" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      0, 0, 'bg',shape_type, shape_levels_str))
              yield {
                'name'      : seqlogo_file,
                'actions'   : ["analysis_scripts/get_logo_pwms.py %s %s" % (
                    shapemer_file, seqlogo_file)],
                'file_dep'  : [shapemer_file],
                'targets'   : [seqlogo_file],
                'clean'     : True,
              }


def task_get_combined_shapemers_logo():
  """ Merge seqlogo pwms for all cycles in both contexts (motif-free and
  motif-containing) in a TF experiment """
  for task in task_infos:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        df = pd.DataFrame(columns=['cycle', 'motif', 'context', 'infile'])
        for cycle in [0] + cycles:
          for motif, dist in izip(task.motifs, task.distances):
              seqlogo_file = "%s/%s.pwm" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      0, 0, 'bg',shape_type, shape_levels_str))
              df = df.append({'cycle':cycle, 'motif':motif, 'context':'bg',
                  'infile':seqlogo_file}, ignore_index=True)
        for lflank, rflank in flank_configs:
          for cycle in [0] + cycles:
            for motif, dist in izip(task.motifs, task.distances):
              seqlogo_file = "%s/%s.pwm" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(cycle, motif, dist,
                      lflank, rflank, 'fg',shape_type, shape_levels_str))
              df = df.append({'cycle':cycle, 'motif':motif, 'context':'fg',
                  'infile':seqlogo_file}, ignore_index=True)
        outfile = "%s/%s.pwm" % (top_data_dir,
            task.tf_info.getContextedShapemerFile(1111, 'all', 0, 0, 0,
                'combined',shape_type, shape_levels_str))
        yield {
          'name'      : outfile,
          'actions'   : [(combine_csvs, [df, outfile])],
          'file_dep'  : df['infile'].tolist(),
          'targets'   : [outfile],
          'clean'     : True,
        }


def task_get_combined_shapemers_logo_all():
  """ Merge seqlogo pwms for all cycles in both contexts (motif-free and
  motif-containing) of all TF experiments """
  for lflank, rflank in flank_configs:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        df = pd.DataFrame(columns=['family', 'tf', 'primer', 'infile'])
        for task in task_infos:
          infile = "%s/%s.pwm" % (top_data_dir,
              task.tf_info.getContextedShapemerFile(1111, 'all', 0, 0, 0,
                  'combined',shape_type, shape_levels_str))
          df = df.append({'family':task.family, 'tf':task.tf,
              'primer':task.primer, 'infile':infile}, ignore_index=True)
        print(df)
        outfile = '../seqlogos/d0/new.seqlogo.%s.allcycles.1.1.csv' % (shape_type)
        yield {
          'name'      : outfile,
          'actions'   : [(combine_csvs, [df, outfile])],
          'file_dep'  : df['infile'].tolist(),
          'targets'   : [outfile],
          'clean'     : True,
        }




def task_plot_seqlogo_pwms():
  """ Plot seqlogo for all cycles in both contexts (motif-free and
  motif-containing) in a TF experiment """
  for en_th in en_thresholds:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for task in task_infos:
            for motif, dist in izip(task.motifs, task.distances):
              outdir = '/'.join([top_seqlogo_dir, levels_type, en_th,
                  task.family, task.tf, task.primer])
              infile = "%s/%s.pwm" % (top_data_dir,
                  task.tf_info.getContextedShapemerFile(1111, 'all', 0, 0, 0,
                      'combined',shape_type, shape_levels_str))
              enrich = '%s/%s/%s/%s/%s/%s/enriched.%s.%s.%s.%d.%s.%d.%d.csv' % (
                  top_results_dir, levels_type, en_th, task.family, task.tf,
                  task.primer, task.tf, task.primer, shape_type, 4, motif,
                  lflank, rflank)
              outfile = "%s/%s" % (outdir, '.'.join(['seqlogo', task.tf,
                  task.primer, shape_type, "allcycles", motif, str(lflank),
                  str(rflank), 'pdf']))
              yield {
                'name'      : outfile,
                'actions'   : ["seqlogo_scripts/plot_seqlogo.R -i %s -e %s -t %s.%s.%s.%s -o %s" % (
                    infile, enrich, task.family, task.tf, task.primer, motif, outfile)],
                'file_dep'  : [infile, enrich],
                'targets'   : [outfile],
                'clean'     : True,
              }


def task_combine_seqlogos():
  """ Merge plots of seqlogos for all cycles in both contexts (motif-free and
  motif-containing) in all TF experiments """
  for en_th in en_thresholds:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          infiles = []
          for task in task_infos:
            for motif, dist in izip(task.motifs, task.distances):
              outdir = '/'.join([top_seqlogo_dir, levels_type, en_th,
                  task.family, task.tf, task.primer])
              infile = "%s/%s" % (outdir, '.'.join(['seqlogo', task.tf,
                  task.primer, shape_type, "allcycles", motif, str(lflank),
                  str(rflank), 'csv']))
              infiles.append("%s/%s" % (outdir, '.'.join(['seqlogo', task.tf,
                  task.primer, shape_type, "allcycles", motif, str(lflank),
                  str(rflank), 'pdf'])))
          outfile = '%s/fig_seqlogo_enriched_shapemers_%s_%s_th%s.pdf' % (
              top_results_dir, levels_type, fg_type, en_th)
          yield {
            'name'      : outfile,
            'actions'   : [(merge_pdfs_exclude_empty, [infiles, outfile])],
            'targets'   : [outfile],
            'file_dep'  : infiles,
            'clean'     : True,
          }

#def task_get_bg_shapemers_median_ic():
#  """ Generate shape mers """
#  for task in task_infos:
#    for shape_type in shapes:
#      for levels_type in ['publish']: #discrete_levels_type:
#        shinfo = shape_info[levels_type][shape_type]
#        shape_levels_str = shinfo.getLevelsStr()
#        for cycle in [4]: #task.cycles:
#          for motif, dist in izip(task.motifs, task.distances):
#              seqlogo_file = "%s/%s.pwm" % (top_data_dir, task.tf_info.getContextedShapemerFile(cycle, motif, dist, 0, 0, 'bg',shape_type, shape_levels_str))
#              median_ic_file = "%s/%s.mic" % (top_data_dir, task.tf_info.getContextedShapemerFile(cycle, motif, dist, 0, 0, 'bg',shape_type, shape_levels_str))
#              yield {
#                'name'      : median_ic_file,
#                'actions'   : ["analysis_scripts/compute_median_ic.R %s %s" %(seqlogo_file, median_ic_file)],
#                'file_dep'  : [seqlogo_file],
#                'targets'   : [median_ic_file],
#                'clean'     : True,
#              }
#
#
#
#
#def task_combine_seqlogo_median_ics():
#  """ Get summary information for enrichment """
#  for lflank, rflank in flank_configs:
#    for cycle in [4]: #task.cycles:
#      for levels_type in ['publish']:  #discrete_levels_type:
#        #shape_levels_str = shape_info[levels_type][shape_type].getLevelsStr()
#        df = pd.DataFrame(index = pd.MultiIndex.from_product([task_infos, shapes], names = ["task", "shape"]))
#        df = df.reset_index()
#        df['tf'] = df['task'].apply(lambda x: x.tf) 
#        df['primer'] = df['task'].apply(lambda x: x.primer) 
#        df['family'] = df['task'].apply(lambda x: x.family) 
#        df['motif'] = df['task'].apply(lambda x: x.motifs[0]) 
#        df['dist'] = df['task'].apply(lambda x: x.distances[0]) 
#        print(df.head())
#        df['infile'] = df.apply(lambda x: "%s/%s.mic" % (top_data_dir, x.task.tf_info.getContextedShapemerFile(4, x.motif, x.dist, 0, 0, 'bg', x['shape'], shape_info[levels_type][x['shape']].getLevelsStr())), axis='columns')
#
#        df = df[['shape', 'family','tf','primer','infile']]
#
#        outfile = "%s/%s/median_ics.%d.l%d.r%d.csv" % (top_results_dir, levels_type, cycle, lflank, rflank)
#        yield {
#          'name'      : outfile,
#          'actions'   : [(combine_csvs, [df, outfile])],
#          'file_dep'  : df['infile'].tolist(),
#          'targets'   : [outfile],
#          'clean'     : True,
#        }
#
