from dodo_common import *


task_infos = []

for i, row in tfs.iterrows():
  tf = row['tf']
  bc = row['primer']
  cycles = [row['final']] #[4] #[int(x) for x in row['cycles'].split(';')]
  #cycles = [1,2,3,4]
  motif = row['motif']
  family = row['family']
  dist = row['distance']
  accession = row['accession']
  task_infos.append(TaskInfo(tf, bc, family, accession, [motif], cycles, [dist]))

tmp = tfs[['family', 'tf', 'primer', 'motif', 'distance', 'final']]
cross = pd.merge(tmp, tmp, on=tmp.assign(key_col=1)['key_col'])
cross = cross.loc[cross['family_x'] != cross['family_y']]
#cross = cross.loc[~(((cross['family_x'] == 'homeodomain') & (cross['family_y'] == 'ETS')) | ((cross['family_y'] == 'homeodomain') & (cross['family_x'] == 'ETS')))]


def combine_data_frames(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (tf1, primer1, family1, tf2, primer2, family2) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['tf.x'] = tf1
    tmp['primer.x'] = primer1
    tmp['family.x'] = family1
    tmp['tf.y'] = tf2
    tmp['primer.y'] = primer2
    tmp['family.y'] = family2
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)


def combine_qvalues(infiles, infos, outfile):
  df = pd.DataFrame()
  for infile, (en_th, shape_type, levels_type, family, ctx) in izip(infiles, infos):
    tmp = pd.read_csv(infile)
    tmp['en_th'] = en_th
    tmp['shape'] = shape_type
    tmp['levels_type'] = levels_type
    tmp['family'] = family
    tmp['ctx'] = ctx
    df = df.append(tmp, ignore_index=True)
  df.to_csv(outfile, index=False)

def task_get_summary():
  """ Get summary information for enrichment """
  for en_th in en_thresholds:
    for task in task_infos:
      for lflank, rflank in flank_configs:
        for shape_type in shapes:
          for levels_type in discrete_levels_type:
            shinfo = shape_info[levels_type][shape_type]
            shape_levels_str = shinfo.getLevelsStr()
            for cycle in task.cycles:
              for motif, dist in izip(task.motifs, task.distances):
                  input_files = [top_data_dir +'/'+ task.tf_info.getContextedShapemerCountFile(cycle, motif, dist, lflank, rflank, ctx, shape_type, shape_levels_str) for ctx in contexts]
                  resdir = '/'.join([top_results_dir, levels_type, en_th, task.family, task.tf, task.primer])
                  file_fisher = '.'.join(['fisher', task.tf, task.primer, shape_type, str(cycle), motif, str(lflank), str(rflank), 'csv'])
                  output_file = '/'.join([resdir, file_fisher])
                  prob_files = [getCycle0ProbFile(top_data_dir, ctx, shape_type, motif, dist, lflank, rflank, shape_levels_str) for ctx in contexts]
                  yield {
                    'name'      : output_file,
                    'actions'   : ["results_scripts/get_summary.R %s %s %d %s %s %d %d %d %s %s %s %s %s" % (task.tf, task.primer, cycle, shape_type, motif, dist, lflank, rflank, en_th, resdir, shape_levels_str, task.family, top_data_dir)],
                    'file_dep'  : input_files + prob_files,
                    'targets'   : [output_file],
                    'clean'     : True,
                  }

def task_get_cross_summary():
  """ Get summary information for enrichment """
  for i, row in cross.iterrows():
    tf1, primer1, motif1, family1, dist1, cycle1 = row[['tf_x', 'primer_x', 'motif_x', 'family_x', 'distance_x', 'final_x']]
    tf2, primer2, motif2, family2, dist2, cycle2 = row[['tf_y', 'primer_y', 'motif_y', 'family_y', 'distance_y', 'final_y']]
    for en_th in en_thresholds:
      for lflank, rflank in flank_configs:
        for shape_type in shapes:
          for levels_type in discrete_levels_type:
            shinfo = shape_info[levels_type][shape_type]
            shape_levels_str = shinfo.getLevelsStr()
            resdir = '%s/%s/%s/%s_%s' % (top_results_dir, levels_type, en_th, min(family1, family2), max(family1, family2))
            info1 = TFInfo(tf1, primer1, family1)
            info2 = TFInfo(tf2, primer2, family2)
            input_files = [top_data_dir +'/'+ info1.getContextedShapemerCountFile(int(cycle1), motif1, int(dist1), lflank, rflank, 'fg', shape_type, shape_levels_str)]
            input_files += [top_data_dir +'/'+ info2.getContextedShapemerCountFile(int(cycle2), motif2, int(dist2), lflank, rflank, 'bg', shape_type, shape_levels_str)]
            file_fisher = '.'.join(['fisher', tf1, primer1, motif1, str(cycle1), tf2, primer2, motif2, str(cycle2), shape_type, str(lflank), str(rflank), 'csv'])
            output_file = '/'.join([resdir, file_fisher])
            prob_files = [getCycle0ProbFile(top_data_dir, 'fg', shape_type, motif1, dist1, lflank, rflank, shape_levels_str)]
            prob_files += [getCycle0ProbFile(top_data_dir, 'bg', shape_type, motif2, dist2, lflank, rflank, shape_levels_str)]
            yield {
              'name'      : output_file,
              'actions'   : ["results_scripts/get_cross_summary.R %s %s %s %s %d %d %s %s %s %s %d %d %s %d %d %s %s %s %s" % (tf1, primer1, family1, motif1, dist1, cycle1, tf2, primer2, family2, motif2, dist2, cycle2, shape_type, lflank, rflank, en_th, resdir, shape_levels_str, top_data_dir)],
              'file_dep'  : input_files + prob_files,
              'targets'   : [output_file],
              'clean'     : True,
            }

def task_same_enriched():
  for en_th in en_thresholds:
      for lflank, rflank in flank_configs:
        for shape_type in shapes:
          for levels_type in discrete_levels_type:
            inputs = ['/'.join([top_results_dir, levels_type, en_th, task.family, task.tf, task.primer, '.'.join(['enriched', task.tf, task.primer, shape_type, str(task.cycles[0]), motif, str(lflank), str(rflank), 'csv'])]) for task in task_infos for motif, dist in izip(task.motifs, task.distances)]
            infos = [(task.tf, task.primer, task.family, task.tf, task.primer, task.family) for task in task_infos for motif, dist in izip(task.motifs, task.distances)]
            target = '%s/%s/%s/denriched_%s.%s.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, 'same', shape_type, lflank, rflank)
            yield {
              'name'      : target,
              #'actions'   : ["results_scripts/fisher.R %s %s %s %d %d %s" % (shape_type, en_th, 'same', lflank, rflank, levels_type)],
              'actions'   : [(combine_data_frames, [inputs, infos, target])],
              'file_dep'  : inputs,
              'targets'   : [target],
              'clean'     : True,
            }


def task_same_fisher():
  for en_th in en_thresholds:
      for lflank, rflank in flank_configs:
        for shape_type in shapes:
          for levels_type in discrete_levels_type:
            inputs = ['/'.join([top_results_dir, levels_type, en_th, task.family, task.tf, task.primer, '.'.join(['fisher', task.tf, task.primer, shape_type, str(task.cycles[0]), motif, str(lflank), str(rflank), 'csv'])]) for task in task_infos for motif, dist in izip(task.motifs, task.distances)]
            infos = [(task.tf, task.primer, task.family, task.tf, task.primer, task.family) for task in task_infos for motif, dist in izip(task.motifs, task.distances)]
            target = '%s/%s/%s/dfisher_%s.%s.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, 'same', shape_type, lflank, rflank)
            yield {
              'name'      : target,
              #'actions'   : ["results_scripts/fisher.R %s %s %s %d %d %s" % (shape_type, en_th, 'same', lflank, rflank, levels_type)],
              'actions'   : [(combine_data_frames, [inputs, infos, target])],
              'file_dep'  : inputs,
              'targets'   : [target],
              'clean'     : True,
            }

def task_cross_fisher():
  combinations = []
  for i, row in cross.iterrows():
    tf1, primer1, motif1, family1, dist1, cycle1 = row[['tf_x', 'primer_x', 'motif_x', 'family_x', 'distance_x', 'final_x']]
    tf2, primer2, motif2, family2, dist2, cycle2 = row[['tf_y', 'primer_y', 'motif_y', 'family_y', 'distance_y', 'final_y']]
    resdir = '%s_%s' % (min(family1, family2), max(family1, family2))
    combinations.append((resdir, tf1, primer1, motif1, family1, dist1, cycle1, tf2, primer2, motif2, family2, dist2, cycle2))
  for en_th in en_thresholds:
      for lflank, rflank in flank_configs:
        for shape_type in shapes:
          for levels_type in discrete_levels_type:
            inputs = ['/'.join([top_results_dir, levels_type, en_th, resdir, '.'.join(['fisher', tf1, primer1, motif1, str(cycle1), tf2, primer2, motif2, str(cycle2), shape_type, str(lflank), str(rflank), 'csv'])]) for resdir, tf1, primer1, motif1, family1, dist1, cycle1, tf2, primer2, motif2, family2, dist2, cycle2 in combinations]
            infos = [(tf1, primer1, family1, tf2, primer2, family2) for resdir, tf1, primer1, motif1, family1, dist1, cycle1, tf2, primer2, motif2, family2, dist2, cycle2 in combinations]
            target = '%s/%s/%s/dfisher_%s.%s.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, 'cross', shape_type, lflank, rflank)
            yield {
              'name'      : target,
              #'actions'   : ["results_scripts/fisher.R %s %s %s %d %d %s" % (shape_type, en_th, 'cross', lflank, rflank, levels_type)],
              'actions'   : [(combine_data_frames, [inputs, infos, target])],
              'file_dep'  : inputs,
              'targets'   : [target],
              'clean'     : True,
            }

def task_compute_qvalue():
  for en_th in en_thresholds:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          for family in ['bHLH', 'ETS', 'homeodomain']:
            for exp_type in ['same', 'cross']:
              infile = '%s/%s/%s/dfisher_%s.%s.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, exp_type, shape_type, lflank, rflank)
              outfile = '%s/%s/%s/dqvalue_%s.%s.%s.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, exp_type, family, shape_type, lflank, rflank)
              yield {
                'name'      : outfile,
                'actions'   : ["results_scripts/compute_qvalue.R %s %s %s %s" % (shape_type, family, infile, outfile)],
                'file_dep'  : [infile],
                'targets'   : [outfile],
                'clean'     : True,
              }


def task_combine_all_qvalues():
  for lflank, rflank in flank_configs:
    infiles = ['%s/%s/%s/dqvalue_%s.%s.%s.l%d.r%d.csv' % (top_results_dir, levels_type, en_th, exp_type, family, shape_type, lflank, rflank)  for en_th in en_thresholds     for shape_type in shapes for levels_type in discrete_levels_type for family in ['bHLH', 'ETS', 'homeodomain'] for exp_type in ['same', 'cross']]
    infos = [(en_th, shape_type, levels_type, family, 'experiment' if (exp_type == 'same') else 'control')   for en_th in en_thresholds     for shape_type in shapes for levels_type in discrete_levels_type for family in ['bHLH', 'ETS', 'homeodomain'] for exp_type in ['same', 'cross']]
    outfile = '%s/dqvalue_l%d.r%d.csv' % (top_results_dir, lflank, rflank)
    yield {
      'name'      : outfile,
      'actions'   : [(combine_qvalues, [infiles, infos, outfile])],
      'file_dep'  : infiles,
      'targets'   : [outfile],
      'clean'     : True,
    }

