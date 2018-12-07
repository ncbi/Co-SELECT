from dodo_common import *

task_infos = []

#tfs = tfs[tfs['tf'] == 'MAX']

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


def task_preprocess():
  """ Unzip fastq files, keep only sequence info of those containing only ACGT """
  for task in task_infos:
    for cycle in task.cycles:
      fastq_file = "%s/%s.fastq.gz" % (download_dir, task.accession)
      seq_file = "%s/%s" % (top_data_dir, task.tf_info.getSequenceFile(cycle))
      count_file = "%s.cnt" % (seq_file)
      ensure_dir(seq_file)
      yield {
        'name'      : seq_file,
        'actions'   : [(unzip_seq_filter_N, [task.tf_info.primer, fastq_file, seq_file, count_file])],
        'file_dep'  : [fastq_file],
        'targets'   : [seq_file, count_file],
        'clean'     : True,
      }

def task_get_shape():
  """ Generate MGW values from input apramer partitions using DNAShape program """
  for task in task_infos:
    for shape_type in shapes:
      for cycle in task.cycles:
          seq_file = task.tf_info.getSequenceFile(cycle)
          input_file = "%s/%s" % (top_data_dir, seq_file)
          output_file = input_file + "." + shape_type
          yield {
            'name'      : output_file,
            'actions'   : ["%s %s %s" % (dnashape_exe, input_file, shape_type)],
            'file_dep'  : [input_file],
            'targets'   : [output_file],
            'clean'     : True,
          }
 
def task_discretize_shape():
  """ Discretize MGW values obtained using DNAShape program """
  for task in task_infos:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        for cycle in task.cycles:
            seq_file = task.tf_info.getSequenceFile(cycle)
            input_file = "%s/%s.%s" % (top_data_dir, seq_file, shape_type)
            output_file = "%s/%s" % (top_data_dir, task.tf_info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
            yield {
              'name'      : output_file,
              'actions'   : [(task.tf_info.discretize_shape, [shinfo, input_file, output_file])],
              'file_dep'  : [input_file],
              'targets'   : [output_file],
              'clean'     : True,
            }


def task_partition():
  """ Partition aptamer sequences into foreground and background based on distance from MOTIF """
  for task in task_infos:
    for cycle in task.cycles:
      seq_file = "%s/%s" % (top_data_dir, task.tf_info.getSequenceFile(cycle))
      for motif, dist in izip(task.motifs, task.distances):
          fg_file = "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'fg'))
          bg_file = "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'bg'))
          yield {
            'name'      : ':'.join([seq_file, motif, str(dist)]),
            'actions'   : [(task.tf_info.partition_aptamers, [seq_file, motif, dist, fg_file, bg_file])],
            'file_dep'  : [seq_file],
            'targets'   : [fg_file, bg_file],
            'clean'     : True,
          }
      

def task_get_fg_parts():
  """ Generate shape mers """
  for task in task_infos:
    for lflank, rflank in flank_configs:
      for cycle in task.cycles:
        seq_file = "%s/%s" % (top_data_dir, task.tf_info.getSequenceFile(cycle))
        for motif, dist in izip(task.motifs, task.distances):
            fg_file = "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'fg'))
            parts_file = "%s/%s" % (top_data_dir, task.tf_info.getFgPartsFile(cycle, motif, dist, lflank, rflank))
            yield {
              'name'      : parts_file,
              'actions'   : [(task.tf_info.gen_fg_parts, [seq_file, fg_file, motif, parts_file])],
              'file_dep'  : [seq_file, fg_file],
              'targets'   : [parts_file],
              'clean'     : True,
            }

def task_get_fg_shapemers():
  """ Generate shape mers """
  for task in task_infos:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in task.cycles:
            seq_file = "%s/%s" % (top_data_dir, task.tf_info.getSequenceFile(cycle))
            shape_file = "%s/%s" % (top_data_dir, task.tf_info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
            count_file = "%s.cnt" % (seq_file)
            for motif, dist in izip(task.motifs, task.distances):
                context_file = "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'fg'))
                parts_file = "%s/%s" % (top_data_dir, task.tf_info.getFgPartsFile(cycle, motif, dist, lflank, rflank))
                shapemer_file = "%s/%s" % (top_data_dir, task.tf_info.getContextedShapemerFile(cycle, motif, dist, lflank, rflank, 'fg',shape_type, shape_levels_str))
                yield {
                  'name'      : shapemer_file,
                  'actions'   : [(task.tf_info.gen_fg_shapemers, [shinfo, task.shape_length, seq_file, shape_file, count_file, context_file, motif, lflank, rflank, parts_file, shapemer_file])],
                  'file_dep'  : [shape_file, count_file, context_file, parts_file],
                  'targets'   : [shapemer_file],
                  'clean'     : True,
                }


def task_get_bg_shapemers():
  """ Generate shape mers """
  for task in task_infos:
    for shape_type in shapes:
      for levels_type in discrete_levels_type:
        shinfo = shape_info[levels_type][shape_type]
        shape_levels_str = shinfo.getLevelsStr()
        for cycle in task.cycles:
          seq_file = "%s/%s" % (top_data_dir, task.tf_info.getSequenceFile(cycle))
          shape_file = "%s/%s" % (top_data_dir, task.tf_info.getDiscreteShapeFile(cycle, shape_type, shape_levels_str))
          count_file = "%s.cnt" % (seq_file)
          for motif, dist in izip(task.motifs, task.distances):
              context_file = "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'bg'))
              shapemer_file = "%s/%s" % (top_data_dir, task.tf_info.getContextedShapemerFile(cycle, motif, dist, 0, 0, 'bg',shape_type, shape_levels_str))
              yield {
                'name'      : shapemer_file,
                'actions'   : [(task.tf_info.gen_bg_shapemers, [shinfo, task.shape_length, shape_file, count_file, context_file, shapemer_file])],
                'file_dep'  : [shape_file, count_file, context_file],
                'targets'   : [shapemer_file],
                'clean'     : True,
              }

def task_count_fg_shapemers():
  """ Count shape mers """
  for task in task_infos:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in task.cycles:
            for motif, dist in izip(task.motifs, task.distances):
                shapemer_file = "%s/%s" % (top_data_dir, task.tf_info.getContextedShapemerFile(cycle, motif, dist, lflank, rflank, 'fg',shape_type, shape_levels_str))
                shapemer_count_file = "%s/%s" % (top_data_dir, task.tf_info.getContextedShapemerCountFile(cycle, motif, dist, lflank, rflank, 'fg',shape_type, shape_levels_str))
                yield {
                  'name'      : shapemer_count_file,
                  'actions'   : ["cat %s | cut -d' ' -f1 | sort | uniq -c > %s" % (shapemer_file, shapemer_count_file)],
                  'file_dep'  : [shapemer_file],
                  'targets'   : [shapemer_count_file],
                  'clean'     : True,
                }

def task_count_bg_shapemers():
  """ Count shape mers """
  for task in task_infos:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in task.cycles:
            for motif, dist in izip(task.motifs, task.distances):
                shapemer_file = "%s/%s" % (top_data_dir, task.tf_info.getContextedShapemerFile(cycle, motif, dist, 0, 0, 'bg',shape_type, shape_levels_str))
                shapemer_count_file = "%s/%s" % (top_data_dir, task.tf_info.getContextedShapemerCountFile(cycle, motif, dist, 0, 0, 'bg',shape_type, shape_levels_str))
                yield {
                  'name'      : shapemer_count_file,
                  'actions'   : ["cat %s | cut -d' ' -f1 | sort | uniq -c > %s" % (shapemer_file, shapemer_count_file)],
                  'file_dep'  : [shapemer_file],
                  'targets'   : [shapemer_count_file],
                  'clean'     : True,
                }

def task_get_fg_coverage():
  """ Get coverage information for each shape mer """
  for task in task_infos:
    for lflank, rflank in flank_configs:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in task.cycles:
            for motif, dist in izip(task.motifs, task.distances):
                shapemer_file = "%s/%s" % (top_data_dir, task.tf_info.getContextedShapemerFile(cycle, motif, dist, lflank, rflank, 'fg',shape_type, shape_levels_str))
                cov_file = shapemer_file + ".cov"
                context_file = "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'fg'))
                count_file = "%s/%s.cnt" % (top_data_dir, task.tf_info.getSequenceFile(cycle))
                yield {
                  'name'      : cov_file,
                  'actions'   : [(getCoverage, [context_file, count_file, shapemer_file, cov_file])],
                  'file_dep'  : [shapemer_file, count_file, context_file],
                  'targets'   : [cov_file],
                  'clean'     : True,
                }

def task_get_bg_coverage():
  """ Get coverage information for each shape mer """
  for task in task_infos:
      for shape_type in shapes:
        for levels_type in discrete_levels_type:
          shinfo = shape_info[levels_type][shape_type]
          shape_levels_str = shinfo.getLevelsStr()
          for cycle in task.cycles:
            for motif, dist in izip(task.motifs, task.distances):
                shapemer_file = "%s/%s" % (top_data_dir, task.tf_info.getContextedShapemerFile(cycle, motif, dist, 0, 0, 'bg',shape_type, shape_levels_str))
                cov_file = shapemer_file + ".cov"
                context_file = "%s/%s" % (top_data_dir, task.tf_info.getContextFile(cycle, motif, dist, 'bg'))
                count_file = "%s/%s.cnt" % (top_data_dir, task.tf_info.getSequenceFile(cycle))
                yield {
                  'name'      : cov_file,
                  'actions'   : [(getCoverage, [context_file, count_file, shapemer_file, cov_file])],
                  'file_dep'  : [shapemer_file, count_file, context_file],
                  'targets'   : [cov_file],
                  'clean'     : True,
                }

