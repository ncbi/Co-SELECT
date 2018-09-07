from dodo_common import *
from doit.tools import run_once

accession = pd.read_table(ena_project_file)
accession = accession[['run_accession', 'fastq_ftp']][0:2]
print accession

def task_downloads():
  for i, r in accession.iterrows():
      fastq_file = "%s/%s.fastq.gz" % (download_dir, r['run_accession'])
      url = r['fastq_ftp']
      yield {
        'name'      : fastq_file,
        'actions'   : ['wget -qO- {} > {}'.format(url, fastq_file)],
        'uptodate'  : [run_once],
        'targets'   : [fastq_file],
        'clean'     : True,
      }

