dnashape_exe              = '../DNAshapeR/src/dnashape'
ena_project_file          = './PRJEB14744.txt'
nonzero_accession_file    = './PRJEB14744_nonzero_cycle.csv'
zero_accession_file       = './PRJEB14744_zero_cycle.csv'
tf_info_file              = './tf_inventory_jolma_ronshamir.csv'
tf_motif_file             = './tf_coremotif.csv'
tf_run_coselect_file      = './tf_run_coselect.csv'

orig_data_dir             = '../data'
seqmer_data_dir           = '../seqmerdata'
download_dir              = '../downloads'

contexts                  = ['fg', 'bg']
fg_types                  = ['d0', 'd1all', 'd1enriched']
fg_type                   = 'd1all'
fg_type                   = 'd0'
fg_type                   = 'd1enriched'

cycles                    = [3, 4] if fg_type == 'd0' else [4]

shapes                    = ["MGW", "ProT", "HelT", "Roll"]

flank_configs             = [(1,1)]

discrete_levels_type      = ['publish', 'other1', 'other2']

en_thresholds             = ['1.05', '1.10', '1.20', '1.50']

discrete_levels_type      = ['publish', 'other1']

top_data_dirs = {
  'd0'         : '../data',
  'd1all'      : 'xxx',
  'd1enriched' : '../derived'
}

probability_dirs = {
  'd0'         : '../data/simple_predict',
  'd1all'      : 'xxx',
  'd1enriched' : '../derived/simple_predict'
}

top_data_dir              = top_data_dirs[fg_type]
probability_dir           = probability_dirs[fg_type]
top_results_dir           = '../results/{}'.format(fg_type)


print top_data_dir
print probability_dir
print top_results_dir

