#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Clear all but the original data files from /mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/' """

import glob
import os
import shutil


data_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted'

#delete the master flats directory
master_flat = data_dir+'/master_flats'
if os.path.exists(master_flat):
	shutil.rmtree(master_flat)

#delete master bias files
bias_fpaths = glob.glob( data_dir + '/Night*/Bias/master_bias.fits')
for fpath in bias_fpaths:
	if os.path.exists(fpath):
		os.remove(fpath) 

#delete any files with 'ovscorr' in the name. Will remove all debiased and flattened files too
ovscorr_files = glob.glob( data_dir + '/Night*/*/*ovscorr*')
for line in ovscorr_files:
	os.remove(line)
	
print 'Files deleted'
	
	
