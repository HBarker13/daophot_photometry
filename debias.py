#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Use the master bias from each night to de-bias the flat and image data"""

import glob
from astropy.io import fits
import numpy as np
import os
import sys
from shutil import copyfile

data_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/'
night_dirs = glob.glob(data_dir+'/Night*')



#loop through the nights
for night in night_dirs:

	
	#check there is a master bias
	bias_fpath = night+'/Bias/master_bias.fits'
	if not os.path.exists(bias_fpath):
		print 'Master bias could not be found'
		print bias_fpath
		print 'Try master_bias.py'
		sys.exit()
		
		
		
	#read in the master bias	
	master_bias = fits.open(bias_fpath)
	bias = np.array(master_bias[0].data)
	master_bias.close()
	print 'Read in ', bias_fpath



	#list all the files that need de-biasing (all except the bias directory)
	img_dirs = glob.glob(night+'/*')
	img_dirs = [dirpath for dirpath in img_dirs if 'Bias' not in dirpath and 'README' not in dirpath]
	
	
	
	#loop though the files
	for dirname in img_dirs:
		
		#only use overscan corrected files, and no debiased files, and break if none are found
		fpaths = glob.glob(dirname+'/*ovscorr.fits')
		fpaths = [line for line in fpaths if 'debiased' not in line]
		if len(fpaths)==0:
			print 'No overscan corrected files found in ', dirname
			print 'Try overscan_correct.py'
			#sys.exit()
			raw_input('Press any key to continue')
			continue

		
		#loop through the files in each directory
		for fpath in fpaths:	
		
			
			#copy file to and update the new one
			newname, _ = fpath.rsplit('.fit')
			newname = newname + '_debiased.fits'  #change to xxx_ovscorr_debiased.fits
			
			#skip if the debiased file already exists
			if os.path.exists(newname):
				#delete the old file
				if os.path.exists(fpath):
					os.remove(fpath)
				continue
			
			copyfile(fpath, newname)
			
			open_img = fits.open(newname, mode='update')
			hdr = open_img[0].header
			img = np.array(open_img[0].data)
			
			#check the master bias and image are the same size:
			if np.shape(bias) != np.shape(img):
				print 'Bias does not have the same dimensions as the img'
				print 'Master bias: ', len(bias[0]), len(bias[1])
				print fpath
				print len(img[0]), len(img[1])
				#sys.exit()
				raw_input('Press to continue')
				continue
			
			#subtract the bias from the image
			debiased = np.subtract(img, bias)
			

			#save the debiased image to the file
			open_img[0].data = debiased
			
			#change header to add notes for overscan correction and trimming
			hdr.append( ('DEBIAS', 'True', 'Debiased'), end=True )
	
			#save the changes
			open_img.flush()
			open_img.close()
			print 'Created: ', newname
			print
				
				
			#delete the old file
			os.remove(fpath)	
		
	

	


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
