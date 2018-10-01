#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Trim overscan corrected images to remove the overscan region, plus a few pixels"""


import os
import glob
from astropy.io import fits
import numpy as np
from shutil import copyfile


#read in files to be processed
data_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/'
fpaths = glob.glob(data_dir+'/Night*/*/*.fits')
fpaths = [f for f in fpaths if 'ovscorr' in f]


for fpath in fpaths:
	
	#open the file in update mode
	openfile = fits.open(fpath, mode='update')
	hdr = openfile[0].header
	img = openfile[0].data
	img_arr = np.array(img)
	
	#skip if there is header flag that the frame has already been trimmed
	if 'TRIM' in hdr:
		openfile.close()
		continue
	
	
	clipped = img_arr[70:len(img_arr[1])-70, 70:len(img_arr[0])-70 ] 

	#save the corrected image to the fits file
	openfile[0].data = clipped

	
	#change header to add notes for overscan correction and trimming
	hdr.append( ('TRIM', 'True', 'Ovscan region trimmed'), end=True )
	

	#save the changes
	openfile.flush()
	openfile.close()
	print 'Trimmed', fpath
	print

















