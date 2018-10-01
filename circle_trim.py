#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Trim the science observations so only the circular filter aperture is left. MUST be run on files that have been trimmed (ie. 50 pix removed from each edge.
Assumes the circular filter is in the same place in all frames."""


import os
import glob
from astropy.io import fits
import numpy as np
from shutil import copyfile


data_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/'

#using Night2/FlatU/T2M3Im-20120516.2000428-0232.fits to deifne the aperture
#as there's a large difference between the circular filter and square frame
frame_fpath = data_dir + '/Night2/Flat_U/T2M3Im-20120516.200428-0232.fits'
frame = fits.getdata(frame_fpath)

#trim 70 pixels from each edge (like trim.py)
clipped = frame[70:len(frame[1])-70, 70:len(frame[0])-70 ]

#the sky (that we want to keep) has a value ~> 2900
#the rest of the frame is ~1600
#make an array of pixel coords with values >1900
clipped = np.array(clipped) > 1900

	
#read in files to be processed
fpaths = glob.glob(data_dir+'/Night*/*/*.fits')
fpaths = [f for f in fpaths if 'ovscorr_debiased_deflat.fits' in f]


for fpath in fpaths:
	
	"""
	#open the file in update mode
	openfile = fits.open(fpath, mode='update')
	hdr = openfile[0].header
	img = openfile[0].data
	img_arr = np.array(img)
	"""
	
	newname = fpath[:-5]+'_circletrim.fits'
	if os.path.exists(newname):
		continue
		
	copyfile(fpath, newname)
	openfile = fits.open(newname, mode='update')
	hdr = openfile[0].header
	img = openfile[0].data
	img_arr = np.array(img)
	
	#trim another 50 pixels
	#img_arr = img_arr[50:len(img_arr[1])-50, 50:len(img_arr[0])-50 ]
	

	circle = np.where(clipped, img_arr, float('nan')) #where clipped is true, keep the frame value
	

	#save the corrected image to the fits file
	openfile[0].data = circle

	
	#change header to add notes for overscan correction and trimming
	#hdr.append( ('TRIM', 'True', 'Ovscan region trimmed'), end=True )
	

	#save the changes
	openfile.flush()
	openfile.close()
	print 'Trimmed', fpath
	print

















