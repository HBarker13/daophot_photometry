#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Takes the overscan region, the two leftmost columns of the CCD, averages the value of these pixels and subtracts this value from all pixels in the CCD. Creates a new, overscan corrected file."""

#from pyraf import iraf
import os
import glob
from astropy.io import fits
import numpy as np
from shutil import copyfile


#read in files to be processed (ie. that haven't got ovscorr in name)
data_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/'
#all files (bias, flat, data) in all the Night dirs 
fpaths = glob.glob(data_dir+'/Night*/*/*.fits')
fpaths = [f for f in fpaths if 'ovscorr' not in f]


for fpath in fpaths:
	#print fpath
	
	#copy file and rename it
	newname, _ = fpath.rsplit('.fit')
	debias_name = newname + '_ovscorr_debiased.fits'
	newname = newname + '_ovscorr.fits'  #change to xxx_ovscorr.fits
	
	#skip any files that have already been corrected:
	if os.path.exists(newname):
		continue
	if os.path.exists(debias_name):
		continue	
	
	copyfile(fpath, newname)

	#open the file in update mode
	openfile = fits.open(newname, mode='update')
	hdr = openfile[0].header
	img = openfile[0].data
	img_arr = np.array(img)
	

	#create an array of the overscan region, the leftmost two columns
	ovsc_region = img_arr[:, [0,1] ]

	
	#calculate the mean of the overscan region
	ovsc_val = np.mean(ovsc_region)

	#subtract the mean from all pixels in the image
	corrected_img = img_arr - ovsc_val
	
	
	#save the corrected image to the fits file
	openfile[0].data = corrected_img

	
	#change header to add notes for overscan correction and trimming
	hdr.append( ('OVERSCAN', 'True', 'Overscan correction'), end=True )
	

	#save the changes
	openfile.flush()
	openfile.close()
	print 'Created: ', newname

print
print

















