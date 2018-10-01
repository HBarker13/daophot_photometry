#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""script to create a master bias frame using overscan corrected bias frames taken on the same night"""

import os
import glob
from astropy.io import fits
import numpy as np


"""
Notes from observing logs:
	N1 bias taken towards the end of the night
	N2 beginning of the night: 'no sky flats -computer problem'
	N3 bias frames x40 at beginning of night and bias frames 'x20 end of night'
"""

data_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/'
night_dirs = glob.glob(data_dir+'/Night*/Bias')

#expected dimensions of the bias arrays
len0 = 2148-140 #70 pixels taken off each edge by trim.py
len1 = 2148-140  

master_bias = []

for night in night_dirs:
	
	
	savename = night+'/master_bias.fits'
	if os.path.exists(savename):
		continue
	
	#use overscan corrected bias files
	bias_fpaths = glob.glob(night+'/*ovscorr.fits')
	if len(bias_fpaths)==0:
		print 'No overscan corrected files found.'
		print 'Try: overscan_correct.py'
		import sys
		sys.exit()
		
	
	#calculate the median
	
	#make a list for each pixel
	holder = [ [ [] for y in xrange(len0)] for x in xrange(len1)]

	#loop through the bias frames and add to the list of values for each pixel
	for num, bias in enumerate(bias_fpaths):
	
		print 'Bias frame', num+1, '/', len(bias_fpaths)
	
		open_bias = fits.open(bias)
		bias_data = open_bias[0].data
		open_bias.close()
		
		for i in range(len1):
			for j in range(len0):
				holder[i][j].append(bias_data[i][j])

	#loop through each pixel in the holder and find the median for each list
	master_bias = [ [ [] for y in xrange(len0)] for x in xrange(len1)]
	
	for i in range(len1):
		for j in range(len0):
			median = np.median( holder[i][j] )
			master_bias[i][j] = median
		
		
		
		
	"""	
	#calculate the mean
	
	#start with empty array for the bias for each night
	summed_bias = []
	
	counter = 0
	for bias in bias_fpaths:
		open_bias = fits.open(bias)
		bias_data = open_bias[0].data
		open_bias.close()
		

			
		#sum all the arrays
		if len(summed_bias)==0:
			summed_bias = bias_data
		else:
			#add the master bias and new bias arrays
			summed_bias = np.add(summed_bias, bias_data)
		
		counter+=1


	#divide the calculated master_bias array by the total number of biases
	#ie. calculate the mean for each pixel
	#divisor_arr = np.full( (x_len, y_len), float(counter) )
	mean_bias = np.divide( np.array(summed_bias), counter)
	
	#check the difference between the mean and meadian biases
	diff = np.subtract(mean_bias, master_bias)
	"""
	
	

	#save the master_bias as a new fits file
	hdu = fits.PrimaryHDU(master_bias)
	hdu.writeto(savename, clobber=True)
	print 'Saved: ', savename
	print

		
		
		
		
		
		
		
		
		
		
		
			
		
		
		

