#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Create master flats for each filter on each night"""

import glob
from astropy.io import fits
import numpy as np
import os
import sys
from shutil import copyfile

data_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted'


"""There aren't flats for every night. I've kept the flats in directories according to the night they were made to avoid more confusion with the filters used on Night 1. I have kept the same labelling as was used in the observing log (ie. with the incorrect filter names.) The corrections are:
	U -> V
	B -> R
	V -> I
	R -> RG delta50 ?
	I -> Bessel z
I have deviated from the wording on the observing log and named all the flatfields 'flat', rather than 'sky'
"""

#master flat directory. 
#NB. In this directory, the flats are named correctly, ie. the Night1 flats have been labelled 
#using the filter ~really used, not the one recorded in the observing notes
master_dir = data_dir + '/master_flats'
if not os.path.exists(master_dir):
	os.makedirs(master_dir)


#expected dimensions of the flat arrays
len0 = 2148 - 140 #70 pixels taken off each edge in trim.py
len1 = 2148 - 140


#loop through each night
for nightnum in range(1,4):
	print 'Night', nightnum
	night_dir = data_dir + '/Night'+str(nightnum)
	
	flat_dirs = [line for line in glob.glob(night_dir+'/*') if 'Flat_' in line]
	

	#loop through each directory containing flat files
	for flat_dir in flat_dirs:
			
		#get the filter, remembering the changes for Night1
		_, filtername = flat_dir.rsplit('Flat_')
		if nightnum==1:
			if filtername=='U': filtername = 'V'
			elif filtername=='B': filtername = 'R'
			elif filtername=='V': filtername = 'I'
			elif filtername=='R': filtername = 'RGdelta50'
			elif filtername=='I': filtername = 'Besselz'
		
		print 'Filter', filtername
		
		#filename the master flat will be saved to
		savepath = master_dir + '/Night'+str(nightnum) + '_masterflat_' + filtername+'.fits'
		
		if os.path.exists(savepath):
			continue
			 
		
		#start with an empty flat array, and sum the arrays in each flat file
		summed_flat = []
		
		#use overscan corrected, debiased flat frames
		flat_fpaths = glob.glob(flat_dir+'/*ovscorr_debiased.fits')
		if len(flat_fpaths)==0:
			print 'No overscan corrected, debiased flat files found'
			print 'Try overscan_correct.py, debias.py'
			#import sys
			#sys.exit()
			raw_input('Press any key to continue')
			continue
		
		
		"""
		#make the master flat by calculating the mean	
		counter = 0
		for flat in flat_fpaths:
			open_flat = fits.open(flat)
			flat_data = open_flat[0].data
			open_flat.close()
		
			#sum the flat arrays
			if len(summed_flat)==0:
				summed_flat = flat_data
				
			else:
				summed_flat = np.add(summed_flat, flat_data)
				
			counter+=1
			
		#divide the calculated master_bias array by the total number of biases
		#ie. calculate the mean for each pixel
		master_flat = np.divide( np.array(summed_flat), counter)
		"""
		
		
		#make the master flat by calculating the median
		
		#make a list for each pixel
		holder = [ [ [] for y in xrange(len0)] for x in xrange(len1)]

		#loop through the flat frames and add to the list of values for each pixel
		for num, flat in enumerate(flat_fpaths):
	
			print 'Flat frame', num+1, '/', len(flat_fpaths)
	
			open_flat = fits.open(flat)
			flat_data = open_flat[0].data
			open_flat.close()
		
			for i in range(len1):
				for j in range(len0):
					holder[i][j].append(flat_data[i][j])

		#loop through each pixel in the holder and find the median for each list
		#Normalise the master flat array (so the average pixel value is 1.0)
		master_flat = [ [ [] for y in xrange(len0)] for x in xrange(len1)]
		
		median_sum = 0.0
		pixel_num = 0.0
	
		for i in range(len1):
			for j in range(len0):
			
				median = np.median( holder[i][j] )
				master_flat[i][j] = median
				
				#calculate the average value only from pixels
				#in a square in the region we care about.
				#Had to brute force it instead
	
				if median>1000:
					median_sum += median
					pixel_num += 1.0
				
		#mean value of the masterflat pixels = sum / number of pixels
		avg = median_sum / pixel_num
		
		print avg
				
		#divide the master flat by this average 
		normalised_master = np.divide( np.array(master_flat), avg)	

		
		#save the master flat
		hdu = fits.PrimaryHDU(normalised_master)
		hdu.writeto(savepath)
		print 'Saved: ', savepath	
		
		
	print	



#the combined I band master flat is made in deflat.py


















