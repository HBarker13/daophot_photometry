#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Use the master flat in each filter to de-flatten each of the science images"""

import glob
from astropy.io import fits
import numpy as np
import os
import sys
from shutil import copyfile



data_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted'
master_flat_dir = data_dir + '/master_flats'
night_dirs = glob.glob(data_dir+'/Night*')

#list of all the master flat files
master_fpaths = glob.glob(master_flat_dir+'/*.fits')



#use the Night1 + Night3 I band flat
combined_I_fpath = master_flat_dir+'/combined_masterflat_I.fits'

if not os.path.exists(combined_I_fpath):
	copyfile( master_flat_dir+'/Night1_masterflat_I.fits', combined_I_fpath)
	copied = fits.open(combined_I_fpath, mode='update')
	master_1 = copied[0].data
	
	master = fits.open( master_flat_dir+'/Night3_masterflat_I.fits')
	master_3 = np.array( master[0].data )
	
	summed = np.add(master_1, master_3)
	newimg = np.divide(summed, 2)
	copied[0].data = newimg
	copied.flush()
	copied.close()	
	
master = fits.open(combined_I_fpath)
master_i = master[0].data
master.close()	



#read in all the master flats: fasters than looping over all the exposures, opening them
#and checking they have the correct filter

for master_flat in master_fpaths:

	#find the filter of the master flat
	_, filtername = master_flat.rsplit('masterflat_')
	filtername = filtername[:-5] #remove '.fits'

	master = fits.open(master_flat)
	master_arr = np.array(master[0].data)
	master.close()
	
	if filtername=='U': master_u = master_arr
	elif filtername=='B': master_b = master_arr
	elif filtername=='V': master_v = master_arr
	elif filtername=='R': master_r = master_arr
	elif filtername=='Bessellz': master_bessell_z = master_arr	

		
	
	

#loop though nights
for night in night_dirs:
	
	obj_dirs = [line for line in glob.glob(night+'/*') if 'Bias' not in line and 'Flat' not in line]
		
	#choose a directory/object
	for obj in obj_dirs:
		
		#list overscan corrected, debiased files
		exposures = glob.glob(obj+'/*ovscorr_debiased.fits')
			
		for fpath in exposures:
			print fpath
		
			flat = None
			
			#get exposure number because the first 70 frames of Night2 are wrong, 
			#like all of night 1
			#format = T2M3Im-date.XXXXX-framenumber_ovscorr_debias.fits
			_, date, expnum = fpath.rsplit('-', 2)
			expnum, _ = expnum.split('_', 1)
			expnum = int(expnum)
			

			
			#copy the file into the deflattened filename
			savepath = fpath[:-5] + '_deflat.fits'
			
			#skip if the deflattened file already exists
			if os.path.exists(savepath):
				continue
			
			copyfile(fpath, savepath)
			
			#open the copy and find out the filter
			exposure = fits.open(savepath, mode='update')
			hdr = exposure[0].header
			img = np.array(exposure[0].data)
			
					
			
			#all data seems to be in header['object']: name_filter_exptime
			filtername = None
			try:				
				name, filtername, exptime = hdr['object'].split('_')
			except:
				pass
				
			try:
				name, filtername, exptime = hdr['object'].split(' ')
			except:
				pass
				
			try:
				name, filtername = hdr['object'].split('_')
			except:
				pass
				
				
			""" Very ugly, but will fix all the wrong header filter labels here. NOTE. The filters for Night1 are as
			they are listed in the observing log. ie. they are incorrect"""
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/st_wd1056_exp34-46/T2M3Im-20120515.112830-0034_ovscorr_debiased.fits':
				filtername='U'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.091312-0019_ovscorr_debiased.fits':
				filtername='V'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.091422-0020_ovscorr_debiased.fits':
				filtername='V'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.092129-0021_ovscorr_debiased.fits':
				filtername='V'
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.103606-0022_ovscorr_debiased.fits':
				filtername='B'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.103652-0023_ovscorr_debiased.fits':
				filtername='B'	
			
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.104142-0025_ovscorr_debiased.fits':
				filtername='I'
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.104333-0026_ovscorr_debiased.fits':
				filtername='U'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.104449-0027_ovscorr_debiased.fits':
				filtername='U'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.104057-0024_ovscorr_debiased.fits':
				filtername='I'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.105638-0028_ovscorr_debiased.fits':
				filtername='I'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.105805-0029_ovscorr_debiased.fits':
				filtername='R'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.111354-0030_ovscorr_debiased.fits':
				filtername='Halpha'	
		
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.111529-0031_ovscorr_debiased.fits':
				filtername='Halpha'
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.112002-0032_ovscorr_debiased.fits':
				filtername='I'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night1/Lo5/T2M3Im-20120515.112133-0033_ovscorr_debiased.fits':
				filtername='I'	
				
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night2/st_wd1056_exp66-88/T2M3Im-20120516.111511-0080_ovscorr_debiased.fits':
				filtername='U'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night2/st_wd1056_exp66-88/T2M3Im-20120516.111223-0078_ovscorr_debiased.fits':
				filtername='U'
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night2/st_wd1056_exp66-88/T2M3Im-20120516.105112-0071_ovscorr_debiased.fits':
				filtername='I'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night2/st_wd1056_exp66-88/T2M3Im-20120516.111626-0081_ovscorr_debiased.fits':
				filtername='U'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night2/st_wd1056_exp66-88/T2M3Im-20120516.110657-0075_ovscorr_debiased.fits':
				filtername='I'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night2/st_wd1056_exp66-88/T2M3Im-20120516.110852-0077_ovscorr_debiased.fits':
				filtername='I'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night2/st_wd1056_exp66-88/T2M3Im-20120516.110816-0076_ovscorr_debiased.fits':
				filtername='I'	
				
			if fpath=='/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Sorted/Night2/st_wd1056_exp66-88/T2M3Im-20120516.111356-0079_ovscorr_debiased.fits':
				filtername='U'	


			#night 1 and the first 70 night 2 exposures are wrong
			#U -> V
			#B -> R
			#V -> I
			#R -> RG delta50 ?
			#I -> Bessel z
			if filtername=='Halpha': continue
			
			if 'Night1' in night:
				if filtername == 'U': flat = master_v
				if filtername == 'B': flat = master_r
				if filtername == 'V': flat = master_i
				if filtername == 'R': continue #master_rgdelta
				if filtername == 'I': flat = master_bessell_z
			
			if 'Night2' in night and expnum<71:
				if filtername == 'U': flat = master_v
				if filtername == 'B': flat = master_r
				if filtername == 'V': flat = master_i
				if filtername == 'R': continue #master_rgdelta
				if filtername == 'I': flat = master_bessell_z
	
			else:
				if filtername == 'U': flat = master_u
				if filtername == 'B': flat = master_b
				if filtername == 'V': flat = master_v
				if filtername == 'R': flat = master_r
				if filtername == 'I': flat = master_i



			#check the master flat and image are the same size:
			if np.shape(flat) != np.shape(img):
				print filtername
				print fpath
			
				print 'Master flat does not have the same dimensions as the img'
				print fpath
				print 'Img: ', len(img[0]), len(img[1])
				print 'Master flat: ', len(flat[0]), len(flat[1])
				#sys.exit()
				raw_input('Press any key to continue')
				
				
			#divide the img by the masterflat
			deflattened = np.divide(img, flat)
			
			#save the image
			exposure[0].data = deflattened	

			#add a header flag
			hdr.append( ('DEFLAT', 'True', 'De-flattened'), end=True )
			
			
			#save the changes
			exposure.flush()
			exposure.close()
			print 'Created: ', savepath
			print
				
				
			#delete the old file
			os.remove(fpath)	
			
			
			
			
			
			
			
			
			
			
			
			
			
