#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#Use photutils (which uses daophot) to get psf photometry


from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.psf import prepare_psf_model
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.io import fits
from astropy.table import Table
import os
import glob
from pyexcel_ods import get_data
import numpy as np

import matplotlib.pyplot as plt

#takes an uncollimated table and converts into recarray
def make_recarray(tab, title_list):
	dtype_list = ['|S20' for item in title_list]
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]

	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]	
		data_array.append(col)
		
	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array	






#top directory paths
working_dir =  '/mirror2/scratch/hbarker/Orsola_2.3m_ANU'
data_dir = working_dir+'/Sorted'

type_choice = 's'
#type_choice = None
type_options = ['s', 'cs']
while type_choice not in type_options:
	type_choice = raw_input('Process standard stars (s) or CS targets (cs)? ')
	type_choice = type_choice.lower()


if type_choice=='s':
	#read in Standards_info.ods to get information about the
	#observing data and loop through all the files
	info_fpath = working_dir + '/Standards_info.ods'
	if not os.path.exists(info_fpath):
		print 'File cannot be found'
		print info_fpath
		sys.exit()
		
	fin_fpath = working_dir + '/photutils_standard_mags.tab'	


elif type_choice=='cs':
	info_fpath = working_dir + '/Targets_info.ods'
	if not os.path.exists(info_fpath):
		print 'File cannot be found'
		print info_fpath
		sys.exit()
		
	fin_fpath = working_dir + '/photutils_cs_mags.tab'



#make a recarray of each sheet
info = get_data(info_fpath)



fintab = []



for sheetname in info:

	if sheetname=='Overview': continue		
	if sheetname=='fake_wd1056': continue	


	if sheetname!='LSE44': continue

	#if sheetname!='121-968': continue

	
	
	
	print sheetname	
	star_data = info[sheetname]
	colnames = [val.encode('utf-8') for val in star_data[0] ]
	data = [line for line in star_data[1:] if len(line)==9] #lines with all entries.
	
	#make a recarray of each sheet
	data_arr = make_recarray(data, colnames)	


	for star_info in data_arr:
	

	
		#if star_info['Night']!='Night3': continue
		#if star_info['Filter']!='U': continue	
		#if float(star_info['Filename'])>200: continue
	
	
		print star_info['Night']
		print 'Frame', star_info['Filename']
		
		
		

		#get frame details from the spreadsheet
		fwhm = float( star_info['fwhm'])
		exptime = float(star_info['Exp_time'])
		
		if star_info['sky_sigma'] =='x': 
			skysigma = 4.0
		else:
			skysigma = float(star_info['sky_sigma'])

		xcoord = float(star_info['x_coord'])
		ycoord = float(star_info['y_coord'])
		
		
		
		#use the info to get to the correct directory and file
		night_dir = data_dir + '/' + star_info['Night']
		frame_num = star_info['Filename']
		#add leading zeros
		if len(frame_num)==3:
			frame_num = '0'+frame_num
		elif len(frame_num)==2:
			frame_num = '00'+frame_num
		
			
		
		#Stars with multiple observations on the same night have multiple directories
		#so we can't know the directory name ahead of time
		if type_choice=='s': dirs = glob.glob(night_dir+'/st_'+sheetname+'*')
		elif type_choice=='cs': dirs = glob.glob(night_dir+'/'+sheetname+'*')
		
		img_fpath = None
		for standard_dir in dirs:
			img_paths = glob.glob(standard_dir+'/*-'+frame_num+'_ovscorr_debiased_deflat_circletrim.fits')
			
			if len(img_paths)==0:
				continue
			elif len(img_paths)==1:
				img_fpath = img_paths[0]
			else:
				print 'Error: many files found'
				for l in img_paths:
					print l
				sys.exit()
				
		if img_fpath==None:
			print 'File not found'
			print
			print sheetname, star_info['Night'], frame_num
			print
			sys.exit()


		image = fits.getdata(img_fpath)

	
		"""
		from photutils.datasets import (make_random_gaussians_table, make_noise_image, make_gaussian_sources_image)
		sigma_psf = 2.0
		sources = Table()
		sources['flux'] = [700, 800, 700, 800]
		sources['x_mean'] = [12, 17, 12, 17]
		sources['y_mean'] = [15, 15, 20, 20]
		sources['x_stddev'] = sigma_psf*np.ones(4)
		sources['y_stddev'] = sources['x_stddev']
		sources['theta'] = [0, 0, 0, 0]
		sources['id'] = [1, 2, 3, 4]
		tshape = (32, 32)
		fwhm = sigma_psf*gaussian_sigma_to_fwhm
		image = (make_gaussian_sources_image(tshape, sources) + make_noise_image(tshape, type='poisson', mean=6., random_state=1) + make_noise_image(tshape, type='gaussian', mean=0., stddev=2., random_state=1))
		"""
		
		#plt.imshow(image, cmap='viridis', aspect=1, interpolation='nearest', origin='lower') 
		#plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04) 
		#plt.show()
		
		
		
		

		#calculate background
		bkgrms = MADStdBackgroundRMS()
		std = bkgrms(image)
		print 'Background level :', std
		
		
		#find stars
		print 'Finding stars'
		iraffind = IRAFStarFinder(threshold=3.5*std, fwhm=fwhm, minsep_fwhm=0.01, roundhi=5.0, roundlo=-5.0, sharplo=0.0, sharphi=2.0)
		
			
		daogroup = DAOGroup(2.0*fwhm)
		
		mmm_bkg = MMMBackground()
		fitter = LevMarLSQFitter()
		psf_model = IntegratedGaussianPRF(sigma=(fwhm/gaussian_sigma_to_fwhm))
		
		#fitshape must be odd
		fitshape = 2*int(fwhm)+1
	
		print 'Performing photometry'
		#photometry = IterativelySubtractedPSFPhotometry(finder=iraffind, group_maker=daogroup, bkg_estimator=mmm_bkg, psf_model=psf_model,fitter=LevMarLSQFitter(), niters=None, fitshape=(fitshape, fitshape), aperture_radius=fwhm)
		#niters = None means continue until all stars subtracted - might recur infinitely
		
		
		from photutils.psf import DAOPhotPSFPhotometry
		photometry =  DAOPhotPSFPhotometry(crit_separation=3*fwhm, threshold=3.5*std, fwhm=fwhm, psf_model=psf_model, fitshape=(fitshape, fitshape), sharplo=0.0, sharphi=2.0, roundlo=-5.0, roundhi=5.0, fitter=iraffind, niters=100, aperture_radius=1.2*fwhm)
		raw_input('Done')
		
		print 'Results table'
		result_tab = photometry(image=image)
		
		raw_input( 'Creating residual image')
		residual_image = photometry.get_residual_image()
		
		print 'Plotting'
		plt.subplot(1, 2, 1)
		plt.imshow(image, cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
		plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
		plt.subplot(1 ,2, 2)
		plt.imshow(residual_image, cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
		plt.title('Residual Image')
		plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
		plt.show()


print
print '---------------Completed---------------'
print











