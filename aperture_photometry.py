#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""A wrapper to loop over all the 2.3m data and get aperture photometry using sextractor. Based on code sent by David Jones."""

#Support for tables
from astropy.table import Table
from astropy.table import Column

#Import fits support
from astropy.io import fits

#Import math functions
import numpy as np
import math
from math import sqrt

#Plotting
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

#Photometry in the style of sextractor
import sep
import os 
from pyexcel_ods import get_data
import glob
from functools import partial
import sys


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
gain=1.1 #gain in e-/adu


type_choice = None
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
	line_length = 9

elif type_choice=='cs':
	info_fpath = working_dir + '/Targets_info.ods'
	if not os.path.exists(info_fpath):
		print 'File cannot be found'
		print info_fpath
		sys.exit()
	line_length = 8


#read in the spreadhseet with the observation information
info = get_data(info_fpath)
fintab = []


#loop over the sheets (named after the standard stars/PN)
for sheetname in info:
	
	if sheetname=='Overview':
		continue
		
		
	if sheetname=='fake_wd1056': continue	
	if sheetname=='HaTr5': continue #No CS
	#if sheetname=='wd1056': continue
		
		
		
	print sheetname
	star_data = info[sheetname]
	#use shorted sheetname used for directory names

	
	colnames = [val.encode('utf-8') for val in star_data[0] ]
	

	#lines with all entries. Line length set in standard/cs choice
	star_data = [line for line in star_data[1:] if len(line)==line_length] 


	#make a recarray of each sheet
	data_arr = make_recarray(star_data, colnames)
	

	#loop over each entry in the sheet. ie. loop of info of each frame
	for obj in data_arr:
		print 'Frame', obj['Filename']


	
		#use the info to get to the correct directory and file
		night_dir = data_dir + '/' + obj['Night']
		frame_num = obj['Filename']
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
		for dirname in dirs:
			img_paths = glob.glob(dirname+'/*-'+frame_num+'_ovscorr_debiased_deflat_circletrim.fits')
			
			if len(img_paths)==0:
				continue
			elif len(img_paths)==1:
				img_fpath = img_paths[0]
			else:
				print 'Error: many files found'
				for l in img_paths:
					print l
				sys.exit()
				
		#try looking for another directory name		
		if img_fpath==None:
			sheetname, _ = sheetname.split('-')
			print sheetname	
			img_paths = glob.glob(dirname+'/*-'+frame_num+'_ovscorr_debiased_deflat_circletrim.fits')
			
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
			print sheetname, obj['Night'], frame_num
			print
			sys.exit()
			
			

		
		#read in the fits file
		hdu = fits.open( img_fpath )
		data = hdu[0].data
		
		#aperture size. Dave recommends 1.5*seeing. 
		fwhm = float( obj['fwhm'] )
		aper = 1.5*fwhm
		
		#Define sky annuli in terms of fwhm
		#ann_in_pix = 2*fwhm  #pixels
		#ann_out_pix = 3*fwhm
		
		
		#reorder the bytes: details at https://sep.readthedocs.io/en/v1.0.x/tutorial.html
		data = data.byteswap().newbyteorder()
		bkg = sep.Background(data) #calculate the background


		#Set a detection threshold for stars to be 1.5x RMS of background
		thresh = 1.5 * bkg.globalrms
		#Find all objects that are above threshold in a frame which is the background subtracted data
		objects = sep.extract(data-bkg.back(), thresh)
		
			
				
		#Table with list of object coords and fluxes
		stars=Table( [objects['x'],objects['y'],objects['flux']], names=['x','y','sepflux'])
		


		#Do the sky subtraction using the average for the frame
		#data_sub = sky subtracted data
		fluxes, fluxerrs, flag = sep.sum_circle( (data - bkg) , stars['x'], stars['y'], aper, err=bkg.globalrms, gain=gain)


		#Do sky subtraction using annuli
		#fluxes, fluxerrs, flag = sep.sum_circle( data, stars['x'], stars['y'], aper, err=bkg.globalrms, gain=gain, bkgann=(ann_in_pix,ann_out_pix) )


		#Add fluxes and uncertainties to table
		stars.add_column( Column( fluxes, name='flux') )
		stars.add_column( Column( fluxerrs, name='flux_err') )


		#Calculate instrumental magnitudes
		expt = float( obj['Exp_time'] )
		mags = [ -2.5*math.log10( line / expt) if line>0 else float('nan') for line in fluxes ]
		stars.add_column( Column( mags, name='mag') )
		
		
		#Calculate instrumental magnitude errors
		up_errs = [ abs( -2.5*math.log10( (line[0]+line[1]) / expt) - line[2] ) if (line[0]+line[1])>0 else float('nan') for line in zip(fluxes, fluxerrs, mags) ]

		low_errs = [ abs( -2.5*math.log10( (line[0]-line[1]) / expt) - line[2] ) if (line[0]-line[1])>0 else float('nan') for line in zip(fluxes, fluxerrs, mags) ]
		
		avg_errs = [ (line[0]+line[1])/2. for line in zip(up_errs, low_errs) ]
		stars.add_column( Column( avg_errs, name='mag_err') )
		


		#find the target stars's coordinates in the stars table	
		stars_coords = [ [ind, line['x'], line['y']] for ind,line in enumerate(stars)]
		target_coords = [ float( obj['x_coord']), float( obj['y_coord']) ]
		
		
		#function to calculate distance between coords
		dist=lambda true, test: ( abs( true[0]-test[1] )**2 + abs( true[1]-test[2] )**2 )
		nearest = min(stars_coords, key=partial(dist, target_coords) ) 


		#print target_coords
		#print nearest

		"""
		#check the identified star is within a reasonable distance of where it should be
		difference_val = 4.
		match_star_flag = False
		if abs(nearest[1]-float( obj['x_coord']))>difference_val or abs(nearest[2]-float( obj['y_coord']))>difference_val:
			print 'Standard star coordinates could not be found'
			newline = [ sheetname, obj['Night'], obj['Filename'], obj['Filter'], 'Coords not found']
			fintab.append(newline)
			match_star_flag = True
			continue
		if match_star_flag==True: continue
		"""



		#find the star with the nearest coordinates in the mags_list to get magnitudes
		object_line = stars[ nearest[0] ]


		"""
		#plot the location of the objects on the image
		fig, ax = plt.subplots()
		im = ax.imshow(data, interpolation='nearest', cmap='gray', vmin= bkg.globalback-bkg.globalrms, vmax=bkg.globalback+bkg.globalrms, origin='lower')
		for i, row in enumerate(objects):
			if row['x']==object_line['x'] and row['y']==object_line['y']:
				e = Ellipse( xy=(objects['x'][i], objects['y'][i]), width=6*objects['a'][i], height=6*objects['b'][i], angle=objects['theta'][i]*180./np.pi)
				e.set_facecolor('none')
				e.set_edgecolor('red')
				ax.add_artist(e)
		plt.show()
		"""
		




		
		#pn, night, frame number, filter, exposure time, airmass, mag, mag_err
		newline = [ sheetname, obj['Night'], obj['Filename'], obj['Filter'],  obj['Exp_time'], obj['Airmass'], fwhm, obj['x_coord'], obj['y_coord'], object_line['mag'], object_line['mag_err'] ]
		fintab.append(newline)
		
		print newline
		print 'Line added'
		print
		
		


#save the table to file
print 'Saving'
if type_choice=='s': fin_fpath = working_dir + '/aper_standard_mags.tab'
elif type_choice=='cs': fin_fpath = working_dir + '/aper_cs_mags.tab'
	
with open(fin_fpath, 'w+') as f:
	f.write('star	night	frame_number	filter	exposure_time	airmass	fwhm	x	y	mag	mag_err \n')
	for line in fintab:
		print_line = ''
		for s in line:
			print_line+=str(s)
			print_line+='	'
		print_line+='\n'
		f.write(print_line)	

print 'File written'
print fin_fpath
print


	
	
	
	
	
	
	
	
