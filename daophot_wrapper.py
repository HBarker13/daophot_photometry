#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""A wrapper to loop over all the 2.3m data and apply run_daophot.sh"""


import os
import glob
import subprocess as sp
import sys
from pyexcel_ods import get_data
import numpy as np
import math
from functools import partial
import time


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

elif type_choice=='cs':
	info_fpath = working_dir + '/Targets_info.ods'
	if not os.path.exists(info_fpath):
		print 'File cannot be found'
		print info_fpath
		sys.exit()
	


info = get_data(info_fpath)


fintab = []


start = time.time()
#loop over the sheets (named after the standard stars/PN)
for sheetname in info:
	
	if sheetname=='Overview':
		continue
		
		
	if sheetname=='fake_wd1056': continue	
		
		

	print sheetname	
	star_data = info[sheetname]
	colnames = [val.encode('utf-8') for val in star_data[0] ]
	data = [line for line in star_data[1:] if len(line)==8] #lines with all entries. Also removes stars with soemthing in line 9 (ie. manual removal)
	
	
	#make a recarray of each sheet
	data_arr = make_recarray(data, colnames)
	
	#loop over each entry in the sheet. ie. loop of info of each frame
	for line in data_arr:


		if line['Night']!='Night1': continue
		
		


		#clear daophot junk
		print 'Calling tidy_daophot.py'
		sp.call(["tidy_daophot.py"])
	
	
		print 'Frame', line['Filename']

	
		#use the info to get to the correct directory and file
		night_dir = data_dir + '/' + line['Night']
		frame_num = line['Filename']
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
			print sheetname, line['Night'], frame_num
			print
			sys.exit()
			
			
			
			
			
		#check the path with the reference star coordinates exists	
		coords_path = standard_dir + '/pick_stars.tab'
		ref_tab_flag = False
		if not os.path.exists(coords_path):
			print "Coordinates path doesn't exist"
			print coords_path
			ref_tab_flag = True
			continue 	
		if ref_tab_flag==True: continue
			
			
		#if float(line['fwhm'])>5.0:
		#	continue	
			
			
			
		#IF FWHM>5.0, need to change the value as daophot can't handle it
		#if float(line['fwhm'])>5.0:
		#	fwhm = '5.0'
		#else:
		#	fwhm = line['fwhm']
		
		
		fwhm = line['fwhm']
					
			
		#call run_daophot.sh: pass the frame number, fwhm measured in iraf and filepath	
		print 'Running daophot'
		sp.call(["run_daophot.sh", frame_num, fwhm, img_fpath])

		daophot_flag=False


		#check the length of the PSF file, and if its >1, run allstar
		psf_tab = []
		psf_fpath = os.getcwd()+'/'+frame_num+'.psf'
		with open(psf_fpath, 'r') as f:
			psf_tab = [l.strip().split() for l in f]


		if len(psf_tab)==0:
			print 'PSF failed'
			print sheetname, line['Night'], line['Filename']
			newline = [sheetname, line['Night'], line['Filename'], line['Filter'], 'PSF failed']
			fintab.append(newline)
			
			"""
			#try re-running daophot, changing the .coo file so only the CS is in there
			coo_fpath = os.getcwd()+'/'+frame_num+'.coo'
			with open(coo_fpath, 'r') as f:
				coo_tab = [l.strip().split() for l in f]
			coo_intro = coo_tab[:3]
			coo_tab = coo_tab[3:] #get rid of colnames
				
			coo_coords = [ [ind, float(l[1]), float(l[2])] for ind,l in enumerate(coo_tab)]
			coords = [ float(line['x_coord']), float(line['y_coord']) ]
		
			#function to calculate distance between coords
			dist=lambda true, test: ( abs( true[0]-test[1] )**2 + abs( true[1]-test[2] )**2 )
			nearest = min(coo_coords, key=partial(dist, coords) ) 
			coo_line = coo_tab[ nearest[0] ]
					
			#rewrite the .coo file
			with open(coo_fpath, 'w+') as f:
				for entry in coo_intro:
					for val in entry:
						f.write(val+'	')
					f.write('\n')
				for val in coo_line:
					f.write(str(val)+'	')
				f.write('\n')
			print 'New .coo file written'
			raw_input('')					
			
			
			print 'Re-running daophot'
			sp.call(["rerun_daophot.sh", frame_num, fwhm, img_fpath])
			
			
			
			#check the length of the PSF file, and if its >1, run allstar
			psf_tab = []
			psf_fpath = os.getcwd()+'/'+frame_num+'.psf'
			with open(psf_fpath, 'r') as f:
				psf_tab = [l.strip().split() for l in f]
			"""
				
			if len(psf_tab)==0: daophot_flag = True
			
			continue
		if daophot_flag==True: continue


			
	
		os.remove(os.getcwd()+'/'+frame_num+'.lst')
		sp.call(["call_allstar.sh", frame_num, img_fpath])


		
		#if an als file isn't made, the allstar must have failed
		als_fpath = os.getcwd()+'/'+frame_num+'.als'
		if not os.path.exists(als_fpath):
			daophot_flag = True
			newline = [sheetname, line['Night'], line['Filename'], line['Filter'], 'Allstar failed']
			fintab.append(newline)
			continue
		if daophot_flag==True: continue
		print 'Checked als file'
		
		
		
		
		#add pixel coords to info file, open the als file and read the magnitude
		with open(als_fpath) as f:
			mags_list = [l.split() for l in f if len(l.split())!=0]
		mags_list = mags_list[1:] #skip the first line with column names
	
			
		#als_colnames = ['ID', 'X', 'Y', 'mag', 'mag_err', 'sky', 'iterations', 'chi', 'sharp']
		


		#find the coordinates in the als file cloest to that of the standard star	
		als_coords = [ [ind, float(l[1]), float(l[2])] for ind,l in enumerate(mags_list)]
		coords = [ float(line['x_coord']), float(line['y_coord']) ]
		
		#function to calculate distance between coords
		dist=lambda true, test: ( abs( true[0]-test[1] )**2 + abs( true[1]-test[2] )**2 )
		nearest = min(als_coords, key=partial(dist, coords) ) 

		
		#check the identified star is within a reasonable distance of where it should be
		difference_val = 10.
		match_star_flag = False
		if abs(nearest[1]-float(line['x_coord']))>difference_val or abs(nearest[2]-float(line['y_coord']))>difference_val:
			print 'Standard star coordinates could not be found'
			print als_fpath
			newline = [ sheetname, line['Night'], line['Filename'], line['Filter'], 'Coords not found']
			fintab.append(newline)
			match_star_flag = True
			continue
		if match_star_flag==True: continue


		#find the star with the nearest coordinates in the mags_list to get magnitudes
		standard_line = mags_list[ nearest[0]]


		#calculate the "normalised" magnitude, ie. take out the exposure time	
		raw_counts = 10**( float(standard_line[3])/2.5 )
		norm_counts = raw_counts / float(line['Exp_time'])
		norm_mag = 2.5*math.log10(norm_counts)

		
		err_count = 10**( ( float(standard_line[3])+float(standard_line[4]) ) / 2.5 ) / float(line['Exp_time'])
		norm_err = 2.5*math.log10(err_count) - norm_mag
		
		
		
		#pn, night, frame number, exposure time, daophot mag, daophot err, normal mag, normal err
		newline = [sheetname, line['Night'], line['Filename'], line['Filter'],  line['Exp_time'], line['Airmass'], standard_line[3], standard_line[4], norm_mag, norm_err]
		fintab.append(newline)
		print newline
		print 'Line added'
		print
		print
		print
		print
		print
		
		


#save the table to file
print 'Saving'
if type_choice=='s': fin_fpath = working_dir + '/wrapper_standard_mags.tab'
elif type_choice=='cs': fin_fpath = working_dir + '/wrapper_cs_mags.tab'
	
with open(fin_fpath, 'w+') as f:
	f.write('star	night	frame_number	filter	exposure_time	airmass	daophot_mag	daophot_err	mag	mag_err \n')
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


print 'star	night	frame_number	filter	exposure_time	airmass	daophot_mag	daophot_err	mag	mag_err'
for l in fintab:
	print l
print

end = time.time()
print 'Time elapsed: ', end-start	
	
	
	
	
	
	
	
	
	
