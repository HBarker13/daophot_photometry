#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#test script using pyraf to use iraf to do daophot reduction
#Using notes by Geert, https://github.com/barentsen/surveytools/blob/master/surveytools/daophot.py


from pyraf import iraf
from iraf import digiphot

import os
import glob
import shutil
import numpy as np
import math

from astropy.io import fits
from astropy import table
from astropy.table import Table

from functools import partial
from pyexcel_ods import get_data




#get the daophot tasks
iraf.digiphot.daophot()
print 'Iraf loaded'
print

# Allow overwriting files
iraf.set(clobber="yes")






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





def set_params():
	#Set daophot parameters
	#http://iraf.noao.edu/scripts/irafhelp?val=datapars
	
	
	# Ensure we start from the iraf defaults
        for module in ['datapars', 'findpars', 'centerpars', 'fitskypars', 'photpars', 'daopars', 'daofind', 'phot']:
        	iraf.unlearn(module)
	
	

	iraf.datapars.scale = 1.0 #set scale as pixels
	iraf.datapars.epadu = 1.1 #gain in electrons per adu
	iraf.datapars.readnoise = 3.278 #readout noise in electrons

	iraf.daopars.function = "gauss" #functional form of the analytic component of the psf model  #USING GAUSSIAN FOR NOW
	iraf.daopars.varorder = 0 #order of variablilty of the psf model. 0 = constant over the image. I don't hav enough references to make a psf that varies
	iraf.daopars.nclean = 5 #number of additional iterations of PSF to compute the PSF look-up tables
	iraf.daopars.saturated = "no" # use saturated stars to computer the psf wings?
	iraf.daopars.recenter = "yes" #compute new positions as well as magnitudes for stars
	iraf.daopars.fitsky = "yes" # compute new sky values for stars in the input list
   
	iraf.centerpars.calgorithm = "none" #can be none if initial potions are acurate ie. if they are calculated by daofind
	iraf.photpars.zmag = 20.0 #zero-point offet for the magnitude scale
	
	
	#now the data-dependant parameters
	#iraf.centerpars.cbox = 3*fwhm #width of the subraster used for object centering in pixels. Recommended values are 2.5 - 4 * fwhm
	
		
	iraf.datapars.sigma = skysigma #standard deviation of the sky in pixels. Should be representative of sky noise. Default=5.75
	
	
	
	
	
	
	iraf.datapars.datamin = -100 #Iain says this should be fine, and means I can use the U and B frames
	#iraf.datapars.datamin = 3*iraf.datapars.sigma #minimum good pixel value (in ADU) Default=3
	iraf.datapars.datamax = 36000 #maximum good pixel value (ADU)
	

	
	
	#-------------DAOPHOT PARAMS---------------------#
	
	
	iraf.datapars.itime = exptime #exposure time in arbitary units
	iraf.datapars.fwhmpsf = fwhm #fwhm in pixels

	iraf.daopars.fitrad = 1.2*fwhm  #fitting radius in pixels. Only pixels within the fitting radius of the center of a star will  contribute  to  the fits  computed  by  the  PEAK, NSTAR and ALLSTAR tasks. For most frames, fitrad should = fwhm, but can be larger for uncrowded frames

	iraf.daopars.psfrad = 4*fwhm  #radius in pixels, within which the psf must be defined

	
	iraf.daopars.sannulus = 4*fwhm #inner sky radius
	iraf.daopars.wsannulus = 1*fwhm #outer sky radius




	#SKY FITTING
	#According to https://www.google.co.uk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0ahUKEwjsl9Pw2JrVAhXDAcAKHbL2DugQFggoMAA&url=ftp%3A%2F%2Firaf.noao.edu%2Fftp%2Fdocs%2Fdaorefman.ps.Z&usg=AFQjCNEm7fBwR19te02kqUMeMMjoBvOp2g
	#use fitskypars
	#if sky contamination is due to nebulosity, salogri = "median", "centroid" or "crosscor"
	iraf.fitskypars.salgorithm = "centroid" # Compute  the  intensity-weighted mean or centroid of the sky pixel histogram.  
	iraf.fitskypars.annulus = 4*fwhm #inner sky radius in pixels
	iraf.fitskypars.dannulus = 1*fwhm #width of sky fitting region



	#list  of  aperture  radii in units of the  scale parameter or the name of the file containing  the  list  of  aperture  radii. 
	#List  elements may be separated by whitespace or commas. A range syntax of the form ap1:apN:apstep is also supported
	#iraf.photpars.apertures = "3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0"
	iraf.photpars.apertures = fwhm
	
	#only want to find the brightest stars
	iraf.findpars.threshold = 5*skysigma #significance threshold for detection. Default=4


	#defaults used
	#print iraf.findpars.roundlo = -1
	#print iraf.findpars.roundhi = 1
	#print iraf.findpars.sharplo = 0.2
	#print iraf.findpars.sharphi = 1






def daofind():
	print 'FIND'
	"""Identifies and computes approximate centroids for small luminous objects in the field. Read out noise and gain in photons (or electrons) must be defined, and a significance level (in standard deviations) given, above which an object is considered real.  Find computes the actual brightness enhancment of an object in data numbers above the local sky brightness.
	
	Searches the image for local density maxima, which have a fwhm = datapars.fwhmpsf and a peak amplitude greater than findpars.threshold * datapars.sigma

	Image = input image
	Output = output name. Default is .coo

	ratio = 1   : The ratio of the sigma of the Gaussian convolution kernel along the minor axis direction to the sigma along the major axis direction. 
	Ratio defaults to 1.0 in which case the image is convolved with a circular Gaussian. 

	theta = 0 : The position of the major axis of the elliptical Gaussian. Theta is measured counter-clockwise from the x axis. 

	nsigma = 1.5 The semi-major axis of the Gaussian convolution kernel used to computed the density enhancement and mean density images in Gaussian sigma. 
	This semi- major axis is equal to min (2.0, 0.42466 * nsigma * datapars.fwhmpsf / datapars.scale) pixels.
	"""

	iraf.daofind(image=img_fpath, output=daofind_out, fwhmpsf=fwhm, ratio=1, theta=0, verify='no', verbose='no')
	print 'Created', daofind_out
	print






def daophot():
	print 'PHOT'
	"""Obtain sky values and concentric aperture photometry for all objects found by Find.
	output default is .mag
	"""
	iraf.daophot.phot(image=img_fpath, coords=daofind_out, output=phot_out, interactive='no', verify='no', verbose='no')
	print 'Performed aperture photometry'
	print 'Created', phot_out
	print
	
	




def pstselect():
	print 'PICK'
	"""Choose candidates for psf fitting
	pstfile = output psf star list file (default .pst)
	maxnpsf = maximum number of psf stars
	
	
	Reads in the photometry file and sorts the list in order of increasing magnitude, after rejecting any stars that have INDEF valued magnitudes, or which lie less than fitrad / scale pixels from the edge of the image . From this list the brightest maxnpsf stars which have no brighter neighbor stars within (psfrad + fitrad ) / scale + 1 pixels are selected as candidate psf stars.
	
	"""
	iraf.daophot.pstselect(image=img_fpath, photfile=phot_out, pstfile=pst_out, maxnpsf=100, verify='no', verbose='yes')

	#check the length of the file
	psf_flag = check_psflist(pst_out)
	
	
	print 'Chosen PSF stars'
	print 'Created', pst_out
	print
	return psf_flag






def psf():
	print 'PSF'
	"""Builds the psf for the image using stars from photfile. 
	
	 Suitable PSF stars are normally selected interactively using the image display and image cursor and matched with the stars in photfile using the cursor position and a tolerance specified by the matchrad parameter in the DAOPARS task. A star must be in the photometry file before it can be used as a PSF star. If a match is found, PSF checks that the candidate star is not too close to the edge of the image and that it contains no bad pixels as defined by datamin and datamax in the DATAPARS task. After selection a mesh, contour, or profile plot of the data subraster around the candidate star is displayed in the graphics window, PSF enters graphics cursor command mode and the user is given the option to accept or reject the star. If the user accepts the star it is added to the PSF star list. Commands in the graphics cursor menu permit the user to manipulate the floor and ceiling levels of the contour plot and the viewing angles for the mesh plot interactively.

Users who know which stars they wish to use as PSF stars ahead of time or who are without access to an image display can also select PSF stars by id number, after which mesh, contour, or radial profile plots will be displayed in the graphics window in the usual way.

If the user does not wish to see any plots of the PSF stars or interact with the fitting process, the image cursor may be redirected to a text file containing cursor commands icommands which specify the PSF stars to be used in the fit. If plotfile is defined contour, mesh, or profile plots of the selected psf stars can be saved in a metacode plot file for later examination
	
	
	
	
	photfile = file with aperture photometry
	pstfile = input psf star list from pstselect (.pst)
	
	psfimage = file to save psf image to, '.fits' is added automatically
	opsfile = output psf star list (default .pst)
	groupfil = output psf star group file (defaul .psg)
	"""
	
	#can't overwrite the psf file, so delete it if it already exists
	if os.path.exists(psfimg_out+'.fits'):
		os.remove( psfimg_out+'.fits')
	

	iraf.daophot.psf(image=img_fpath, photfile=phot_out, pstfile= pst_out , psfimage=psfimg_out, groupfil=psg_out, opstfile=psf_pst_out, verify='no', verbose='no', interactive='no')
	

	
	print 'Performed psf fitting'
	print 'PSF image:', psfimg_out
	print 'PSF fitting stars:', psf_pst_out
	print 

	
	
	
	
def allstar():
	print 'Allstar'
	"""gets psf photometry"""
	
	#need to remove fits files if they exist, they can't be clobbered
	if os.path.exists(subimage+'.fits'):
		os.remove(subimage+'.fits')
	
	
	iraf.daophot.allstar(image=img_fpath, photfile=phot_out, psfimage=psfimg_out, allstarfile=allstar_out, rejfile=None, subimage=subimage, verify='no', verbose='no')
	print 'PSF photometry calculated'
	print 'Photometry file:', allstar_out
	print 'Star subtracted image:', subimage+'.fits'
	print	
	
	
	
	
	
	
#check the psf reference stars list to check its length	
def check_psflist(pstselect_output_path):

	with open(pstselect_output_path, 'r') as f:
		tbl = [ line for line in f if line[0]!='#']
		
	if len(tbl)>1:
		psf_flag = True
	else:
		psf_flag = False
	
	return psf_flag
		


	
	
	
#copied from Geert's github code	
#Seems to do a good job getting rid of junk detections	
def pstselect_prune(pstselect_output_path, new_path):
	print 'PRUNING'
	print

	#rename the old file
	shutil.copyfile(pstselect_output_path, pstselect_output_path+'.old')
	tbl = Table.read(pstselect_output_path, format='daophot')
	

	# We prune objects using sigma-clipping on the sky estimate.
	# We will try increasing values of sigma, until less 40% of the stars
	# are rejected.
	from astropy.stats import sigma_clip
	for sigma in [2., 3., 4., 5., 10.]:
		bad_mask = sigma_clip(tbl['MSKY'].data, sigma=sigma, iters=None).mask
	        if bad_mask.sum() < 0.4*len(tbl):  # stop if <40% rejected
	            break
         

	# Now write the new file without the pruned objects to disk
	fh = open(new_path, 'w+')
	fh.write("#N ID    XCENTER   YCENTER   MAG         MSKY\n"
	         "#U ##    pixels    pixels    magnitudes  counts\n"
	        "#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-15.7g\n")
	for row in tbl[~bad_mask]:
	        fh.write('{0:<9d}{1:<10.3f}{2:<10.3f}{3:<12.3f}{4:<15.7g}\n'.format(
	                  row['ID'], row['XCENTER'], row['YCENTER'],
	                  row['MAG'], row['MSKY']))
	fh.close()	
	
	
	#check the length of the file
	psf_flag = check_psflist(pst_out)
	
	return psf_flag
	
	
	
	
	

#remove any stars near the edges of the frame
def pstselect_remove_outliers(pstselect_output_path, new_path):
	print 'REMOVING OUTLIERS'
	print 

	#rename the old file
	shutil.copyfile(pstselect_output_path, pstselect_output_path+'.old')
	tbl = Table.read(pstselect_output_path, format='daophot')
	
	#get the coordinates of the image to find the centre of the frame
	#frame size should be ~2008, 2008
	template = fits.getdata(img_fpath)
	x_size, y_size = template.shape
	x_mid = int(x_size/2.)
	y_mid = int(y_size/2.)
	
	#remove stars >350 pixels from the centre of the frame
	fh = open(new_path, 'w+')
	fh.write("#N ID    XCENTER   YCENTER   MAG         MSKY\n"
	         "#U ##    pixels    pixels    magnitudes  counts\n"
	        "#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-15.7g\n")
	for row in tbl:
		if math.sqrt( abs( row['XCENTER']-x_mid )**2 + abs( row['YCENTER']-y_mid )**2  ) <350.:
			        fh.write('{0:<9d}{1:<10.3f}{2:<10.3f}{3:<12.3f}{4:<15.7g}\n'.format(
	                  row['ID'], row['XCENTER'], row['YCENTER'],
	                  row['MAG'], row['MSKY']))
	fh.close()
	
	#check the length of the file
	psf_flag = check_psflist(pst_out)
	
	return psf_flag	
		




def pstselect_prechosen(pstselect_output_path):
	#rename the old file
	shutil.copyfile(pstselect_output_path, pstselect_output_path+'.old')

	
	#read in the list of stars I've pre-chosen
	star_dir, _ = img_fpath.rsplit('/', 1)
	pick_stars = star_dir + '/pick_stars_'+ star_info['Filter']+'.tab'
	if not os.path.exists( pick_stars ):
		print 'pick_stars table not found'
		print pick_stars
		raw_input('Paused')
		
	with open(pick_stars, 'r') as f:
		pst_coords = [line.strip().split() for line in f]
	

	#read in the lst file
	with open(pstselect_output_path, 'r') as f:
		tab = Table.read(pstselect_output_path, format='daophot')

		
	#find stars in the lst file that are within 10 pixels
	#of the ones I've pre-chosen
	counter = 0
	newtab = []
	for star in tab:		
		for coords in pst_coords:
			if abs( star['XCENTER']-float(coords[0]) )<15.0 and abs( star['YCENTER']-float(coords[1]) )<15.0:
				newtab.append(star)
				counter+=1
				print 'Reference stars found: ', counter, '/', len(pst_coords)
	


	# Now write the new file without the pruned objects to disk
	fh = open(pstselect_output_path, 'w+')
	fh.write("#N ID    XCENTER   YCENTER   MAG         MSKY\n"
	         "#U ##    pixels    pixels    magnitudes  counts\n"
	        "#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-15.7g\n")
	for row in newtab:
	        fh.write('{0:<9d}{1:<10.3f}{2:<10.3f}{3:<12.3f}{4:<15.7g}\n'.format(
	                  row['ID'], row['XCENTER'], row['YCENTER'],
	                  row['MAG'], row['MSKY']))
	fh.close()	

	
	#check the length of the file
	psf_flag = check_psflist(pst_out)	
	
	return psf_flag





def plot_daofind():
	#Plot the coordinates of the stars identified by daofind in a
	#fits file
	
	#read in the field created by daofind
	with open(daofind_out, 'r') as f:
		daofind_tab = [line.strip().split() for line in f if line[0]!='#']
		
	# x = column number, y = row number	
	y_coords = [ int( float(line[0]) ) for line in daofind_tab]
	x_coords = [ int( float(line[1]) ) for line in daofind_tab]


	#open the original image to use it as a base
	template = fits.getdata(img_fpath)
	canvas = np.zeros(template.shape)
	
	
	for x,y in zip(x_coords, y_coords):
		#make the filled in squares larger
		#so they're easier to see on a ds9 RGB image
		for xval in range(x-5, x+5):
			for yval in range(y-5, y+5):
				canvas[ xval ][ yval ] = 1 
				
		
	#save as a new fits file
	savepath = root + '_daofind-plot.fits'
	hdu = fits.PrimaryHDU( data=canvas )
	hdu.writeto( savepath, overwrite=True )
	
	print 'Plotted daofind stars: ', savepath
	print
	
	
	
	
	
	
	
def plot_psf():
	#Plot the coordinates of the stars identified by psf in a
	#fits file
	
	
	#read in the field created by daofind
	with open(psf_pst_out, 'r') as f:
		psf_tab = [line.strip().split() for line in f if line[0]!='#']
		
	# x = column number, y = row number	
	y_coords = [ int( float(line[1]) ) for line in psf_tab]
	x_coords = [ int( float(line[2]) ) for line in psf_tab]


	#open the original image to use it as a base
	template = fits.getdata(img_fpath)
	canvas = np.zeros(template.shape)
	
	for x,y in zip(x_coords, y_coords):
		#make the filled in squares larger
		#so they're easier to see on a ds9 RGB image
		for xval in range(x-5, x+5):
			for yval in range(y-5, y+5):
				canvas[ xval ][ yval ] = 1 
		
	#save as a new fits file
	savepath = root + '_psf-plot.fits'
	hdu = fits.PrimaryHDU( data=canvas )
	hdu.writeto( savepath, overwrite=True )
	
	print 'Plotted psf stars: ', savepath
	print	







#get the allstar table line for the target star
def get_target():

	#open the allstar file
	tbl = Table.read(allstar_out, format='daophot')
	als_coords = [ [ind, float(line['XCENTER']), float(line['YCENTER'])] for ind,line in enumerate(tbl) ]
	
	#target coordinates
	coords = [ xcoord, ycoord ]


	#find the star in the allstar file with the nearest coordinates
	#to the target's
	dist=lambda true, test: ( abs( true[0]-test[1] )**2 + abs( true[1]-test[2] )**2 )
	nearest = min(als_coords, key=partial(dist, coords) ) 
	target_line = tbl[ nearest[0] ]
	success_flag = True

	#check the identified star is within a reasonable distance of where it should be
	difference_val = 10.
	match_star_flag = False
	if abs( target_line['XCENTER']-float(xcoord) )>difference_val or abs( target_line['YCENTER']-float( ycoord) )>difference_val:
		print 'Standard star coordinates could not be found'
		print als_fpath
		newline = [ name, night, frame_num, 'Coords not found']
		fintab.append(newline)
		success_flag = False

	return success_flag, target_line







#----------------Call the functions----------------#


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
		
	fin_fpath = working_dir + '/pyraf_standard_mags.tab'	


elif type_choice=='cs':
	info_fpath = working_dir + '/Targets_info.ods'
	if not os.path.exists(info_fpath):
		print 'File cannot be found'
		print info_fpath
		sys.exit()
		
	fin_fpath = working_dir + '/pyraf_cs_mags.tab'



#make a recarray of each sheet
info = get_data(info_fpath)



fintab = []

with open(fin_fpath, 'w+') as f:
	f.write('star	night	frame_number	filter	exposure_time	airmass	fwhm	fi	ps	is	os	mag	mag_err	norm_mag	norm_err \n')




for sheetname in info:

	if sheetname=='Overview': continue		
	if sheetname=='fake_wd1056': continue	
	
	
	
	
	
	#if sheetname!='LSE44': continue

	
	
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



		#define the filenames that will be used by the daophot commands
		root = working_dir + '/pyraf_docs/'+ frame_num
		daofind_out = root + '_daofind-out.coo'
		phot_out = root + '_phot-out.mag'
		pst_out = root + '_pstselect-out.pst'
		psfimg_out = root + '_psfimg'  #daophot adds .fits
		psf_pst_out = root + '_psf-pst-out.psg'
		psg_out = root + '_psg-out.psg'
		allstar_out = root+'_allstar-out.txt'
		subimage = root+'_subimage' #.fits added by allstar


		#set the daophot parameters
		set_params()
		
		
		#run daofind and plot the coordinates of detected objects in an image
		daofind()
		plot_daofind()


		#run daophot and check more than one psf reference star is found
		daophot()
		pst_flag = pstselect()
		if pst_flag==False:
			print 'PSF star selection failed'
			newline = [sheetname, star_info['Night'], frame_num, star_info['Filter'], exptime, star_info['Airmass'], 'pstselect_failed']
			fintab.append(newline)
			raw_input('Paused')
			print
			continue
			
			
	
			
		"""
		#only use psf reference stars I've pre-chosen			
		pst_flag = pstselect_prechosen(pst_out)
		if pst_flag==False:
			print 'PSF star selection failed'
			newline = [sheetname, star_info['Night'], frame_num, star_info['Filter'], exptime, star_info['Airmass'], 'PSF star selection failed'  ]
			fintab.append(newline)
			raw_input('Paused')
			print
			continue
		"""	
		
				
		#remove psf reference stars using sigma-clipping on the sky estimate
		pst_flag = pstselect_prune(pst_out, pst_out)
		if pst_flag==False:
			print 'Prune: No PSF stars left in the reference list.'
			newline = [sheetname, star_info['Night'], frame_num, star_info['Filter'], exptime, star_info['Airmass'], 'PSF star selection failed'  ]
			fintab.append(newline)
			continue
		

		
		
		#remove psf reference stars that aren't in the centre of the frame
		#if there are enough stars
		with open(pst_out, 'r') as f:
			psftbl = [ line for line in f if line[0]!='#']
		if len(psftbl)>10:
		
			pst_flag = pstselect_remove_outliers(pst_out, pst_out)
			if pst_flag==False:
				print 'Outliers: No PSF stars left in the reference list.'
				newline = [sheetname, star_info['Night'], frame_num, star_info['Filter'], exptime, star_info['Airmass'], 'PSF star selection failed'  ]
				fintab.append(newline)
				
				raw_input('')
				#continue
						


		#create the psf for the image and plot the coordinates of the
		#objects used to create the psf
		psf()
		
		#skip if psf fails
		if not os.path.exists(psf_pst_out):
			continue
		
		plot_psf()
			
		
		allstar()
		flag,target = get_target()
	
	
	
		#If the standard star couldn't be found, continue
		if flag==False:
			continue
	
		fi = iraf.daopars.fitrad	
		ps = iraf.daopars.psfrad
		inner_sky = iraf.daopars.sannulus
		outer_sky = iraf.daopars.wsannulus





		#calculate the "normalised" magnitude, ie. take out the exposure time	
		raw_counts = 10**( float( target['MAG'] )/2.5 )
		norm_counts = raw_counts / float( exptime )
		norm_mag = 2.5*math.log10(norm_counts)
	
		err_count = 10**( ( float( target['MAG'] )+float( target['MERR']  ) ) / 2.5 ) / float( exptime )
		norm_err = 2.5*math.log10(err_count) - norm_mag



		newline = [sheetname, star_info['Night'], frame_num, star_info['Filter'], exptime, star_info['Airmass'], fwhm, fi, ps, inner_sky, outer_sky, float(target['MAG']), float(target['MERR']), norm_mag, norm_err  ]
		fintab.append(newline)
		print newline
		
		
		print_line = ''
		for s in newline:
			print_line+=str(s)
			print_line+='	'
		print_line+='\n'
		with open(fin_fpath, 'a') as f:
			f.write(print_line)	
		
		#raw_input('PAUSED')






"""
#save the table to file
print 'Saving'

with open(fin_fpath, 'w+') as f:
	f.write('star	night	frame_number	filter	exposure_time	airmass	fwhm	fi	ps	is	os	mag	mag_err	norm_mag	norm_err \n')
	for line in fintab:
		print_line = ''
		for s in line:
			print_line+=str(s)
			print_line+='	'
		print_line+='\n'
		f.write(print_line)	

"""
print 'File written'
print fin_fpath
print



print
print '---------------Completed---------------'
print











