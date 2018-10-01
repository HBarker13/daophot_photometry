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
     
     
#Parameters for SEP  
#arcsec_rad=1.5 #aperture radius in arcsec, should be roughly 1.5x the seeing   
#pixscale=0.338 #arcsec per pixel, from Imager website: http://rsaa.anu.edu.au/observatories/instruments/imager
#aper=arcsec_rad/pixscale #in pixels

gain=1.1 #gain in e-/adu

fwhm = 3.3 #pixels
aper = 1.5*fwhm


#Open file
file="test.fits"
hdu=fits.open(file)
header=hdu[0].header
data=hdu[0].data

#Sep commands
data = data.byteswap().newbyteorder()
bkg = sep.Background(data) #calculate the background


#looks fine on my test frame
#print bkg.globalback
#print bkg.globalrms

#subtract the background from the image to see what it looks like
data_sub = data - bkg
"""
sub_path = os.getcwd()+'/skysub.fits'        
skysub = fits.PrimaryHDU( data_sub )
skysub.writeto(sub_path, clobber=True)
"""


#Set a detection threshold for stars to be 1.5x RMS of background
thresh = 1.5 * bkg.globalrms
#Find all objects that are above threshold in a frame which is the background subtracted data
objects = sep.extract(data-bkg.back(), thresh)
print 'Number of objects detected:', len(objects)

"""
#plot the location of the objects on the image
fig, ax = plt.subplots()
im = ax.imshow(data, interpolation='nearest', cmap='gray', vmin= bkg.globalback-bkg.globalrms, vmax=bkg.globalback+bkg.globalrms, origin='lower')
for i in range( len(objects) ):
	e = Ellipse( xy=(objects['x'][i], objects['y'][i]), width=6*objects['a'][i], height=6*objects['b'][i], angle=objects['theta'][i]*180./np.pi)
	e.set_facecolor('none')
	e.set_edgecolor('red')
	ax.add_artist(e)
plt.show()
"""




#Table with list of object coords and fluxes
stars=Table( [objects['x'],objects['y'],objects['flux']], names=['x','y','sepflux'])
#NOTE this can get a little messed up when the background is a bit weird,
#like here where you have sky only in a circular region 
#so better to extract things yourself and just use the object detections



#This calculates aperture fluxes, uncertainties and flags for all objects.
#data_sub = sky subtracted data
#fluxes, flux_err, flag = sep.sum_circle( data_sub, stars['x'], stars['y'], aper, err=bkg.globalrms, gain=gain)


#Can also do sky subtraction using annuli if you don't trust sep to calculate the background
ann_in_pix = 3*fwhm  #pixels
ann_out_pix = 4*fwhm
fluxes, fluxerrs, flag = sep.sum_circle( data, stars['x'], stars['y'], aper, err=bkg.globalrms, gain=gain, bkgann=(ann_in_pix,ann_out_pix) )



#Add fluxes and uncertainties to table
stars.add_column( Column( fluxes, name='flux') )
stars.add_column( Column( fluxerrs, name='flux_err') )



#Calculate instrumental magnitudes
expt = 20.
mags = [ -2.5*math.log10( line / expt) if line>0 else float('nan') for line in fluxes ]


print mags








