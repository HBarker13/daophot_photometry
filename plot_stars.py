#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""read in the .ads file created by daophot/allstar, and plot the coordinates of detected
stars onto the original image"""


import os
import glob
from astropy.io import fits
import numpy as np

print
#assume working from the current directory
srt_fpath = os.getcwd()+'/'+raw_input('.srt filename: ')
if srt_fpath[-4: ]!= '.srt':
	srt_fpath+='.srt'


#read in the table, skipping header material
with open(srt_fpath, 'r') as f:
	tab = [line.strip().split() for i,line in enumerate(f) if i>2]
	
#round coordinates to inetegers	
x_coords = [int( float(line[1]) ) for line in tab]
y_coords = [int( float(line[2]) ) for line in tab]


#use the original image for the pixel grid template
img_fpath = os.getcwd()+'/'+raw_input('img name: ')
if img_fpath[-5:]!='.fits':
	img_fpath +='.fits'
print img_fpath
	
openimg = fits.open(img_fpath)
hdr = openimg[0].header
data = openimg[0].data
openimg.close()
canvas = np.zeros(data.shape)
	


#make the canvas pixel value 1.0 wherever there is a star
for pair in zip(x_coords, y_coords):
	canvas[pair[1]][pair[0]]=1
	
	

#save the file
savepath = os.getcwd()+'/star_points.fits'
hdu = fits.PrimaryHDU(data=canvas, header=hdr)
hdu.writeto(savepath, clobber=True)
print 'Saved', savepath
print












