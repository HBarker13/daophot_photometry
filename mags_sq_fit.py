#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Use a linear least-squares method to calculate observed magntiudes from instrumental magnitudes"""
#NB. I've erronously used pn as a variable name instead of star. This script
#processes data from standard STARS, not PN

 
from pyexcel_ods import get_data
import json
import numpy as np
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit
import os
import sys



#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):
	dtype_list = ['|S20' if item!='Airmass' and item!='mag' and item!='mag_err' else '>f4' for item in title_list]
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]

	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]	
		data_array.append(col)
		
	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array	




#read in ods file with standard star magnitudes and airmasses
working_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU'
standard_fpath = working_dir + '/Standards.ods'

#working_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/Gaia_reduction'
#standard_fpath = working_dir + '/Standards_gaia.ods'

if not os.path.exists(standard_fpath):
	print 'File cannot be found'
	print standard_fpath
	sys.exit()


all_data = get_data(standard_fpath)


tab = dict()
print 'Reading in data'
for sheetname in all_data:
	
	if sheetname=='Overview':
		continue	
		

	pn_data = all_data[sheetname]
	

	#skip line[0], contains unnecessary titles
	colnames = [val.encode('utf-8') for val in pn_data[1] ]

	
	pn = pn_data[2:]
	
	
	#skip any lines shorter than 12 (lines with a remove flag are 13)
	cut_table = []
	for line in pn:
		if len(line)==len(colnames):
			cut_table.append(line)
			
	
	pn_arr = make_recarray(cut_table, colnames)
	
	
	#remove rows witha  remove flag
	#pn_arr = pn_arr[ pn_arr['remove']!='x']
	#pn_arr = pn_arr[ pn_arr['remove']!='?']
	
		
	tab[sheetname]=pn_arr
	

		
#pn name: [U, B, V, I ]
true_colours = { 'MCT2019': [12.185, 13.397, 13.685, 13.946], '111-1925':[13.045, 12.783, 12.387, 11.904], 'wd1056':[12.775, 13.86, 14.047, 14.350], 'G93-48':[11.942, 12.732, 12.743, 12.938], 'LSE44':[11.042, 12.194, 12.459, 12.604], '121-968':[9.16, 10.071, 10.256, 10.429]  }	

 
#need to solve the following equations for each filter on each night
#
# B = O_B + b + C_B*(B-V) - KB*Z   
#
# O = instrumental offsets
# C = colour terms
# K = extinction coefficients
# Z = airmass
# B-V is the colour on the standard system


filternames = ['U', 'B', 'V', 'I']
nights = ['Night1', 'Night2', 'Night3']
for night in nights:
	
	for filtername in filternames:	
	
		filter_data = []
		for pn in tab:
		
			
			#get the 'true' observed mags from the Overview sheet
			if filtername =='U':
				true_mag = true_colours[pn][0]
				colour = true_colours[pn][0] - true_colours[pn][1] #U-B
				
			elif filtername=='B':
				true_mag = true_colours[pn][1]	
				colour = true_colours[pn][1] - true_colours[pn][2] #B-V
				
			elif filtername=='V':
				true_mag = true_colours[pn][2]
				colour = true_colours[pn][1] - true_colours[pn][2] #B-V
				
			elif filtername=='I':
				true_mag = true_colours[pn][3]
				colour = true_colours[pn][2] - true_colours[pn][3] #V-I	
			
			
	
			#choose the correct night and filter
			pn_mags = tab[pn]			
			cut_mags = pn_mags[ pn_mags['Night']==night ] 
			cut_mags = cut_mags[ cut_mags['Filter']==filtername ]

			
			if len(cut_mags)==0:
				continue
			
			for line in cut_mags:
				print pn, line['Filter'], line['Airmass'], line['mag']
				filter_data.append( [line['Airmass'], line['mag'], line['mag_err'], true_mag, colour] )
			print	
				
				 
		if len(filter_data)==0:
			continue		

		
		#Calculate the transformation coefficients O, C, K using a least-squares method
		def mag_transform(x, O, C, K):
			#true_mag = O + instrumental_mag + C*true_colour - K*airmass
			return O + x[1] + (C*x[4]) - (K*x[0])
			

		true_mags = np.array([line[3] for line in filter_data])
		popt, pcov = curve_fit(mag_transform, np.transpose( filter_data ), true_mags)
		
		O = popt[0]
		C = popt[1]
		K = popt[2]
		
		#Errors
		perr = np.sqrt(np.diag(pcov))
		Oerr = perr[0]
		Cerr = perr[1]
		Kerr = perr[2]

		
		print night, ' Filter: ', filtername
		print 'O: ', O, ' +/- ', Oerr
		print 'K: ', K, ' +/- ', Kerr
		print 'C :', C, ' +/- ', Cerr



		#plot figures like paper 1, fig 1

		#Y1 = true_mag - O - instrumental_mag -C*true_colour
		#Y2 = true_mag - O - instrumental_mag + K*airmass

		Y1 = [line[3] - O - line[1] - C*line[4] for line in filter_data]		
		airmass = [line[0] for line in filter_data]

				
		
		#best fit line
		# Y1 = -K*airmass
		x = np.linspace(min(airmass), max(airmass), 100)
		y = [ -K*val for val in x  ]
 	
		plt.figure()
		plt.plot( airmass, Y1, 'o')
		plt.plot(x, y, 'k--')
		plt.ylabel('Y1')
		plt.xlabel('Airmass')
		plt.title(night+','+filtername)
		plt.show()


 
 
 
 
 
 
 
 
 
 
 
 
 


