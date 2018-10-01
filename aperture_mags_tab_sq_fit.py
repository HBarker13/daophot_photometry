#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Use a linear least-squares method to calculate observed magntiudes from instrumental magnitudes"""
"""uses the table created by aperture_photometry.py as input"""
 
import numpy as np
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit
import os
import sys



#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):
	dtype_list = ['>f4' for item in title_list]
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]

	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]	
		data_array.append(col)
		
	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array	




#read in ods file with standard star magnitudes and airmasses
working_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU'
standard_fpath = working_dir + '/aper_standard_mags.tab'


if not os.path.exists(standard_fpath):
	print 'File cannot be found'
	print standard_fpath
	sys.exit()


print 'Reading in data'
in_tab = np.genfromtxt(standard_fpath, dtype=None, comments='#', names=True)

#remove lines where photometry failed
#in_tab = in_tab[ in_tab['mag']!='PSF star selection failed']

#skip wd1056 - its very close to another star
#in_tab = in_tab[ in_tab['star']!='wd1056']


#put data into a dictionary, with target names as keys
print 'Sorting'
tab = dict()
target_names = set(in_tab['star'])

for target in target_names:
	target_arr = in_tab[ in_tab['star']==target]
	tab[target] = target_arr
print
	

		
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

savetab=[]
print 'Calculating transformation coefficients'
filternames = ['U', 'B', 'V', 'I']
nights = ['Night1', 'Night2', 'Night3']
for night in nights:
	print night
	
	
	counter = 0
	for filtername in filternames:	
		print filtername
	
		filter_data = []
		for target in tab:
			print target
			
				

			
			#get the 'true' observed mags from the Overview sheet
			if filtername =='U':
				true_mag = true_colours[target][0]
				colour = true_colours[target][0] - true_colours[target][1] #U-B
				
			elif filtername=='B':
				true_mag = true_colours[target][1]	
				colour = true_colours[target][1] - true_colours[target][2] #B-V
				
			elif filtername=='V':
				true_mag = true_colours[target][2]
				colour = true_colours[target][1] - true_colours[target][2] #B-V
				
			elif filtername=='I':
				true_mag = true_colours[target][3]
				colour = true_colours[target][2] - true_colours[target][3] #V-I	
			
			
	
			#choose the correct night and filter
			target_mags = tab[target]
	
			cut_mags = target_mags[ target_mags['night']==night ] 
			cut_mags = cut_mags[ cut_mags['filter']==filtername ]

			
			if len(cut_mags)==0:
				continue
			
			for line in cut_mags:
			
				#remove bad lines
				if night=='Night1' and filtername=='V' and target=='wd1056' and line['frame_number']==88: continue
				if night=='Night1' and filtername=='I' and target=='wd1056' and line['frame_number']==93: continue
				if night=='Night1' and filtername=='I' and target=='wd1056' and line['frame_number']==94: continue
				if night=='Night1' and filtername=='I' and target=='wd1056' and line['frame_number']==95: continue
				if night=='Night2' and filtername=='V' and target=='wd1056' and line['frame_number']==67: continue
				if night=='Night2' and filtername=='V' and target=='wd1056' and line['frame_number']==68: continue
				if night=='Night2' and filtername=='V' and target=='wd1056' and line['frame_number']==69: continue
			
			
				counter+=1
				print counter, line['frame_number'],  float(line['airmass']), float(line['mag']), float(line['mag_err'])
				filter_data.append( [ float(line['airmass']), float(line['mag']), float(line['mag_err']), float(true_mag), float(colour), counter] )
			print	
				
				 
		if len(filter_data)==0:
			print 'No stars'
			print
			continue
			
		
		num_stars = len(filter_data)	
		print 'Number of stars:', num_stars	
		print 
		print	
		
		#plt.figure()
		#airmasses = [ l[0] for l in filter_data]
		#mags = [ l[1]-l[3] for l in filter_data]
		#errs = [ l[2] for l in filter_data]
		#plt.errorbar(airmasses, mags, yerr=errs, fmt='o' )
		#plt.show()
		
		filter_tab = make_recarray(filter_data, ['airmass', 'mag', 'mag_err', 'true_mag', 'colour', 'counter'])
		
		#Calculate the transformation coefficients O, C, K using a least-squares method
		def mag_transform(x, O, C, K):
			return O + x['mag'] + (C*x['colour']) - (K*x['airmass'])	

		true_mags = filter_tab['true_mag']
		popt, pcov = curve_fit(mag_transform, filter_tab, true_mags, sigma=filter_tab['mag_err'])

		
		O = popt[0]
		C = popt[1]
		K = popt[2]
		
		#Errors
		perr = np.sqrt(np.diag(pcov))
		Oerr = perr[0]
		Cerr = perr[1]
		Kerr = perr[2]

		
		#print night, ' Filter: ', filtername
		print 'O: ', O, ' +/- ', Oerr
		print 'K: ', K, ' +/- ', Kerr
		print 'C :', C, ' +/- ', Cerr
		




		#plot figures like paper 1, fig 1

		#Y1 = true_mag - O - instrumental_mag -C*colour
		#Y2 = true_mag - O - instrumental_mag + K*airmass

		Y1 = [line['true_mag'] - O - line['mag'] - C*line['colour'] for line in filter_tab]		
		airmass = filter_tab['airmass']
		inds = filter_tab['counter']

				
		
		#best fit line
		# Y1 = -K*airmass
		x = np.linspace(min(airmass), max(airmass), 100)
		y = [ -K*val for val in x  ]
 	
		plt.figure()
		plt.plot( airmass, Y1, 'o')
		plt.plot(x, y, 'k--')
		
		for i,txt in enumerate(inds):
			plt.annotate( txt, (airmass[i], Y1[i]) )
		
		
		plt.ylabel('Y1')
		plt.xlabel('Airmass')
		plt.title(night+','+filtername)
		plt.show()
		#print

		savetab.append(night+'	'+filtername+'	'+str(O)+'	'+str(Oerr)+'	'+str(K)+'	'+str(Kerr)+'	'+str(C)+'	'+str(Cerr)+'	'+str(num_stars)+'	\n')
		
		

		
#write to file
print 'Saving'
savepath = working_dir+'/calib_coeffs.tab'
with open(savepath, 'w+') as f:
	#colnames = Night Filter O Oerr K Kerr C Cerr
	f.write('Night	Filter	O	O_err	K	K_err	C	C_err	num_stars	\n')
	for line in savetab:
		f.write(line)
		
 
 
 
 
 
 
 
 
 
 
 
 
 


