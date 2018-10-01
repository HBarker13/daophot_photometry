#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Calculate average magnitudes of objects in each filter using the results from daophot_wrapper.py  ie. wrapper_mags.tab files"""


import math
import numpy as np


working_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU'
filternames = ['U', 'B', 'V', 'I']
nights = ['Night1', 'Night2', 'Night3']


#choose which set of objects to process
type_choice = None
while type_choice!='s' and type_choice!='t':
	type_choice = raw_input('Averge magnitudes for standard stars (s) or CS targets (t)? ')
	type_chocice = type_choice.lower()
	
	
#read in the mags, files created by daophot_wrapper.py	
if type_choice=='s':
	mags_fpath = working_dir+'/wrapper_standard_mags.tab'
elif type_choice=='t':
	mags_fpath = working_dir+'/wrapper_cs_mags.tab'


with open(mags_fpath, 'r') as f:
	all_mags = [line.strip().split() for line in f]
	

#skip the first line (column names)
#star	night	frame_number	filter	exposure_time	airmass	daophot_mag	daophot_err	mag	mag_err 
all_mags = all_mags[1:]
all_mags = [line for line in all_mags if len(line)==10]


#get the names of all the objects
names = set( [line[0] for line in all_mags] )


fintab = []

#loop through the object names
for name in names:
	#print name
	
	
	#loop through nights
	for night in nights:
		#print night
		night_mags = [line for line in all_mags if line[1]==night]
	
	
		#get all the lines in the table with the same object name
		obj_mags = [line for line in night_mags if line[0]==name]
	
	
	
		#loop through the filters
		for filtername in filternames:
			#print filtername
		
			#get all the lines in obj_mags with the same filter
			filter_lines = [line for line in obj_mags if line[3]==filtername]
			if len(filter_lines)==0: continue


			#calculate average airmass
			airmasses = [ float(line[5]) for line in filter_lines]
			avg_airmass = np.mean( airmasses )


		
			#calculate the average magnitude value
			#need to convert into counts, can't just sum magnitudes. mags have been altered so exposure time=1
			raw_mags = [float(line[8]) for line in filter_lines]
			raw_counts = [ 10**( line/2.5 ) for line in raw_mags]
			avg_counts = sum(raw_counts)/len(raw_counts)
			avg_mag = 2.5 * math.log10(avg_counts)
		

			#calculate the average error
			raw_err = [ [float(line[8]), float(line[9])] for line in filter_lines]

			#convert to counts
			raw_counts_upper = [ (10**( line[1]/2.5 )) + (10**( line[0]/2.5 )) for line in raw_err]
			raw_counts_lower = [ (10**( line[1]/2.5 )) - (10**( line[0]/2.5 )) for line in raw_err]
		
		
			#err = math.sqrt( err1^2 + err2^2 + err3^2 ) / 3
			sq_upper_err = [ line**2 for line in raw_counts_upper]
			upper_counts_err = math.sqrt( sum(sq_upper_err)  /len(sq_upper_err) )
		
			sq_lower_err = [ line**2 for line in raw_counts_lower]
			lower_counts_err = math.sqrt( sum(sq_lower_err)  /len(sq_lower_err) )
		

			upper_err = ( 2.5*math.log10( upper_counts_err ) - avg_mag)
			lower_err = ( 2.5*math.log10( lower_counts_err ) - avg_mag)
		
			#average the two, the values are so close transforming to counts is unnecessary
			avg_err = (upper_err + lower_err)/2
		
		
		
			newline = [ name, night, filtername, avg_airmass, avg_mag, avg_err ]
			fintab.append(newline)
	
	
	
	
#save the new magnitudes
if type_choice=='s':
	savepath = working_dir + '/avg_standard_mags.tab'
elif type_choice=='t':
	savepath = working_dir + '/avg_cs_mags.tab'
	

with open(savepath, 'w+') as f:
	f.write('Star	Night	Filter	Airmass	Mag	Mag_err	\n')
	for line in fintab:
		f.write( str(line[0])+'	'+str(line[1])+'	'+str(line[2])+'	'+str(line[3])+'	'+str(line[4])+'	'+str(line[5])+'\n')
		
print 'Saved', savepath	


	
	
	
	
	
	
	
	
	
	
	
	
	








