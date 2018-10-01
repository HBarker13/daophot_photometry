#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Use the calibration factors made by mags_tab_sq_fits.py to calculate the true CS magntidues from the instrumental ones."""


import numpy as np
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit
import os
import sys
import math


#read in instrumental mags, and calibration values
working_dir = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU'
calib_path = working_dir + '/calib_coeffs.tab'
cs_mags_path = working_dir + '/aper_cs_mags.tab'


#read in calibration values
#Night	Filter	O	O_err	K	K_err	C	C_err	
calib_tab = np.genfromtxt(calib_path, dtype=None, comments='#', names=True)	

#read in instrumental magnitudes
raw_mags = np.genfromtxt(cs_mags_path, dtype=None, comments='#', names=True)	

			
#get object names
names = set(raw_mags['star'])	

nights = ['Night1', 'Night2', 'Night3']	
filtername = ['U', 'B', 'V', 'I']
	

#table to save results to
savetab = []
savetab.append(['name',	'night', 'U', 'B', 'V', 'I'])

	
#loop through objects
for name in names:
	if name!='K1-22': continue
	print name
	
	for night in nights:
		#if night!='Night2':continue
		print night
		
		#get the calibration details for the night
		O_U = None
		K_U = None
		C_U = None
		
		O_B = None
		O_B = None
		O_B = None
		
		O_V = None
		O_V = None
		O_V = None
		
		O_I = None
		O_I = None
		O_I = None		

		#find the correct calibration values
		calib = calib_tab[calib_tab['Night']==night]		
		for l in calib:
				if l['Filter']=='U':
					O_U = float( l['O'] )
					K_U = float( l['K'] )
					C_U = float( l['C'] )
					
					O_U_err = float( l['O_err'] )
					K_U_err = float( l['K_err'] )
					C_U_err = float( l['C_err'] )
					
				elif l['Filter']=='B':
					O_B = float( l['O'] )
					K_B = float( l['K'] )
					C_B = float( l['C'] )	
					
					O_B_err = float( l['O_err'] )
					K_B_err = float( l['K_err'] )
					C_B_err = float( l['C_err'] )
		
				elif l['Filter']=='V':
					O_V = float( l['O'] )
					K_V = float( l['K'] )
					C_V = float( l['C'] )
										
					O_V_err = float( l['O_err'] )
					K_V_err = float( l['K_err'] )
					C_V_err = float( l['C_err'] )
					
				elif l['Filter']=='I':
					O_I = float( l['O'] )
					K_I = float( l['K'] )
					C_I = float( l['C'] )
					
					O_I_err = float( l['O_err'] )
					K_I_err = float( l['K_err'] )
					C_I_err = float( l['C_err'] )	
		
		
		
		
		
		
		
		
		
		
		
		
		#instrumental mags of the target CS
		inst_mags = raw_mags[ raw_mags['star']==name]
		inst_mags = inst_mags[ inst_mags['night']==night]



		U=None
		B=None
		V=None
		I=None
		

		
		#averaging
		_mags = inst_mags[ inst_mags['filter']=='U']
		u = np.mean(_mags['mag'])
		if len(_mags)!=0:
			u_err = math.sqrt( sum( [ l['mag_err']*l['mag_err'] for l in _mags ]) )/len(_mags)
		else: u_err = None
		U_airmass = np.mean(_mags['airmass'])
		
		_mags = inst_mags[ inst_mags['filter']=='B']
		b = np.mean(_mags['mag'])
		if len(_mags)!=0:
			b_err = math.sqrt( sum( [ l['mag_err']*l['mag_err'] for l in _mags ]) )/len(_mags)
		else: b_err = None
		B_airmass = np.mean(_mags['airmass'])
		
		_mags = inst_mags[ inst_mags['filter']=='V']
		v = np.mean(_mags['mag'])
		if len(_mags)!=0:
			v_err = math.sqrt( sum( [ l['mag_err']*l['mag_err'] for l in _mags ]) )/len(_mags)
		else: v_err = None
		V_airmass = np.mean(_mags['airmass'])
		
		_mags = inst_mags[ inst_mags['filter']=='I']
		i = np.mean(_mags['mag'])
		if len(_mags)!=0:
			i_err = math.sqrt( sum( [ l['mag_err']*l['mag_err'] for l in _mags ]) )/len(_mags)
		else: i_err = None
		I_airmass = np.mean(_mags['airmass'])						
		
		"""
		for line in inst_mags:
			if line['filter']=='U': U = line
			elif line['filter']=='B': B = line
			elif line['filter']=='V': V = line
			elif line['filter']=='I': I = line
			else:
				print 'Error: filtername problem'
				print line
				sys.exit()
				
			
		
		try:
			u = float(U['mag'])
			u_err = float(U['mag_err'])
			U_airmass = float( U['airmass'] )
		except:
			u = None
			U_airmass = None
		
		try:
			b = float(B['mag'])
			b_err = float(B['mag_err'])
			B_airmass = float( B['airmass'] )
		except:
			b = None
			B_airmass = None
		
		try:	
			v = float(V['mag'])
			v_err = float(V['mag_err'])
			V_airmass = float( V['airmass'] )
		except:
			v = None
			V_airmass = None
			
		try:	
			i = float(I['mag'])
			i_err = float(I['mag_err'])
			I_airmass = float( I['airmass'] )
		except:
			i = None
			I_airmass = None
			
		"""	
			

	
				
				
		#first run		
		U = u
		B = b
		V = v
		I = i
		U_err = None
		B_err = None
		V_err = None
		I_err = None	

			
		#loop until the values converge	
		test = 100.	
		counter = 0
		while(test>0.000000001):	
		
			U_prev = U
			B_prev = B
			V_prev = V
			I_prev = I
			
			print U,B,V,I
			raw_input('')
			
			#skip calculations if there isn't a measurments for a target in a particular filter
			if not np.isnan(U) and not np.isnan(B):
				U = ( 1 / (1-C_U) ) * ( O_U + u - (K_U*U_airmass) - (C_U*B ) )			
				U_err_squared = ( O_U_err**2 +  u_err**2 + ( (U-B)**2 * C_U_err**2 ) ) / (1 - (2*C_U**2) ) 
				U_err = math.sqrt(U_err_squared)
				
			if not np.isnan(B) and not np.isnan(V):	
				B = ( 1 / (1-C_B) ) * ( O_B + b - (K_B*B_airmass) - (C_B*V ) )
				B_err_squared = ( O_B_err**2 +  b_err**2 + ( (B-V)**2 * C_B_err**2 ) ) / (1 - (2*C_B**2) ) 
				#+ ( B_airmass_err**2 * Kerr**2 ) + apcor_err**2 
				B_err = math.sqrt(B_err_squared)
				
				
				V = ( 1 / (1+C_V) ) * ( O_V + v - (K_V*V_airmass) + (C_V*B ) ) #not a typo
				V_err_squared = ( O_V_err**2 +  v_err**2 + ( (B-V)**2 * C_V_err**2 ) ) / (1 - (2*C_V**2) ) 
				V_err = math.sqrt(V_err_squared)
				
			if not np.isnan(V) and not np.isnan(I):	
				I = ( 1 / (1+C_I) ) * ( O_I + i - (K_I*I_airmass) + (C_I*V ) )
				I_err_squared = ( O_I_err**2 +  i_err**2 + ( (V-I)**2 * C_I_err**2 ) ) / (1 - (2*C_I**2) ) 
				I_err = math.sqrt(I_err_squared)
			
			
			#calculate error
		
			#Orsola calulated errors:
			#I assume airmass error=0 and apcorr_err = 0
						

		
			#where B-V_err is some complicated thing  
			
			#test = math.sqrt( ( (U-U_prev)**2 + (B-B_prev)**2 + (V-V_prev)**2 + (I-I_prev)**2 ) )
			test = 0
			if U!=None: test+= (U-U_prev)**2
			if B!=None: test+= (B-B_prev)**2
			if V!=None: test+= (V-V_prev)**2
			if I!=None: test+= (I-I_prev)**2
			
			test = math.sqrt(test)
		
			
		
			counter+=1
			if counter>1000:
				print 'Counter >1000'
				print 'Convergence failed:', night, name
				raw_input('')
		

		newline = [ name, night, U, U_err, B, B_err, V, V_err, I, I_err]
		print newline
		savetab.append(newline)
		print
		print
		
		
for line in savetab:
	print line

		
		
			

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		


 
 
 
 
 


