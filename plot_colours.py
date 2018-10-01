#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Calculate the colours of CS at different reddenings (based on calc_reddened_cs_synthetic_colours.py  ) and plot the observed CS colours to see if they are consistent"""


import os
import numpy as np
import glob
import pysynphot as S
from astropy.io import ascii
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt


os.environ['PYSYN_CDBS']

	
magsystem = 'vegamag'

#read in the 100kK, log(g)=7 txt files made by running prepare_tubingen.py
tubigen_path = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels'
spectrum_fpath = glob.glob(tubigen_path + '/100kK_7.0_solar/*.txt')
spectrum_fpath = spectrum_fpath[0]

wl = []
flx = []
with open(spectrum_fpath, 'r') as f:
	for line in f:
		line = line.split()
		#[Angstrom, eg/s/cm2/cm]	
		wl.append( float(line[0]) )
		flx.append( float(line[1])*1e-8  ) #convert to per angstrom
print 'Spectrum read in'		
		
#create pysynphot spectrum
spectrum = S.ArraySpectrum(np.array(wl), np.array(flx), 'angstrom', 'flam')		

synth_colours = []

#use pysynphot to redden the spectrum 
#redden the spectrum  using Cardelli 1989, RV=3.1
E_BmV = [0.0, 0.4, 0.6, 0.8, 1.0, 1.3, 1.5, 2.0]
for e in E_BmV:

	print 'E(B-V):', e
	red_spectrum = spectrum *  S.Extinction(e, 'mwavg')
	
	#plt.figure()
	#plt.plot(red_spectrum.wave, red_spectrum.flux)
	#plt.xlabel(red_spectrum.waveunits)
	#plt.ylabel(red_spectrum.fluxunits)
	#plt.xlim(3000, 10000)
	#plt.ylim(0, 0.6e10)
	#plt.show()
		
	obs_U = S.Observation(red_spectrum, S.ObsBandpass('U'))
	obs_B = S.Observation(red_spectrum, S.ObsBandpass('B'))
	obs_V = S.Observation(red_spectrum, S.ObsBandpass('V'))
	obs_I = S.Observation(red_spectrum, S.ObsBandpass('I'))


	
	U_min_B = obs_U.effstim(magsystem) - obs_B.effstim(magsystem)
	B_min_V = obs_B.effstim(magsystem) - obs_V.effstim(magsystem)
	V_min_I = obs_V.effstim(magsystem) - obs_I.effstim(magsystem)
	

	#paper1: corrections that were supposed to be applied	
	B_min_V+=0.010
	V_min_I-=0.002
	
	newline = [e, U_min_B, B_min_V, V_min_I]
	synth_colours.append(newline)


print 'Synthetic colours:'
for line in synth_colours:
	print line


#read in observed colours
obs_fpath = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/avg_cs_mags.tab'
with open(obs_fpath, 'r') as f:
	obs = [line.strip().split() for line in f]
	
#skip the first line with column names
colnames = obs[0]
obs = obs[1:]

#calculate observed U-B and V-I colours
obs_colours = []

#get list of object names
obj_names = set( [line[0] for line in obs] )
for name in obj_names:
	
	matches = [line for line in obs if line[0]==name]
	obs_u = [line for line in matches if line[2]=='U']
	obs_b = [line for line in matches if line[2]=='B']
	obs_v = [line for line in matches if line[2]=='V']
	obs_i = [line for line in matches if line[2]=='I']
	
	if len(obs_u)==0: continue
	newline = [name, float(obs_u[0][4])-float(obs_b[0][4]), float(obs_v[0][4])-float(obs_i[0][4]) ]
	
	#if len(obs_b)==0:continue
	#if len(obs_v)==0:continue
	#if len(obs_i)==0:continue
	#newline = [name, float(obs_b[0][4])-float(obs_v[0][4]), float(obs_v[0][4])-float(obs_i[0][4]) ]
	
	obs_colours.append(newline)


obs_umb = [line[1] for line in obs_colours]
obs_vmi = [line[2] for line in obs_colours]


umb = [line[1] for line in synth_colours]
vmi = [line[3] for line in synth_colours]

plt.figure()
plt.plot(vmi, umb, 'k--')
plt.plot(obs_vmi, obs_umb, 'o')
plt.xlabel('V-I')
plt.ylabel('U-B')
plt.show()
	























