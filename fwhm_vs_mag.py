#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Plot daophot instrumental mags as a function of airmass / fwhm"""

from matplotlib import pyplot as plt
import numpy as np

print 'here'

#read in table of photometry
#fpath = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/wrapper_standard_mags.tab'
#fpath = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/daophot_params_test/standard_mags.tab'
#fpath = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/aper_standard_mags.tab'
fpath = '/mirror2/scratch/hbarker/Orsola_2.3m_ANU/pyraf_standard_mags.tab'

tab = np.genfromtxt(fpath, dtype=None, comments='#', names=True)

for l in tab:
	print l

#choose a filter
tab = tab[ tab['filter']=='U']
#tab = tab[ tab['filter']=='B']
#tab = tab[ tab['filter']=='V']
#tab = tab[ tab['filter']=='I']





#tab = tab[ tab['star']=='MCT2019']
#tab = tab[ tab['star']=='111-1925']
#tab = tab[ tab['star']=='LSE44']
#tab = tab[ tab['star']=='G93-48']
#tab = tab[ tab['star']=='121-968']
tab = tab[ tab['star']=='wd1056']




mag = tab['mag']
airmass = tab['airmass']


plt.figure()
plt.plot(airmass, mag, 'o')
plt.xlabel('Airmass')
plt.ylabel('Mag')
plt.show() 	

