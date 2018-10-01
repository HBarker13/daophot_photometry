# process_telescope


All these scripts can be called at once using reduction_wrapper.py

1. overscan_correct.py.
Reads in files with format T2M3Im-yymmdd.time-00##.fits (bias, flat and science frames).
Looks at the overscan region (the leftmost two columns), calculates the mean pixel value, and subtracts this from all pixels. 
Makes newfiles with the format T2M3Im-yymmdd.time-00##_ovscorr.fits


2. trim.py
Trims the overscan region


3. master_bias.py
Reads in overscan corrected bais frames from each night and calculates the median value for each pixel, and saves this as a master bias. 


4. debias.py
Reads in the master bias for each night, and debiases the overscan corrected science and flat files.
Creates files with the format T2M3Im-yymmdd.time-00##_ovscorr_debiased.fits

NB. debias.py removes the T2M3Im-yymmdd.time-00##_ovscorr.fits file. This can be commented out


5. master_flat.py
Creates a master flat in each filter (per night) from the median of the overscan corrected, debiased flat files. NOTE: flat frames in each filter were not taken every night.

Night1 flats: I, R and V band
Night2 flats: U and B band
Night 3: B, bessel z and I band

The two I band master flats (nights 1 and 3) are combined into a single flat.


6. deflat.py
Flat-field corrects frames of with the matching filter (from any night). 


7. circle_trim.py
Takes the fully processed image and sets any pixels outside the circular aperture to NaN. This can help stop pick/psf in daophot mistaking the noise as stars. The circle is not a perfect fit, but works well enough.



Misc: clear_all.py
Deletes all processed files (ie. overscan corected, debiased, flat corrected and master flat and bias frames) so the Sorted directory contains only the raw data.
