#!/bin/bash
#Run daophot photometry
#needs to be run from the directory
#when called by daophot wrapper:  $1 = frame number, $2=fwhm, $3=filepath



# Get rid of any junk files which DAOPHOT has left lying around.
rm -f *jnk.sdf	

in_num=$1   #frame number
fwhm=$2     #measured fwhm
f=$3        #filepath



#convert each fits file to a ndf file
#for f in $fitsfiles
#do

	#full filenames are too long and need shortening
        #eg. T2M3Im-20120515-174021-0151.fits -> 0151.fits
	filename=$(basename "$f")
	fileparts="${filename#*.}"
	framenum="${fileparts#*-}"
	framenum="${framenum%%_*}"
	
	
	if [ "$in_num" != "$framenum" ]; then
	continue
	fi
	
	#copy the file to a new shorter name
	newname=$framenum".fits"
	cp $f $newname

	#convert the fits file to a ndf
	echo $f
	"/star/bin/convert/fits2ndf" "$newname" "${newname/.fits/.sdf}"
	echo

	#delete the copied fits file (with the shorter name)
	rm $newname

	#print the new file name
	filename="${newname/.fits/.sdf}"
	echo "Created " $filename


	#Get filename without the .sdf extension.
	extension="${filename##*.}"
	fileroot="${filename%.*}"


	# Get rid of any junk files which DAOPHOT has left lying around.
	rm -f *jnk.sdf	 
	rm -f "$fileroot.ap"
	rm -f "$fileroot.coo"
	rm -f "$fileroot.lst"
	rm -f "$fileroot.psf"
	    

	
	#bash can't handle non-integer math, need to use bc
	fwhm0_1=$(echo $fwhm*0.1 | bc)
	fwhm0_3=$(echo $fwhm*0.3 | bc)
	fwhm0_5=$(echo $fwhm*0.5 | bc)
	fwhm0_7=$(echo $fwhm*0.7 | bc)
	fwhm2=$(echo $fwhm*2.0 | bc)
	fwhm2_5=$(echo $fwhm*2.5 | bc)
	fwhm3=$(echo $fwhm*3.0 | bc)
	fwhm4=$(echo $fwhm*4.0 | bc)
	fwhm5=$(echo $fwhm*5.0 | bc)
	fwhm6=$(echo $fwhm*6.0 | bc)
	fwhm7=$(echo $fwhm*7.0 | bc)
	fwhm8=$(echo $fwhm*8.0 | bc)
	fwhm9=$(echo $fwhm*9.0 | bc)
	fwhm10=$(echo $fwhm*10.0 | bc)
	

	
	#create fwhm commands to feed into daophot
	#maximum allowed values:
	#fw < 15
	#fi <10
	#ps<35
	#is < OS
	#os < n/a	
	

	fwhm_cmd='FWHM='$fwhm         #FW (pixels) - for which FIND is optimised
	inner_rad='IS='$fwhm4         #IS - inner sky radius
	outer_rad='OS='$fwhm5         #OS - outer sky radius


	if [ $(echo "$fwhm2>10.0"|bc) -eq 1 ]; then
		fit_cmd='FIT=10.0'
	else
		fit_cmd='FIT='$fwhm2
	fi
	
	if [ $(echo "$fwhm6>35.0"|bc) -eq 1 ]; then 
		psf_cmd='PSF=35.0'
	else
		psf_cmd='PSF='$fwhm4
	fi
	


#need to write these to the files
	daophot_opt=$PWD/daophot.opt
	rm -f $daophot_opt
	echo $fwhm_cmd >> $daophot_opt
	echo $fit_cmd >> $daophot_opt
	echo $psf_cmd >> $daophot_opt
	echo "READ=3.278" >> $daophot_opt
	echo "GAIN=1.1" >> $daophot_opt
	echo "TH=5.0" >> $daophot_opt
	echo "AN=1" >> $daophot_opt          #1=Gaussian, 2=Moffat
	echo "LOWBAD=5" >> $daophot_opt
	echo "HIBAD=50000" >> $daophot_opt
	echo "WATCH=-1" >> $daophot_opt      #-1 non-interactive and the star is used anyway, but hopefully is removed by bad pixel rejection later on, 0 asks about stars with bad pixels, 1 needs manual input
	echo "VAR=0" >> $daophot_opt        #0=fixed psf, 1=varies linearly with the position, 2=varies in degrees with position
	echo "EX=5">> $daophot_opt #extra cleaning passes to remove bad pixels from psf tables
	
	echo "Written " $daophot_opt
	echo
	
	photo_opt=$PWD/photo.opt
	rm -f $photo_opt
	echo "A1="$fwhm0_2 >> $photo_opt
	echo "A2="$fwhm0_5 >> $photo_opt
	echo "A3="$fwhm >> $photo_opt
	echo "A4="$fwhm2 >> $photo_opt
	echo "A5="$fwhm2_5 >> $photo_opt
	echo "A6="$fwhm3 >> $photo_opt
	echo "A7"$fwhm4 >> $photo_opt
	echo "A8="$fwhm5 >> $photo_opt
	echo "A9=" >> $photo_opt
	echo "AA=" >> $photo_opt
	echo "AB=" >> $photo_opt
	echo "AC=" >> $photo_opt
	echo $inner_rad >> $photo_opt
	echo $outer_rad >> $photo_opt
	
	echo "Written " $photo_opt
	echo
	
	#FOR ALLSTAR, the manual says OS and IS need not, and probably shouldn't be the same as those used in daophot. Here, OS should be large compared to the fitting radius, FI (probably close to PS). IS can be small (~2 pixels)
	inner_rad_allstar='IS='$fwhm
	outer_rad_allstar='OS='$fwhm6
	
	allstar_opt=$PWD/allstar.opt
	rm -f $allstar_opt
	echo $fit_cmd >> $allstar_opt
	echo $inner_rad_allstar >> $allstar_opt
	echo $outer_rad_allstar >> $allstar_opt
	echo "WATCH=1.0" >> $allstar_opt
	echo "redet=1" >> $allstar_opt
	
	echo "Written " $allstar_opt
	echo
     

#  Write filename to screen.
         echo "***"
         echo "*** Processing file $filename ***"
         echo "***"
         

         
#	start daophot
	source /star/bin/daophot/daophot.csh


#  Run DAOPHOT with the appropriate commands. LEAVE UNINDENTED
# For pick, I want 50 objects no fainter than 25 mags
	/star/bin/daophot/daophot <<__DAOPHOT-END__
ATTACH $fileroot
FIND
1,1
$fileroot
Y
PHOT
YES
A1=3.0
A2=4.0
A3=5.0
A4=6.0
A5=7.0
A6=8.0
A7=9.0
A8=10.0
$inner_rad
$outer_rad
YES
$fileroot.coo
$fileroot.ap
PICK
$fileroot.ap
50,25
$fileroot.lst
__DAOPHOT-END__


#filter the lst file using a python script
	python /mirror2/scratch/hbarker/Orsola_2.3m_ANU/scripts/filter_pick.py -f $fileroot -p $f


echo 'Calculating PSF'
#resume daophot with PSF
	/star/bin/daophot/daophot <<__DAOPHOT-END__
ATTACH $fileroot
PSF
$fileroot.ap
$fileroot.lst
$fileroot.psf
EXIT
__DAOPHOT-END__




























