#!/bin/bash
#Run daophot photometry: the same as run_daophot except with no FIND command
#needs to be run from the directory
#when called by daophot wrapper:  $1 = frame number, $2=fwhm, $3=filepath



# Get rid of any junk files which DAOPHOT has left lying around.
rm -f *jnk.sdf	

in_num=$1
fwhm=$2
f=$3

#manual frame number input
#echo "Enter frame number: "
#read in_num
#echo

#convert each fits file to a ndf file
#for f in $fitsfiles
#do

	#filenames are often too long and need shortening

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

	echo $f
	"/star/bin/convert/fits2ndf" "$newname" "${newname/.fits/.sdf}"
	echo

	#delete the copied original file
	rm $newname


	filename="${newname/.fits/.sdf}"
	echo "Created " $filename



#  Get filename without the .sdf extension.
	extension="${filename##*.}"
	fileroot="${filename%.*}"
	    
         
	#Prompt user for FWHM in manual input
	#echo -n "Enter the FWHM in pixels and press [ENTER]:"
	#read fwhm
	
	
	#bash can't handle non-integer math, need to use bc
	fwhm2=$(echo $fwhm*2.0 | bc)
	fwhm3=$(echo $fwhm*3.0 | bc)
	fwhm4=$(echo $fwhm*4.0 | bc)
	fwhm5=$(echo $fwhm*5.0 | bc)
	fwhm6=$(echo $fwhm*6.0 | bc)
	fwhm7=$(echo $fwhm*7.0 | bc)
	fwhm8=$(echo $fwhm*8.0 | bc)
	fwhm9=$(echo $fwhm*9.0 | bc)
	fwhm10=$(echo $fwhm*10.0 | bc)
	
	
#create fwhm commands to feed into daophot
	fwhm_cmd='FWHM='$fwhm
	fit_cmd='FIT='$fwhm2
	psf_cmd='PSF='$fwhm4
	
	inner_rad='IS='$fwhm4
	outer_rad='OS='$fwhm5




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
	echo "A1=3.0" >> $photo_opt
	echo "A2=4.0" >> $photo_opt
	echo "A3=5.0" >> $photo_opt
	echo "A4=6.0" >> $photo_opt
	echo "A5=7" >> $photo_opt
	echo "A6=8" >> $photo_opt
	echo "A7=9" >> $photo_opt
	echo "A8=10" >> $photo_opt
	echo "A9=" >> $photo_opt
	echo "AA=" >> $photo_opt
	echo "AB=" >> $photo_opt
	echo "AC=" >> $photo_opt
	echo $inner_rad >> $photo_opt
	echo $outer_rad >> $photo_opt
	
	echo "Written " $photo_opt
	echo
	
	
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
#	source /star/bin/daophot/daophot.csh


#  Run DAOPHOT with the appropriate commands. LEAVE UNINDENTED
# For pick, I want 50 objects no fainter than 25 mags
	/star/bin/daophot/daophot <<__DAOPHOT-END__
ATTACH $fileroot
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
2000
2000
3.278
1.1
5
50000
$fileroot.ap
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




























