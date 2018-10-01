#!/bin/bash
#Run daophot photometry
#needs to be run from the directory
#when called by daophot wrapper:  $1 = frame number, $2=fwhm


# Get rid of any junk files which DAOPHOT has left lying around.
rm -f *jnk.sdf	

#read in files
fitsfiles=$PWD/*circletrim.fits


#manual frame number input
echo "Enter frame number: "
read in_num
echo

#convert each fits file to a ndf file
for f in $fitsfiles
do

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


	# Get rid of any junk files which DAOPHOT has left lying around.
	rm -f *jnk.sdf	 
	rm -f "$fileroot.ap"
	rm -f "$fileroot.coo"
	rm -f "$fileroot.lst"
	rm -f "$fileroot.psf"
	    
         
	#Prompt user for FWHM in manual input
	echo -n "Enter the FWHM in pixels and press [ENTER]:"
	read fwhm

	
	#bash can't handle non-integer math, need to use bc
	fwhm2=$(echo $fwhm*2.0 | bc)
	fwhm4=$(echo $fwhm*4.0 | bc)
	fwhm5=$(echo $fwhm*5.0 | bc)
	
	
#create fwhm commands to feed into daophot
#create fwhm commands to feed into daophot
#maximum allowed values:
#fw < 15
#fi <10
#ps<35
#is < OS
#os < n/a	
	
	
#create fwhm commands to feed into daophot
	fwhm_cmd='FWHM='$fwhm
	inner_rad='IS='$fwhm4
	outer_rad='OS='$fwhm5


	if [ $(echo "$fwhm2>10.0"|bc) -eq 1 ]; then
		fit_cmd='FIT=10.0'
	else
		fit_cmd='FIT='$fwhm2
	fi
	
	if [ $(echo "$fwhm4>35.0"|bc) -eq 1 ]; then 
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
	echo "AN=1" >> $daophot_opt          #Gaussian
	echo "LOWBAD=5" >> $daophot_opt
	echo "HIBAD=50000" >> $daophot_opt
	echo "WATCH=1" >> $daophot_opt      #0 uses all stars with no bad pixels, 1 needs manual input
	echo "VAR=1" >> $daophot_opt        #0=fixed psf, 1=varies linearly with the position, 2=varies in degrees with position
	
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
	
	
	allstar_opt=$PWD/allstar.opt
	rm -f $allstar_opt
	echo $fit_cmd >> $allstar_opt
	echo $inner_rad >> $allstar_opt
	echo $outer_rad >> $allstar_opt
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
EXIT
__DAOPHOT-END__


#PSF
#$fileroot.ap
#$fileroot.lst
#$fileroot.psf
#EXIT
#__DAOPHOT-END__

#	Run allstar
#	/star/bin/daophot/allstar <<__ALLSTAR-END__
#y
#$fileroot
#$fileroot.psf
#$fileroot.ap
#$fileroot.als
#y
#__ALLSTAR-END__
	

#	Sort als file by increasing x coordinate, not renumbering stars
#/star/bin/daophot/daophot <<__SORT-END__
#SORT
#2
#$fileroot.als
#$fileroot.srt
#$fileroot.srt
#n
#EXIT
#__SORT-END__
	

#open the sorted file in gedit
#/usr/bin/gedit $fileroot.srt
	

break

done
































