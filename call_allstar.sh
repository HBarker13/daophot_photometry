#!/bin/bash
#Run daophot photometry
#needs to be run from the directory
#when called by daophot wrapper:  $1 = frame number, $2=fwhm, $3=filepath


in_num=$1
f=$2


#filenames are often too long and need shortening
filename=$(basename "$f")
fileparts="${filename#*.}"
framenum="${fileparts#*-}"
framenum="${framenum%%_*}"
	
fileroot=$framenum	


#redo pick to get the original lst file
#NOT SURE I NEED TO DO THIS
#	/star/bin/daophot/daophot <<__DAOPHOT-END__
#ATTACH $fileroot
#PICK
#$fileroot.ap
#50,25
#$fileroot.lst
#$fileroot.lst
#__DAOPHOT-END__


echo 'Running allstar'
echo $fileroot
	/star/bin/daophot/allstar <<__ALLSTAR-END__
y
$fileroot
$fileroot.psf
$fileroot.ap
$fileroot.als
$fileroot_sub
__ALLSTAR-END__
	


