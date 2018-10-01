#!/bin/bash
#Convert star subtracted sdf files to fits


#launch convert
convert

#list all *s.sdf files
files=/mirror2/scratch/hbarker/Orsola_2.3m_ANU/daophot_params_test/working_dir/*s.sdf

for f in $files
do
	echo $f
	
	filename=$(basename "$f")
	new_filename=${filename/.sdf/.fits}
	new_path=/mirror2/scratch/hbarker/Orsola_2.3m_ANU/daophot_params_test/subtracted_imgs/$new_filename
	/star/bin/convert/ndf2fits $f $new_path

done

















