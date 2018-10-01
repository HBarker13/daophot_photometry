#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""script called by run_daophot.sh (which is called by daophot_wrapper.py). Takes the .lst file created by daophot > PICK
and removes all stars except those in pick_stars.tab (user defined stars to use as references for doing photometry"""



import argparse
import os
import sys
import time


#read the frame being processed
parser = argparse.ArgumentParser(description="Collect the frame number")
parser.add_argument('-f','--frame_num', help="frame number", required=True)
parser.add_argument('-p','--filepath', help="filepath", required=True)
args= parser.parse_args()

frame_num = args.frame_num
fpath = args.filepath


#get the paths to the lst file and coordinate list
dirpath, _ = fpath.split('/T2', 1)
lst_path = os.getcwd() + '/' + frame_num + '.lst'
if not os.path.exists(lst_path):
	print "lst path doesn't exist"
	print lst_path
	sys.exit()
	
coords_path = dirpath + '/pick_stars.tab'
if not os.path.exists(coords_path):
	print "Coordinates path doesn't exist"
	print coords_path
	sys.exit() 	




#read in the files
with open(lst_path, 'r') as f:
	#number, x, y, ...
	lst = [line.strip().split() for line in f]
	colnames = lst[0]
	info_line = lst[1]
	#skip empty line
	lst = lst[3:] #skip the first lines with the column names, empty lines...


with open(coords_path, 'r') as f:
	#x, y
	coords = [line.strip().split() for line in f]
	coords = [line for line in coords if len(line)==2]

	
#remove any lines from the lst file that aren't what we want
completed = []
new_tab = []
for line in lst:
	x = float( line[1] )
	y = float( line[2] )
	
	for ref in coords:	
		
		if len(ref)!=2:
			print 'Error with pick_stars file'
			print coords_path
			for l in coords:
				print l
			raw_input('EXIT')

		
		if abs( x - float(ref[0]) )<10.0 and abs( y -float(ref[1]))<10.0:
			new_tab.append(line)
			completed.append(ref)
			


if len(new_tab)!=len(coords):
	print "Not all reference stars were found"
	print fpath
	print frame_num
	for c in coords:
		if c not in completed:
			print c
	
	
	

#save the new file: daophot is very picky about the format	
with open(lst_path, 'w+') as f:
	f.write(' ')
	for name in colnames:
		f.write(name)
		if name!=colnames[-1]:
			f.write('   ')
	f.write('\n')
		
	f.write('  '+info_line[0]+' '+info_line[1]+' '+info_line[2]+'    '+info_line[3]+' '+info_line[4]+'   '+info_line[5]+'    '+info_line[6]+'    '+info_line[7]+'    '+info_line[8]+'    '+info_line[9])
	f.write('\n')	
	f.write('\n')

	for line in new_tab:
		f.write('  '+line[0]+' '+line[1]+' '+line[2]+'   '+line[3]+'   '+line[4])
		f.write('\n')	
print 'Created', lst_path


	
	
	
	
	
	

	 

