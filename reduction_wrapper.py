#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""wrapper to call the scripts needed to reduce Orsola's 2.3m data"""

import subprocess as sp
import os
import sys
from datetime import datetime


startTime = datetime.now()

print "Overscan correcting..."
os.system("overscan_correct.py")
print
print "Trimming overscan region."
os.system("trim.py")
print
print "Creating master biases"
os.system("master_bias.py")
print
print "De-biasing frames"
os.system("debias.py")
print
print "Creating master flat"
os.system("master_flat.py")
print
print "Flat correcting..."
os.system("deflat.py")
print
print "Trimming images"
os.system("circle_trim.py")
print



finTime = datetime.now() - startTime
print "Elapsed time: %s" %finTime














