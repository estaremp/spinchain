import sys
import numpy as np
from scipy import stats
import math

#Look for max value for each realisation (over 100) for each disorder and average them
lines=[]

#open file list of realisations for a noise value
with open(sys.argv[1]) as f:
    for line in f:
        if not (line.startswith(" #")):
            cols = line.split()
            lines.append(float(cols[1])) #append second column (first is time) as list

total_real = int(sys.argv[2])
lines_per_test = len(lines)/total_real #number of points the dynamics

maxval=[]

time_win=str(sys.argv[3])
offnoise=str(sys.argv[4])
diagnoise=str(sys.argv[5])

#look for all the max values at each realisation
for j in range(total_real):
    new_lines = lines[(j*lines_per_test):((j+1)*lines_per_test)]
    maxval.append(max(new_lines))
#print maxval

mean=np.mean(maxval)
stdev=np.std(maxval)
sterror=stats.sem(maxval,ddof=0)

#and save
file = open("averagedMaxEOF.data","a")

#print headings
file.write("#TIME WINDOW//OFF-DIAG NOISE//DIAG NOISE//AVERAGED MAX EOF//STANDARD DEVIATION//STANDARD ERROR \n")

#and save
file.write(str(time_win) + '  ' + str(offnoise) + '  ' + str(diagnoise) + '  ' + str(mean) + '  ' +  str(stdev) + '  ' +  str(sterror) + '\n')
file.close()
