import sys
import numpy as np
from scipy import stats
import math

###############################################
###### AVERAGE EOF ############################
###############################################

#Calculate averaged EOF value at t_A

#Open files and take lines not starting with '#'
with open(sys.argv[1]) as e:
    lines=[]
    for line in e:
        if not (line.startswith(" #")):
            lines.append(float(line.split()[1]))

tA=str(sys.argv[2])
offnoise=str(sys.argv[3])
diagnoise=str(sys.argv[4])

mean=np.mean(lines)
stdev=np.std(lines)
sterror=stats.sem(lines,ddof=0)

#Open file where to save result average
file=open("averagedEOFtA.data","a")

file.write("#TIME//OFF-DIAG NOISE//DIAG NOISE//AVERAGED EOF//STANDARD DEVIATION//STANDARD ERROR \n")

#and save
file.write(str(tA) + '  ' + str(offnoise) + '  ' + str(diagnoise) + '  ' + str(mean) + '  ' +  str(stdev) + '  ' +  str(sterror) + '\n')

file.close()





