import sys
import numpy as np
from scipy import stats
import math


###############################################
###### AVERAGE EIGENVALUES ####################
###############################################


#Open files and take lines not starting with '#'
with open(sys.argv[1]) as e:
    lines=[]
    for line in e:
        if not (line.startswith(" #")):
            lines.append(float(line.strip()))

#Number of lines per realisation
results_per_real = int(sys.argv[2])

#Number of realisations
total_real = int(sys.argv[3])
result = [ 0 for y in range(results_per_real) ]

for j in range(total_real):
    for m in range(results_per_real):
        result[m] += lines[j*results_per_real + m]

result = map(lambda y: y/total_real, result)

#Open file where to save result average
file=open("averagedEigenvalues.data","w")

file.write("#AVERAGED EIGENVALUES//STANDARD DEVIATION//ERROR OF THE MEAN\n")

for line in result:
    file.write(str(line)+"\n")

file.close()
