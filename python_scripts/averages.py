import sys
import numpy as np
from scipy import stats
import math

with open(sys.argv[1]) as e:
    raw = e.read()
    lines = map(lambda y: float(y.strip()),  raw.split('\n')[:-1])

results_per_test = int(sys.argv[2])-1
total_tests = len(lines)/results_per_test

result = [ 0 for y in range(results_per_test) ]

for j in range(total_tests):
    for m in range(results_per_test):
        result[m] += lines[j*results_per_test + m]

result = map(lambda y: y/total_tests, result)

file=open("averagedEigenvalues.data","w")

for line in result:
    file.write(str(line)+"\n")

file.close()
