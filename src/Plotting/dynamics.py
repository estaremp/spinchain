import matplotlib.pyplot as plt
import numpy as np
import re
import sys
from pylab import *

#CHANGE FORMAT FOR THE COMPLEX NUMBERS FORTRAN->PYTHON
input=open('dynamics.data','r')
output=open('dynamics_formatted.data','a')
for line in input:
    line=re.sub(r'\(([^,\)]+),([^,\)]+)\)', r'\1+\2j', line.rstrip())
    line=re.sub(r'\+-',r'-',line)
    output.write(line+'\n')
output.close()

#READ THE DATA FROM FILE
fidelity = np.loadtxt('dynamics_formatted.data',dtype=complex,comments='#')

#SET SIZE FIGURE
plt.figure(figsize=(10,7))

#ARGUMENTS PASSED
totaltime=float(sys.argv[1])
injection=int(sys.argv[2])
N=int(sys.argv[3])

#CREATE PLOT
ax1=plt.subplot()

#DEFINE RANGES
ax1.set_xlim([0,totaltime])
ax1.set_ylim([0,1])

#PLOT FIDELITY AGAINST INITIAL STATE (CHANGE WHENEVER)
f=(np.absolute((fidelity[:,2])))**2
plt.plot(np.absolute(fidelity[:,0]),f,color='gray',lw=2,label=r'$|\langle\Psi(t)\vert \psi_{o}\rangle|^2$')

#PLOT FIDELITY AGAINST DIFFERENT STATE
fa=(np.absolute((fidelity[:,4])))**2
plt.plot(np.absolute(fidelity[:,0]),fa,color='black',ls=':',lw=2,label=r'$|\langle\Psi(t)\vert \psi_{A}\rangle|^2$')

#SET SIZE OF AXIS TICKS
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)

#SET LABEL NAMES
plt.xlabel('$\mathrm{time \cdot J_{max}}$',fontsize=25,color='black')
plt.ylabel('${\cal{F}}(t)$',fontsize=25,color='black')

#LEGEND
l=legend(loc=1,frameon=False,borderaxespad=0.,fontsize=20)

#SAVE AS PNG PICTURE
savefig('dynamics.png',transparent=False)
