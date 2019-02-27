import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import re

#CHANGE FORMAT FOR THE COMPLEX NUMBERS FORTRAN->PYTHON
input=open('dynamics.data','r')
output=open('dynamicse_formatted.data','a')
for line in input:
    line=re.sub(r'\(([^,\)]+),([^,\)]+)\)', r' \1+\2j ', line.replace(' ',''))
    line=re.sub(r'\+-',r'-',line)
    output.write(line+'\n')
output.close()

#READ THE DATA FROM FILE
fidelity = np.loadtxt('dynamicse_formatted.data',dtype=complex,comments='#')
entanglement=np.loadtxt('eof.data',comments='#')

#SET SIZE FIGURE
plt.figure(figsize=(10,7))

#ARGUMENTS PASSED
totaltime=float(sys.argv[1])
injection=int(sys.argv[2])

#CREATE PLOT
ax=plt.subplot()

#DEFINE RANGES
ax.set_xlim([0,totaltime])
ax.set_ylim([0,1])

#DEFINE STATE TO MEASURE FIDELITY AGAINST TO
f=(np.absolute((fidelity[:,9])+(fidelity[:,10])+(fidelity[:,14])+(fidelity[:,15]))/2.)**2 

#PLOT FIDELITY AGAINST ANY STATE F
plt.plot(np.absolute(fidelity[:,0]),f,color='red',lw=2,ls='-',markevery=500,label=r'$|\langle\Psi(0)|\Psi(t)\rangle|^2$')

#PLOT EOF
plt.plot(entanglement[:,0],entanglement[:,1],color='limegreen',lw=2,label='$EOF$')

ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)

#SET LABEL NAMES
plt.xlabel('$\mathrm{time \cdot J_{max}}$',fontsize=25,color='black')

#LEGEND
l=legend(loc=1,frameon=False,borderaxespad=0.,fontsize=20)

#SAVE AS PNG PICTURE
savefig('eof.png',transparent=False)
