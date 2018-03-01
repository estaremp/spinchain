import matplotlib.pyplot as plt
import numpy as np
from pylab import *

#READ THE DATA FROM FILE
fidelity=np.loadtxt('dynamics.data')

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

#PLOT FIDELITY AGAINST INITIAL STATE
plt.plot(fidelity[:,0],fidelity[:,injection],color='gray',lw=2,label=r'$|\langle\Psi(t)\vert \psi_{o}\rangle|^2$')

#PLOT FIDELITY AGAINST MIRROR STATE
plt.plot(fidelity[:,0],fidelity[:,N-injection+3],color='black',ls=':',lw=2,label=r'$|\langle\Psi(t)\vert \psi_{M}\rangle|^2$')

#SET SIZE OF AXIS TICKS
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)

#SET LABEL NAMES
plt.xlabel('$\mathrm{time \cdot J}$',fontsize=25,color='black')
plt.ylabel('${\cal{F}}(t)$',fontsize=25,color='black')

#LEGEND
l=legend(loc=1,frameon=False,borderaxespad=0.,fontsize=20)

#SAVE AS PNG PICTURE
savefig('dynamics.png',transparent=False)
