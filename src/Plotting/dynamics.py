import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pylab import *

fidelity2=np.loadtxt('dynamics.data')
fidelity3=np.loadtxt('dynamics.data')

#fidelity4=np.loadtxt('dynamics4.data')
#fidelity5=np.loadtxt('dynamics5.data')
#fidelity6=np.loadtxt('dynamics6.data')
#entanglement=np.loadtxt('dynamics.data')
#analytical=np.loadtxt('EoFanaly.data')

plt.figure(figsize=(20,7))

totaltime=float(sys.argv[1])
injection=int(sys.argv[2])
ax1=plt.subplot()

ax1.set_xlim([0,totaltime])
ax1.set_ylim([0,1])

ax1.spines['bottom'].set_color('black')
ax1.spines['top'].set_color('black')
ax1.spines['left'].set_color('black')
ax1.spines['right'].set_color('black')
ax1.tick_params(axis='x', colors='black')
ax1.tick_params(axis='y', colors='black')

plt.plot(fidelity2[:,0],fidelity2[:,injection],color='red',lw=2,label=r'$|\langle10100|\Psi(t)\rangle|^2$')
#plt.plot(fidelity3[:,0],fidelity3[:,16],color='limegreen',lw=2,label=r'$|\langle00011|\Psi(t)\rangle|^2$')
#plt.plot(fidelity4[:,0],fidelity4[:,1],color='blue',ls='--',lw=2,label=r'$|\langle\Psi(4)|\Psi(t)\rangle|^2$')
#plt.plot(fidelity5[:,0],fidelity5[:,1],color='purple',ls=':',lw=2,label=r'$|\langle\Psi(5)|\Psi(t)\rangle|^2$')
#plt.plot(fidelity6[:,0],fidelity6[:,1],color='orange',ls='.',lw=2,label=r'$|\langle\Psi(6)|\Psi(t)\rangle|^2$')
#plt.plot(entanglement[:,0],entanglement[:,1],color='limegreen',lw=2,label='$EoF_{N}$')
#plt.plot(analytical[:,0],analytical[:,1],color='black',lw=3,ls=':',label='$EoF_{A}$')

ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
plt.xlabel('$\mathrm{time \cdot \Delta}$',fontsize=25,color='black')
plt.ylabel('$\mathrm{Fidelity}$',fontsize=25,color='black')

l=legend(bbox_to_anchor=(0., 1., .85, 0.),ncol=1,frameon=False,borderaxespad=0.,fontsize=20)

savefig('dynamics.png',transparent=True)
