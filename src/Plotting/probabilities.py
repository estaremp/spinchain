import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter
import numpy as np
from pylab import *

#eigenstates probabilities
eigenstates=np.loadtxt('probabilities.data',skiprows=1)

fig = plt.figure(figsize=(12,7))

ax=plt.subplot()

nvecs=int(sys.argv[1])

i=list(range(1,nvecs))

ax.set_xlim([1-(0.1*nvecs),nvecs+(0.1*nvecs)])
ax.set_ylim([0,1])
ax.set_ylabel(r'$|c_{i,m}|^2$',fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_xlabel('|i>',fontsize=25)
ax.plot(i,eigenstates)


plt.title('Eigenstate Site-Probabilities')

fig.set_tight_layout(True)
savefig('probabilities.png',transparent=True)
