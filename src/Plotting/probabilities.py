import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter
import numpy as np
from pylab import *

#READ DATA FROM FILE
eigen=np.loadtxt('probabilities.data',skiprows=2,comments='#')

#SET SIZE FIGURE
fig = plt.figure(figsize=(10,7))

ax=plt.subplot()

#SET LIMITS
ax.set_xlim([0,len(eigen)+1])
ax.set_ylim([0,1])

#SET LABEL NAMES
ax.set_ylabel(r'${\cal{P}}_{i,k}$',fontsize=25)
ax.set_xlabel('site number, $i$',fontsize=25)

#SET SIZE OF THE AXIS TICKS
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)

#SET X TICKS TO BE INTEGERS
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

#PLOT
ax.plot(np.arange(1,len(eigen)+1),eigen)

fig.tight_layout()

#SAVE FIGURE IN A PNG FILE
savefig('eigenstates.png',transparent=False)
