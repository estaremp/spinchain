import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pylab import *

en1=np.loadtxt('eigenvalues.data')

fig = plt.figure(figsize=(10,7))

ax1=plt.subplot()

totalstates=int(sys.argv[1])

ax1.set_xlim([-2,totalstates])
ax1.set_ylim([-10,10])
ax1.set_ylabel('Energy',fontsize=25)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)

ax1.plot(en1, color='red', lw=3, marker='8')

plt.savefig('eigenvalues.png',transparent=True)

