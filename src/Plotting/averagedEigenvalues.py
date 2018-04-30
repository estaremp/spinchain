import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
from pylab import *
import numpy as np

#NOTE: ax1 and ax3 are commented, as well as
#the code drawing diagonals. Uncomment when
#adjusting energy levels graph with three subplots
#i.e. to higlight areas of the spectrum by rescaling it

#LOAD DATA FROM FILE
data=np.loadtxt('averagedEigenvalues.data')

#SET SIZE OF THE FIGURE
fig = plt.figure(figsize=(20,10))

#SET GRID OF SUBPLOTS
G = gridspec.GridSpec(3, 3)

#NUMBER OF STATES AS ARGUMENT
totalstates=int(sys.argv[1])

#ARRAY GOING FROM 1 TO NUMBER OF STATES
#EXCLUDING VACUUM STATE
x = iter(xrange(totalstates))

ax = plt.subplot(G[:,0:2])

ax.set_xlim([0,totalstates+1])
ax.set_ylim(-1*(np.amax(data)+0.1),(np.amax(data)+0.1))

ax.xaxis.set_major_locator(MaxNLocator(integer=True))

ax.set_ylabel('Averaged Energy',fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)

ax.plot(np.arange(1,totalstates+1),data, color='salmon', lw=2, marker='8')

#GENERATE THREE SUBPLOTS TO
#REPRESENT DIFFERENT REGIONS
#OF THE SPECTRUM
#ax1 = plt.subplot(G[0,2])
ax2 = plt.subplot(G[:,2])
#ax3 = plt.subplot(G[2,2])

#SET AXIS LABEL
ax2.set_ylabel('Averaged Energy',fontsize=25)

#INSIVISBLE THE X AXIS
#ax1.axes.get_xaxis().set_visible(False)
ax2.axes.get_xaxis().set_visible(False)
#ax3.axes.get_xaxis().set_visible(False)

#PLOT ALL LEVELS
for i in x:
    #ax1.hlines(data[i],1,2,lw=2)
    ax2.hlines(data[i],1,2,lw=2)
    #ax3.hlines(data[i],1,2,lw=2)

#**MODIFY DEPENDING ON SYSTEM*
#ADJUST THESE LIMTIS
#ax3.set_ylim(-1*(np.amax(data)+0.01),-0.05)ax2.set_ylim(-1*(np.amax(data)+0.01),(np.amax(data)+0.01))
#ax1.set_ylim(0.05,np.amax(data)+0.01)

#SET X AXIS LIMITS

#ax1.set_xlim(0.8,2.2)
ax2.set_xlim(0.8,2.2)
#ax3.set_xlim(0.8,2.2)

#INVISIBLE AXIS
#ax1.spines['bottom'].set_visible(False)
#ax1.spines['top'].set_visible(False)
#ax1.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['right'].set_visible(False)
#ax3.spines['top'].set_visible(False)
#ax3.spines['bottom'].set_visible(False)
#ax3.spines['right'].set_visible(False)

#INVISIBLE AXIS TICKS
#ax1.xaxis.tick_top()
ax2.xaxis.tick_top()
ax2.tick_params(labeltop='off')  # don't put tick labels at the top
#ax3.xaxis.tick_top()
#ax3.tick_params(labeltop='off')  # don't put tick labels at the top

#MOVE TICK POSITIONS TO LEFT
#ax1.yaxis.set_ticks_position('left')
ax2.yaxis.set_ticks_position('left')
#ax3.yaxis.set_ticks_position('left')

#SET MAXIMUM RANGE FOR AXIS TICKS
#ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
#ax3.yaxis.set_major_locator(MaxNLocator(integer=True))

#SET SIZE TICKS
#ax1.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
#ax3.tick_params(axis='y', labelsize=20)

#DRAW LITTLE DIAGONAL LINES SEPARATING SUBPLOTS
#d = .015  # how big to make the diagonal lines in axes coordinates
## arguments to pass plot, just so we don't keep repeating them
#kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
#ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
#
#
#kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
#ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
#
#kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
#ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
#
#kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
#ax3.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

plt.tight_layout()

#SAVE FIGURE AS A PNG FILE
plt.savefig('AveragedEigenvalues.png',transparent=False)

