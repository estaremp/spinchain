import numpy as np
from pylab import *

#LOAD DATA FILES
data = np.loadtxt('exmap.data',comments='#')

#ARGUMENTS PASSED
totaltime=float(sys.argv[1])

#SET SIZE FIGURE
plt.figure(figsize=(10,7))

ax = plt.subplot()

def Colormap(ax,lst):

    #STRUCTURE DATA FOR COLORMAP
    intensity = np.array(lst[:,1:])
    
    x, y = intensity.shape
    
    x1 = lst[:,0]
    y1 = range(0,y+1)
    x2,y2 = np.meshgrid(x1,y1)
    
    #COLORMAP
    mappable = plt.pcolormesh(x2,y2,np.swapaxes(intensity,0,1),cmap='rainbow',vmin=0,vmax=1)
    c = colorbar(mappable)
    c.ax.tick_params(labelsize=20)

    #SET X-AXIS TICKS AND LABEL
    ax.set_xlim([0,totaltime])
    ax.set_xlabel('$\mathrm{time \cdot J_{max}}$',fontsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.tick_params(axis='x', labelsize=20)
    
    #SET Y-AXIS TICKS AND LABEL
    ax.set_ylim([0,y])
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(0.5, end+0.5, 1))
    labels = [item.get_text() for item in ax.get_yticklabels()]
    i = 0
    while (i<(len(data[0,1:]))):
        labels[i]='spin'+str(i+1)
        i=i+1
    ax.set_yticklabels(labels,fontsize=20)

Colormap(ax,data)
#SAVE FIG
savefig('exmap.png',transparent=False)

