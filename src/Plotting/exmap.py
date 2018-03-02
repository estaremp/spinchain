import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter



data = np.loadtxt('exmap.data')


def Colormap(lst):
    
    fig, ax = plt.subplots(figsize=(20,7))
    
    intensity = np.array(lst)
    
    x, y = intensity.shape
    
    x1 = range(x) # changed this also
    y1 = range(y+1) # changed this also
    
    x2,y2 = np.meshgrid(x1,y1)
    

    plt.pcolormesh(x2,y2,np.swapaxes(intensity,0,1),cmap='rainbow',vmin=0,vmax=1) # Transpose of intensity
    #plt.imshow(lst[x2,y2], cmap='RdBu', vmin=0, vmax=1,
    #           extent=[x2.min(), x2.max(), y2.min(), y2.max()],
    #           interpolation='nearest', origin='lower')
    plt.colorbar()
    labels_x = []
    
    for i in range(0,len(data[0,1:])+1,1):
        labels_x.append('%.2E'%(i*(data[-1,0])/len(data[0,1:])))
    ax.set_xticklabels(labels_x)
    ax.set_xlabel('time $\Delta$',fontsize=20)
    ax.set_ylabel('spin',fontsize=20)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)
    plt.savefig('exmap.png',transparent=False)


Colormap(data[:,1:])
