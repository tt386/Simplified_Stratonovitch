import numpy as np

from scipy import special
from scipy import integrate

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import os
import shutil

def PAP_Dist(x,t,L,Phi):
    Dist = np.zeros(len(x))

    Dist[x<0] = (1-Phi)*special.erf((-x[x<0])/np.sqrt(4*t))

    Dist[x>L] = (1-Phi)*special.erf((x[x>L]-L)/np.sqrt(4*t))

    return Dist


def Diffusion_Kernel(r,t):
    return (1./np.sqrt(4*np.pi*t)*np.exp(-r**2/(4*t)))


def Post_PAP_Dist(xlist,P,L,t,Phi,PAPDist):

    #np.pad(PAP_Dist(xlist,P,L,Phi),100,'constant',constant_values=(1-Phi))

    Diffusion = np.pad(Diffusion_Kernel(xlist,t),100,'constant',constant_values=0)

    Dist = np.convolve(PAPDist,Diffusion,mode='same')

    """
    for x in xlist:
        combinedfun = (lambda r:
            PAP_Dist(r,P,L) * Sense(x-r,t))
        #Split integral into upper and lower to eradiate error
        upper = integrate.quad(combinedfun,x,np.inf)[0]#integrate.quad(combinedfun,x,np.inf)[0]
        lower = integrate.quad(combinedfun,-np.inf,x)[0]#integrate.quad(combinedfun,-np.inf,x)[0]
        Dist.append(upper+lower) 
    """
    return np.asarray(Dist)




def AddRects(ax,ylow,xlow,L,t,P):
    #Left Rect
    # Create a Rectangle patch
    rect = patches.Rectangle((xlow, ylow), abs(xlow), ylow*1.1, linewidth=3, edgecolor='k', facecolor='white',zorder=1)
    # Add the patch to the Axes
    ax.add_patch(rect)

    #Right Rect
    #Left Rect
    # Create a Rectangle patch
    rect = patches.Rectangle((L, ylow), abs(xlow), ylow*1.1, linewidth=3, edgecolor='k', facecolor='white',zorder=1)
    # Add the patch to the Axes
    ax.add_patch(rect)


    Middlecol = '#808080'
    if t>P:
        Middlecol = 'white'

    #Middle Ret
    # Create a Rectangle patch
    rect = patches.Rectangle((0,ylow), L, ylow*1.1, linewidth=1, edgecolor='k', facecolor=Middlecol,zorder=1)
    # Add the patch to the Axes
    ax.add_patch(rect)




P = 0.2
L = 3
Phi = 1e-4

SaveDirName = "SavedFiles/L_%d_LogPhi_%d_P_%0.2f"%(L,np.log10(Phi),P)

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copy("Script.py",SaveDirName)

#########################################################################
xlist = np.linspace(-20,20+L/2,10000)
tlist = np.linspace(0,1,100)

SDist = np.ones(len(xlist)) * (1-Phi)
RDist = np.ones(len(xlist)) * Phi

PAPDist = np.pad(PAP_Dist(xlist,P,L,Phi),100,'constant',constant_values=(1-Phi))

#PAPDist = PAP_Dist(xlist,P,L,Phi)

PCount = 0
PPCount = 0 
for t in tlist:
    print(t)

    fig,ax = plt.subplots()
    AddRects(ax,1e-5,min(xlist),L,t,P)

    if t < P:
        PCount += 1    

        Dist = PAP_Dist(xlist,t,L,Phi)

        plt.semilogy(xlist,Dist,linewidth=5)
        plt.semilogy(xlist,RDist,linewidth=5)


    else:
        PPCount += 1
        Dist =  Post_PAP_Dist(xlist,P,L,t-P,Phi,PAPDist)

        plt.semilogy(xlist+1/2,Dist*(xlist[1]-xlist[0]),linewidth=5)
        plt.semilogy(xlist,RDist,linewidth=5)

    Dist[Dist<1e-5] = 0

    plt.ylim(1e-5,1.1)

    xwidth = 10
    plt.xlim(-xwidth,xwidth+L)

    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.tight_layout()

    plt.savefig(SaveDirName+"/t_%0.3f.png"%(t),bbox_inches='tight')

    if t<P:
        plt.savefig(SaveDirName+"/PAP_%0.3f.png"%(t),bbox_inches='tight')

    else:
        plt.savefig(SaveDirName+"/PostPAP_%0.3f.png"%(t),bbox_inches='tight')

    plt.close()

