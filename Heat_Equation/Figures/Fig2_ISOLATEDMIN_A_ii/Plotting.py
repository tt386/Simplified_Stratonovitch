import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as patches
#plt.rcParams['text.usetex'] = True

import numpy as np
import time

from scipy.stats import linregress

from scipy.optimize import curve_fit

import sys

import scipy.integrate as integrate
import scipy.special as special

from scipy.signal import argrelextrema

from scipy.signal import savgol_filter

starttime = time.time()
################################
##ArgParse######################
################################
import os.path

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The Directory %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

from argparse import ArgumentParser

parser = ArgumentParser(description='Plotting')
parser.add_argument('-d','--directory',help='The directory of the data')
args = parser.parse_args()

###############################
##Extract Data#################
###############################
filename = 'datafile.npz'
SlopeMatrix = []
OSRRatioList = []
minmatrix = []

#Find list of all the datafiles
tempdirlist = os.listdir(args.directory)
dirlist = []
for i in tempdirlist:
    if os.path.isdir(os.path.join(args.directory,i)):
        dirlist.append(os.path.join(args.directory,i))



etaList = 0
timetaken = 0
PAP = 0
OSRNumList = []

Phi = 0
dx=0

OSRNum = 0

LeftMostMinima = []
RightMostMinima = []

######################
###Downloading Data###
######################
with np.load(os.path.join(str(args.directory), filename)) as data:
    SlopeMatrix.append(data['RDifferenceList'])
    OSRNumList.append(data['OSRNum'])
    etalist = data["etalist"]
    timetaken = data["timetaken"]
    PAP = data["PAP"]
    Phi = data["Phi"]
    OSRSeparation = data["OSRSeparation"]
    
    xlist = data['xlist']
    PAPMatrix = data["PAPMatrix"]
    YearMatrix = data["YearMatrix"]
    
    timetaken = data['timetaken']
    print("TIME TAKEN:",timetaken)

################
###Find Minima##
################

for p in range(len(etalist)):
    eta = etalist[p]
    PAPylist = PAPMatrix[p]
    Yearylist = YearMatrix[p]

    rlist = Phi * np.ones(len(xlist))

    Initlist = 1-rlist

    plt.figure()
    plt.semilogy(xlist,Initlist,label='Initial S')
    plt.semilogy(xlist,rlist,label='Initial R')
    plt.semilogy(xlist,PAPylist,label='PAP')
    plt.semilogy(xlist,Yearylist,label='Year')
    plt.semilogy(xlist,rlist/(rlist+Yearylist)-rlist,label='Integrand')
    plt.legend(loc='upper left')
    plt.title("Diffusion Coeff: "+ '{0:.5f}'.format(eta))
    
    plt.xlim(-1.5,2.5)

    plt.xlabel("Space")
    plt.grid()
    plt.savefig(str(args.directory) + "/Eta_" + '{0:.5f}'.format(eta) + ".png")
    plt.close()    
    

    L = 1/np.sqrt(etalist[0])

    xlist *= L

    xlow = -1.5 * L
    xupp = 2.5 *L

    yupp = 10**0.5
    ylow = Phi/10#**2

    xticks = []
    for x in np.arange(int(xlow),int(xupp)+1):
        if x%10 == 0:
            xticks.append(x)


    def AddRects(ax):
        #Left Rect
        # Create a Rectangle patch
        rect = patches.Rectangle((xlow, np.log10(ylow)), abs(xlow), np.log10(10), linewidth=3, edgecolor='k', facecolor='white',zorder=1)
        # Add the patch to the Axes
        ax.add_patch(rect)

        #Right Rect
        #Left Rect
        # Create a Rectangle patch
        rect = patches.Rectangle((L, np.log10(ylow)), abs(xlow), np.log10(10), linewidth=3, edgecolor='k', facecolor='white',zorder=1)
        # Add the patch to the Axes
        ax.add_patch(rect)

        #Middle Ret
        # Create a Rectangle patch
        rect = patches.Rectangle((0, np.log10(ylow)), L, np.log10(10), linewidth=1, edgecolor='k', facecolor='#808080',zorder=1)
        # Add the patch to the Axes
        ax.add_patch(rect)

    def Setup(plt,ax,xlist,S,R):
        #plt.plot(xlist,np.log10(S),linewidth=10)
        plt.plot(xlist,np.log10(R),linewidth=10,color='orange',zorder=0)

        plt.xlim(xlow,xupp)

        plt.ylim(np.log10(ylow),np.log10(yupp))

        AddRects(ax)

        #plt.xticks(xticks)

        #plt.grid()
    
        #plt.xticks(fontsize=30)
        plt.yticks(fontsize=50)

        ax.set_yticks([0,-2,-4,-6])
        ax.set_yticklabels(
            [r'$0$',r'$-2$',r'$-4$',r'$-6$'])


        #plt.axis("off")

        #ax.spines.right.set_visible(False)
        #ax.spines.top.set_visible(False)
        #ax.spines.bottom.set_visible(False)

        """
        plt.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            left = False,
            right = False, 
            labelbottom=False) # labels along the bottom edge are off
        """
        ax.axes.xaxis.set_visible(False)
        #ax.axes.yaxis.set_visible(False)

        #plt.figure().set_size_inches(8.000000/2.54, 4.500000/2.54, forward=True)#8.200000/2.54, 4.500000/2.54, forward=True)
        #plt.figure().axes[0].set_position([0.25, 0.24, 0.70, 0.7])#[0.18, 0.24, 0.8, 0.7])


        plt.tight_layout()


    #Initial Dist
    fig,ax = plt.subplots()
    #plt.semilogy(xlist,Initlist,label='Initial S')
    #plt.semilogy(xlist,rlist,label='Initial R')    
    Setup(plt,ax,xlist,Initlist,rlist)
    #plt.savefig(str(args.directory) + "1.pdf",bbox_inches='tight')
    plt.savefig(str(args.directory) + "L_%0.1f_1.eps"%(L),bbox_inches='tight')
    plt.close()


    #PAP
    fig,ax = plt.subplots()
    #plt.semilogy(xlist,PAPylist,label='PAP S')
    #plt.semilogy(xlist,rlist,label='Initial R')
    PAPylist[PAPylist<ylow/10] = ylow/10
    Setup(plt,ax,xlist,PAPylist,rlist)
    #plt.savefig(str(args.directory) + "2.pdf",bbox_inches='tight')
    plt.savefig(str(args.directory) + "L_%0.1f_2.eps"%(L),bbox_inches='tight')
    plt.close()

    #Pre Breeding
    fig,ax = plt.subplots()
    #plt.semilogy(xlist,Yearylist,label='Pre Breed S')
    #plt.semilogy(xlist,rlist,label='Initial R')
    Setup(plt,ax,xlist,Yearylist,rlist)
    plt.plot(xlist,np.log10(Yearylist),linewidth=10,color='blue',zorder=0)
    #plt.savefig(str(args.directory) + "3.pdf",bbox_inches='tight')
    plt.savefig(str(args.directory) + "L_%0.1f_3.eps"%(L),bbox_inches='tight')
    plt.close()

    #Post Breeding
    fig,ax = plt.subplots()
    #plt.semilogy(xlist,Yearylist/(rlist+Yearylist),label='Post Breed S')
    #plt.semilogy(xlist,rlist/(rlist+Yearylist),label='Post Breed R')
    Setup(plt,ax,xlist,Yearylist/(rlist+Yearylist),rlist/(rlist+Yearylist))
    #plt.plot(xlist,np.log10(Yearylist/(rlist+Yearylist)),linewidth=10,color='blue',zorder=0)
    #plt.savefig(str(args.directory) + "4.pdf",bbox_inches='tight')
    plt.savefig(str(args.directory) + "L_%0.1f_4.eps"%(L),bbox_inches='tight')
    plt.close()


    #Difference in R
    fig,ax = plt.subplots()
    plt.plot(xlist,np.log10(rlist/(rlist+Yearylist)-rlist),color='orange',linewidth=10)
    plt.fill_between(xlist,np.log10(rlist/(rlist+Yearylist)-rlist),-10,facecolor='#FFD580')
    plt.xlim(xlow,xupp)
    plt.ylim(np.log10(ylow),np.log10(yupp))
    AddRects(ax)
    plt.xticks(xticks)
    #plt.grid()
    plt.xticks(fontsize=30,fontname='Arial')
    plt.yticks(fontsize=30,fontname='Arial')
    ax.axes.xaxis.set_visible(False)
    plt.text(L/2,np.log10(10**-8),r"$\Delta R \left( %d \right)$"%(L),fontsize=30,ha='center')
    plt.tight_layout() 
    plt.savefig(str(args.directory) + "L_%0.1f_RDiff.png"%(L))
    plt.savefig(str(args.directory) + "L_%0.1f_RDiff.eps"%(L))
    plt.close()

