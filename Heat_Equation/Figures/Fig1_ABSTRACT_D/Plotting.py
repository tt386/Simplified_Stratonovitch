import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as patches
plt.rcParams['text.usetex'] = True

import numpy as np
import time

from scipy.stats import linregress

from scipy.optimize import curve_fit

import sys

import scipy.integrate as integrate
import scipy.special as special

from scipy.signal import argrelextrema

from scipy.signal import savgol_filter

from matplotlib.transforms import Bbox

starttime = time.time()
################
###Find Minima##
################


InfiniteWidth = 3

IsolatedWidth = 2

MultipleWidth = 1
MultipleSep = 1

PeriodicWidth = 8
Periodic_w = 0.4
w_width = PeriodicWidth*Periodic_w
c_width = PeriodicWidth*(1-Periodic_w)


Height = 0.5
totalwidth = 10



def AddRects(ax,SysType):
    if SysType == 'Isolated':
        #Left Infinite Wild
        # Create a Rectangle patch
        rect = patches.Rectangle((0, 0), InfiniteWidth, Height, linewidth=1, edgecolor='k', facecolor='k')
        # Add the patch to the Axes
        ax.add_patch(rect)


        #Isolated Cropping region
        # Create a Rectangle patch
        rect = patches.Rectangle((InfiniteWidth, 0), IsolatedWidth, Height, linewidth=1, edgecolor='k', facecolor='White',hatch='//')
        # Add the patch to the Axes
        ax.add_patch(rect)


        #Right Infinite Wild
        # Create a Rectangle patch
        rect = patches.Rectangle((InfiniteWidth+IsolatedWidth, 0), InfiniteWidth, Height, linewidth=1, edgecolor='k', facecolor='k')
        # Add the patch to the Axes
        ax.add_patch(rect)


        ########################
        ###Annotations##########
        ########################
        plt.annotate(s='', xy=(InfiniteWidth+IsolatedWidth,2*Height), xytext=(InfiniteWidth,2*Height),
            arrowprops=dict(arrowstyle='<->,head_width=1, head_length=0.5',linewidth=4))

        plt.annotate(r'$L$',[InfiniteWidth+0.5*IsolatedWidth,3*Height],fontsize=28,ha='center')

        
        #Rightmost infinite Arrow
        plt.annotate(s='',xy=(2*InfiniteWidth+IsolatedWidth,2*Height),xytext=(1.5*InfiniteWidth+IsolatedWidth,2*Height),
            arrowprops=dict(arrowstyle='->,head_width=1, head_length=0.5',linewidth=4))

        plt.annotate(r'$\infty$',[1.75*InfiniteWidth+IsolatedWidth,3*Height],fontsize=28,ha='center')

        #Leftmost infinite arrow
        plt.annotate(s='',xy=(0.5*InfiniteWidth,2*Height),xytext=(0,2*Height),
            arrowprops=dict(arrowstyle='<-,head_width=1, head_length=0.5',linewidth=4))

        plt.annotate(r'$-\infty$',[0.25*InfiniteWidth,3*Height],fontsize=28,ha='center')



    if SysType == 'Multiple':
        #Left Infinite Wild
        # Create a Rectangle patch
        rect = patches.Rectangle((0, 0), InfiniteWidth, Height, linewidth=1, edgecolor='k', facecolor='k')
        # Add the patch to the Axes
        ax.add_patch(rect)

        ###############################
        #First Multiple Cropping
        # Create a Rectangle patch
        rect = patches.Rectangle((InfiniteWidth, 0), MultipleWidth, Height, linewidth=1, edgecolor='k', facecolor='White',hatch='//')
        # Add the patch to the Axes
        ax.add_patch(rect)

        #First Separation Wild
        # Create a Rectangle patch
        rect = patches.Rectangle((InfiniteWidth + MultipleWidth, 0), MultipleSep, Height, linewidth=1, edgecolor='k', facecolor='k')
        # Add the patch to the Axes
        ax.add_patch(rect)

        ###############################
        #Second Multiple Cropping
        # Create a Rectangle patch
        rect = patches.Rectangle((InfiniteWidth + MultipleWidth + MultipleSep, 0), MultipleWidth, Height, linewidth=1, edgecolor='k', facecolor='White',hatch='//')
        # Add the patch to the Axes
        ax.add_patch(rect)

  
        #Second Separation Wild
        # Create a Rectangle patch
        rect = patches.Rectangle((InfiniteWidth + 2*MultipleWidth + MultipleSep, 0), MultipleSep, Height, linewidth=1, edgecolor='k', facecolor='k')
        # Add the patch to the Axes
        ax.add_patch(rect)
        #################################

        #Third Multiple Cropping
        # Create a Rectangle patch
        rect = patches.Rectangle((InfiniteWidth + 2*MultipleWidth + 2*MultipleSep, 0), MultipleWidth, Height, linewidth=1, edgecolor='k', facecolor='White',hatch='//')
        # Add the patch to the Axes
        ax.add_patch(rect)
        ################################


        #Right infinite Wild
        # Create a Rectangle patch
        rect = patches.Rectangle((InfiniteWidth + 3*MultipleWidth + 2*MultipleSep, 0), InfiniteWidth, Height, linewidth=1, edgecolor='k', facecolor='k')
        # Add the patch to the Axes
        ax.add_patch(rect)


        ########################
        ###Annotations##########
        ########################
        #Width of cropping region
        plt.annotate(s='', xy=(InfiniteWidth+MultipleWidth,2*Height), xytext=(InfiniteWidth,2*Height),
            arrowprops=dict(arrowstyle='<->,head_width=1, head_length=0.5',linewidth=4))

        plt.annotate(r'$L$',[InfiniteWidth+0.5*MultipleWidth,3*Height],fontsize=28,ha='center')


        #Separation Width
        plt.annotate(s='', xy=(InfiniteWidth+2*MultipleWidth+2*MultipleSep,2*Height), xytext=(InfiniteWidth+2*MultipleWidth+MultipleSep,2*Height),
            arrowprops=dict(arrowstyle='<->,head_width=1, head_length=0.5',linewidth=4))

        plt.annotate(r'$\delta$',[InfiniteWidth+2*MultipleWidth+1.5*MultipleSep,3*Height],fontsize=28,ha='center')

        #Rightmost infinite Arrow
        plt.annotate(s='',xy=(2*InfiniteWidth+3*MultipleWidth + 2*MultipleSep,2*Height),xytext=(1.5*InfiniteWidth+3*MultipleWidth + 2*MultipleSep,2*Height),
            arrowprops=dict(arrowstyle='->,head_width=1, head_length=0.5',linewidth=4))

        plt.annotate(r'$\infty$',[1.75*InfiniteWidth+3*MultipleWidth + 2*MultipleSep,3*Height],fontsize=28,ha='center')

        #Leftmost infinite arrow
        plt.annotate(s='',xy=(0.5*InfiniteWidth,2*Height),xytext=(0,2*Height),
            arrowprops=dict(arrowstyle='<-,head_width=1, head_length=0.5',linewidth=4))

        plt.annotate(r'$-\infty$',[0.25*InfiniteWidth,3*Height],fontsize=28,ha='center')



    if SysType == 'Periodic':
        #Left Infinite Wild
        # Create a Rectangle patch
        rect = patches.Rectangle((0, 0), w_width, Height, linewidth=1, edgecolor='k', facecolor='k')
        # Add the patch to the Axes
        ax.add_patch(rect)

        ###############################
        #First Multiple Cropping
        # Create a Rectangle patch
        rect = patches.Rectangle((w_width, 0), c_width, Height, linewidth=1, edgecolor='k', facecolor='White',hatch='//')
        # Add the patch to the Axes
        ax.add_patch(rect)
        

        ##################
        ###Annotate#######
        ##################
        #Overall System Width
        plt.annotate(s='', xy=(PeriodicWidth,2*Height), xytext=(0,2*Height),
            arrowprops=dict(arrowstyle='<->,head_width=1, head_length=0.5',linewidth=4))

        plt.annotate(r'$l$',[0.5*PeriodicWidth,3*Height],fontsize=28,ha='center')

    
        #Wild width
        
        ax.annotate(s='', xy=(w_width,-0.5*Height), xytext=(0,-0.5*Height),
            arrowprops=dict(arrowstyle='<->,head_width=1, head_length=0.5',linewidth=4))

        ax.annotate(r'$l\omega$',[w_width*0.5,-1.5*Height],fontsize=28,ha='center')
        
def Setup(plt,ax,fig,SysType):
    bbox = 0

    if SysType == 'Isolated':
        plt.xlim(0,2*InfiniteWidth+IsolatedWidth)

        x0,x1 = -0.5, 2*InfiniteWidth+IsolatedWidth + 0.5
        y0,y1 = 0, 2*Height

        bbox = Bbox([[x0,y0],[x1,y1]])
        bbox = bbox.transformed(ax.transData).transformed(fig.dpi_scale_trans.inverted())

    if SysType == 'Multiple':
        plt.xlim(0,2*InfiniteWidth + 3*MultipleWidth + 2*MultipleSep)

        x0,x1 = -0.5, 2*InfiniteWidth + 3*MultipleWidth + 2*MultipleSep +0.5
        y0,y1 = 0, 2*Height

        bbox = Bbox([[x0,y0],[x1,y1]])
        bbox = bbox.transformed(ax.transData).transformed(fig.dpi_scale_trans.inverted())


    if SysType == 'Periodic':
        plt.xlim(0,PeriodicWidth)

        x0,x1 = -0.5, PeriodicWidth + 0.5
        y0,y1 = 0, 2*Height

        bbox = Bbox([[x0,y0],[x1,y1]])
        bbox = bbox.transformed(ax.transData).transformed(fig.dpi_scale_trans.inverted())

    plt.ylim(0,2*Height)

    AddRects(ax,SysType)

    #plt.xticks(xticks)

    #plt.grid()

    #plt.xticks(fontsize=30)
    #plt.yticks(fontsize=30,rotation=45)

    plt.axis("off")

    plt.axis("equal")

    plt.tight_layout()

    #fig.set_size_inches(18.000000/2.54, 4.500000/2.54, forward=True)

    return bbox


fig,ax = plt.subplots()
bbox = Setup(plt,ax,fig,'Isolated')
plt.savefig("Isolated.png",bbox_inches=bbox,transparent=True)
plt.close()


fig,ax = plt.subplots()
bbox = Setup(plt,ax,fig,'Multiple')
plt.savefig("Multiple.png",bbox_inches=bbox,transparent=True)
plt.close()


fig,ax = plt.subplots()
bbox = Setup(plt,ax,fig,'Periodic')
plt.savefig("Periodic.png",bbox_inches=bbox,transparent=True)
plt.close()
