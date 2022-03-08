import numpy as np
import matplotlib.pyplot as plt

from scipy import integrate


import os
import time

import copy
import sys

from Params import *

starttime = time.time()
##############################################################################
##############################################################################
##############################################################################
#For the very first step
def A(n,SRatio):
    """Returns the Fourier Coefficients for the Dicichlet Heat Eqn"""
    if n%2 == 0:
        return 0
    else:
        return (1-SRatio)*4/(np.pi * n)

def u1(x,t,l,k,Phi):
    """Numerical result for Dirichlet Heat Eqn"""
    if  (x < l):
        tot = 0
        for n in range(1,100):
            tot += A(n,Phi) * np.sin(n*np.pi*x/l) * np.exp(-k*t*(n*np.pi/l)**2)
        return tot
    else:
        return 0
##############################################################################
##############################################################################
##############################################################################
#For the second step
def A0(L,init,xlist):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    return (1/(2*L))*integrate.simps(init, xlist)

def An(n,init,xlist,L):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    tempinit = np.asarray(copy.copy(init))

    for x in range(len(xlist)):
        tempinit[x] *= np.cos(n*np.pi*xlist[x]/L)

    return (1/L)*integrate.simps(tempinit, xlist)

def Bn(n,init,xlist,L):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    tempinit = np.asarray(copy.copy(init))

    for x in range(len(xlist)):
        tempinit[x] *= np.sin(n*np.pi*xlist[x]/L)

    return (1/L)*integrate.simps(tempinit, xlist)

def u2(t,L,k,init,xlist):
    """Numerical result for Periodic Heat Eqn"""
    A0Val = A0(L,init,xlist)
    AnList = []
    BnList = []

    #print("A0",A0Val)

    for n in range(1,200):
        AnList.append(An(n,init,xlist,L))
        BnList.append(Bn(n,init,xlist,L))

    resultslist = []
    for x in xlist:
        tot = A0Val

        for n in range(1,len(AnList)):
            tot += (np.exp(-k*t*(np.pi*n/L)**2) *
                (AnList[n-1] * np.cos(n*np.pi*x/L) +
                BnList[n-1] * np.sin(n*np.pi*x/L)))
        resultslist.append(tot)
    return resultslist
##############################################################################
##############################################################################
##############################################################################
from argparse import ArgumentParser
parser = ArgumentParser(description='Different OSR Values')
parser.add_argument('-O','--OSRRatio',type=float,required=True,help='OSRRatio')
args = parser.parse_args()

OSR_Proportion = args.OSRRatio


OSRSaveDirName = SaveDirName + "/OSRRatio_" + '{:.1E}'.format(OSR_Proportion)

if not os.path.isdir(OSRSaveDirName):
    os.mkdir(OSRSaveDirName)
    print("Created Directory for OSR Proportion",OSR_Proportion)
 



#SlopeMatrix: matrix of slope values
SlopeList = np.zeros(len(SystemSizeList))
#RawMatrix = np.zeros((len(OSRRatioList),len(SystemSizeList)))
##############################################################################
##############################################################################
##############################################################################

Refuge_Proportion = 1- OSR_Proportion
R_Ratio_List = []

slopelist = []

for j in range(len(SystemSizeList)):
    SystemSize = SystemSizeList[j]
    print("(OSRRatio,System Size) =",OSR_Proportion,SystemSize)

    #Create system
    OSRWidth = (1-Refuge_Proportion)*SystemSize
    RefugeWidth = Refuge_Proportion*SystemSize
    xlist = np.arange(int(SystemSize/dx))*dx
    ylist = []


    #Step 1
    #Solve the analytical result for the system of Refuge with Dirichlet
    #Boundary Conditions
    for x in range(len(xlist)):
        ylist.append(u1(xlist[x],PAP,RefugeWidth,k,Phi))

    #Find area, similar to the number of susceptible pests
    step1area = integrate.simps(ylist,xlist)


    #Steo 2
    #Solve the analytical result for the system fo Refuge, OSR and 
    #Periodic Boundary conditions
    L = SystemSize/2
    yprimelist = np.asarray(u2(T,L,k,ylist,xlist))

    #Ensure no negative values
    yprimelist[yprimelist<0] = 0

    #Again find area, similar to the number of susceptible pests
    step2area = integrate.simps(yprimelist,xlist)
    
    #print("1st area / 2nd area:",step1area/step2area)


    ##########################################################################
    ##########################################################################
    InitialCondition = np.ones(int(SystemSize/dx))*Phi
    InitialTot = np.sum(InitialCondition)#MAYBE I SHOULD INTEGRATE

    RatioAfterBreeding = InitialCondition/(InitialCondition+yprimelist)
    np.set_printoptions(threshold=sys.maxsize)
    #print(RatioAfterBreeding)

    FinalTot =np.sum(RatioAfterBreeding)

    #print(InitialTot)
    #print(FinalTot)

    SlopeList[j] = np.log(FinalTot)- np.log(InitialTot)

    ##########################################################################
    ##########################################################################
    ##########################################################################

endtime = time.time()

print("OSR Proportion",OSR_Proportion,"Time taken:",endtime-starttime)

timetaken = endtime-starttime

##############################################
###Output files###############################
##############################################
OutputDatafilename = OSRSaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    PAP=PAP,
    Year=Year,
    Phi=Phi,
    k=k,
    OSR_Proportion=OSR_Proportion,
    SlopeList=SlopeList,
    #RawMatrix,
    OSRRatioList=OSRRatioList,
    SystemSizeList=SystemSizeList,
    timetaken=timetaken)

