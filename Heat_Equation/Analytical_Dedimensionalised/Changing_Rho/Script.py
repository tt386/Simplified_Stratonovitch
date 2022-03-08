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
#Function for PAP Dist
def A(n,SRatio):
    """Returns the Fourier Coefficients for the Dicichlet Heat Eqn"""
    if n%2 == 0:
        return 0
    else:
        return (1-SRatio)*4/(np.pi * n)

def u1(x,eta,rho,PAP,Phi,N):
    """Numerical result for Dirichlet Heat Eqn"""
    if  (x < rho):
        tot = 0
        for n in range(1,N):
            tot += A(n,Phi) * np.sin(n*np.pi*x/rho) * np.exp(-eta*PAP*(n*np.pi/rho)**2)
        return tot
    else:
        return 0
#####################################################
#####################################################
#####################################################
def Val(yprimelist):
    InitialCondition = np.ones(int(1/dx))*Phi
    InitialTot = np.sum(InitialCondition)#MAYBE I SHOULD INTEGRATE

    RatioAfterBreeding = InitialCondition/(InitialCondition+yprimelist)
    np.set_printoptions(threshold=sys.maxsize)
    #print(RatioAfterBreeding)

    FinalTot =np.sum(RatioAfterBreeding)
    
    return FinalTot    
    
##########################################################
###For the second step####################################
##########################################################
def A0(init,xlist):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    return integrate.simps(init, xlist)

def An(n,init,xlist):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    tempinit = np.asarray(copy.copy(init))

    for x in range(len(xlist)):
        tempinit[x] *= np.cos(2*n*np.pi*xlist[x])

    return 2*integrate.simps(tempinit, xlist)

def Bn(n,init,xlist):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    tempinit = np.asarray(copy.copy(init))

    for x in range(len(xlist)):
        tempinit[x] *= np.sin(2*n*np.pi*xlist[x])

    return 2*integrate.simps(tempinit, xlist)

def u2(eta,PAP,init,xlist,N):
    """Numerical result for Periodic Heat Eqn"""
    A0Val = A0(init,xlist)
    AnList = []
    BnList = []

    #print("A0",A0Val)

    for n in range(1,N):
        AnList.append(An(n,init,xlist))
        BnList.append(Bn(n,init,xlist))

    

    resultslist = []
    
    #ratiolist = []
    for x in xlist:
        tot = A0Val

        for n in range(1,len(AnList)):
            tot += (np.exp(-eta*(1-PAP)*(2*np.pi*n)**2) *
                (AnList[n-1] * np.cos(2*n*np.pi*x) +
                BnList[n-1] * np.sin(2*n*np.pi*x)))

        resultslist.append(max(tot,0))

    #Val = Val(resultslist)
    
    return resultslist#, Val

##############################################################################
##############################################################################
##############################################################################
from argparse import ArgumentParser
parser = ArgumentParser(description='Different OSR Values')
parser.add_argument('-r','--rho',type=float,required=True,help='rho')
args = parser.parse_args()

rho = args.rho


rhoSaveDirName = SaveDirName + "/rho_" + '{:.1E}'.format(rho)

if not os.path.isdir(rhoSaveDirName):
    os.mkdir(rhoSaveDirName)
    print("Created Directory for rho",rho)
 



#SlopeMatrix: matrix of slope values
SlopeList = np.zeros(len(etaList))
#RawMatrix = np.zeros((len(OSRRatioList),len(SystemSizeList)))
##############################################################################
##############################################################################
##############################################################################

R_Ratio_List = []
minlist = []
for j in range(len(etaList)):
    SystemSize = 1.
    eta = etaList[j]
    print("(rho,eta) =",rho,eta)

    #Create system
    OSRWidth = (1-rho)*SystemSize
    RefugeWidth = rho*SystemSize
    xlist = np.arange(int(SystemSize/dx))*dx
    ylist = []


    #Step 1
    #Solve the analytical result for the system of Refuge with Dirichlet
    #Boundary Conditions
    for x in xlist:
        ylist.append(u1(x,eta,rho,PAP,Phi,N))

    #Find area, similar to the number of susceptible pests
    step1area = integrate.simps(ylist,xlist)


    #Steo 2
    #Solve the analytical result for the system fo Refuge, OSR and 
    #Periodic Boundary conditions
    L = SystemSize/2
    yprimelist = np.asarray(u2(eta,PAP,ylist,xlist,N))

    #Ensure no negative values
    yprimelist[yprimelist<0] = 0

    #Again find area, similar to the number of susceptible pests
    step2area = integrate.simps(yprimelist,xlist)
    
    #print("1st area / 2nd area:",step1area/step2area)

    minlist.append(min(yprimelist))

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

print("rho",rho,"Time taken:",endtime-starttime)

timetaken = endtime-starttime

##############################################
###Output files###############################
##############################################
OutputDatafilename = rhoSaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    PAP=PAP,
    Phi=Phi,
    rho=rho,
    N=N,
    SlopeList=SlopeList,
    minlist=minlist,
    rhoList=rhoList,
    etaList=etaList,
    timetaken=timetaken)
print("rho",rho,"Time taken:",endtime-starttime)
