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

def u1(x,t,l,k,Phi,N):
    """Numerical result for Dirichlet Heat Eqn"""
    if  (x < l):
        tot = 0
        for n in range(1,N):
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

def u2(t,L,k,init,xlist,N):
    """Numerical result for Periodic Heat Eqn"""
    A0Val = A0(L,init,xlist)
    AnList = []
    BnList = []

    #print("A0",A0Val)

    for n in range(1,N):
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
#SlopeMatrix: matrix of slope values

AnalyticalList = []

NumericalList = []

SlopeList = np.zeros(len(SystemSizeList))
#RawMatrix = np.zeros((len(OSRRatioList),len(SystemSizeList)))
##############################################################################
##############################################################################
##############################################################################

T = Year*(1.-PAP)

Refuge_Proportion = 1- OSR_Proportion
R_Ratio_List = []

slopelist = []

for j in range(len(SystemSizeList)):

    #########################################################################
    #########################################################################
    #########################################################################
    #Analytical


    SystemSize = SystemSizeList[j]
    print("(Year,System Size) =",Year,SystemSize)

    #Create system
    OSRWidth = (1-Refuge_Proportion)*SystemSize
    RefugeWidth = Refuge_Proportion*SystemSize
    xlist = np.arange(int(SystemSize/dx))*dx
    ylist = []


    #Step 1
    #Solve the analytical result for the system of Refuge with Dirichlet
    #Boundary Conditions
    for x in range(len(xlist)):
        ylist.append(u1(xlist[x],PAP*Year,RefugeWidth,k,Phi,SumTerms1))
        #if xlist[x] < RefugeWidth: ylist.append((1-Phi))
        #else: ylist.append(0)
    #Find area, similar to the number of susceptible pests
    step1area = integrate.simps(ylist,xlist)


    #Steo 2
    #Solve the analytical result for the system fo Refuge, OSR and 
    #Periodic Boundary conditions
    L = SystemSize/2
    yprimelist = np.asarray(u2(T,L,k,ylist,xlist,SumTerms2))

    #Ensure no negative values
    yprimelist[yprimelist<0] = 0

    AnalyticalList.append([copy.deepcopy(ylist),copy.deepcopy(yprimelist)])

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
    ###Sim
    days = int(Year/dt)
    MigrationsPerDay = 1

    F = k * (dt)/(dx)**2

    systemsites = int(SystemSize/dx)

    OSRWidth = int(systemsites * OSR_Proportion)
    RefugeSize = int(systemsites * (1-OSR_Proportion))

    systemsites = OSRWidth + RefugeSize

    #Making lists
    S_Domain = np.ones(systemsites)*(1-Phi)
    R_Domain = np.ones(systemsites)*(Phi)

    PesticideDomain = np.ones(OSRWidth)
    PesticideDomain = np.pad(PesticideDomain,(0,RefugeSize),mode='constant')

    PesticideDomain = np.flip(PesticideDomain)

    EndOfPAPDomain = 0

    ####Process#######################
    for t in range(days):
        #print(t)
        if t < PAP*days:
            S_Domain[PesticideDomain == 1] = 0
            EndOfPAPDomain = copy.copy(S_Domain)
        for tm in range(MigrationsPerDay):

            S_Domain = (S_Domain +
                        F*(np.roll(S_Domain,-1) - 2 * S_Domain +
                           np.roll(S_Domain,1)))
        if t < PAP*days:
            EndOfPAPDomain = copy.copy(S_Domain)

    NumericalList.append([copy.deepcopy(EndOfPAPDomain),copy.deepcopy(S_Domain)])


endtime = time.time()

print("Year",Year,"Time taken:",endtime-starttime)

timetaken = endtime-starttime

##############################################
###Output files###############################
##############################################
OutputDatafilename = SaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    PAP=PAP,
    Year=Year,
    Phi=Phi,
    k=k,
    OSR_Proportion=OSR_Proportion,
    SlopeList=SlopeList,
    SystemSizeList=SystemSizeList,
    timetaken=timetaken,
    AnalyticalList=AnalyticalList,
    NumericalList=NumericalList,
    dx=dx,
    dt=dt)

