from Params import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import copy
import sys

import time

starttime = time.time()

#########################################
###Main Process##########################
#########################################

#Data Storage
DimensionList = []


IntegralList = []

PAPSystemMatrix = []
EndOfYearSystemMatrix = []

for eta in etalist:
    print("Eta:",eta)
    dt = min(dr**2 / (3*eta),1e-4)#min(dr**2 / (3*eta),1e-4)
    timelist = np.arange(0,1,dt)
    
    Dimension = int(1/dt)#100
    OSRHalfWidth = int(0.5/dr)
    
    print("dt:",dt)
    print("Dimension Size:",Dimension)
    
    
    DimensionList.append(Dimension)
    
    #######################################
    ###System Setup########################
    #######################################
    
    System = np.ones(int(Dimension)) * (1-Phi)#/dr))
    RSystem =  np.ones(int(Dimension)) * (Phi)
    """
    for j in range(len(System)):
        System[j] = (1-Phi) * j * (dr)**2
        RSystem[j] = (Phi) * j * (dr)**2
    """
    CSystem = np.ones(int(Dimension))#/dr))
    CSystem[:OSRHalfWidth] = 0
    #print(CSystem)
    
    
    #######################################
    ###Main System#########################
    #######################################
    jlist =np.arange(len(System))
    
    """
    plt.figure()
    plt.title("Eta %0.5f"%(eta))
    """
    for t in range(len(timelist)):
        #print("Time:",timelist[t])
        
        L = len(System)
        prod = System * jlist

        System[1:L-1] += (eta*dt/dr**2) * (
            np.roll(prod,-1) - 2*prod + np.roll(prod,1))[1:L-1]/jlist[1:L-1]

        System[0] += (eta*dt/dr**2) * (2*prod[1] - 2*prod[0])
         
        """
                        (eta*dt/dr**2) * (
            (np.roll(System,-1) - System)[1:L-1]/jlist[1:L-1] + 
            (np.roll(System,1) - 2*System + np.roll(System,-1))[1:L-1])
        
        System[0] += (eta*dt/dr**2) * (2*System[1] - 2*System[0])
        """
        if timelist[t] < PAP:
            System *= CSystem
        
        if timelist[t] < PAP and (timelist[t+1]>=PAP):
            PAPSystemMatrix.append(copy.deepcopy(System))
            
    EndOfYearSystemMatrix.append(copy.deepcopy(System))
    
    #Next year ratio
    Integrand = []
    RatioList = Phi/(Phi+System) - Phi
    for i in range(len(System)):
        Integrand.append(4*np.pi * (i*dr)**2 * RatioList[i])
        
    print("Eta:",eta,"NextYearRatio:",simps(Integrand,np.arange(int(Dimension))*dr))
    IntegralList.append(simps(Integrand,np.arange(int(Dimension))*dr))



endtime = time.time()
timetaken = endtime-starttime

print("Time:",endtime-starttime)

OutputDatafilename = SaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    PAP=PAP,
    Phi=Phi,
    dr=dr,
    DimensionList=DimensionList,
    etalist=etalist,
    IntegralList=IntegralList,
    PAPSystemMatrix=PAPSystemMatrix,
    EndOfYearSystemMatrix=EndOfYearSystemMatrix,
    timetaken=timetaken)

