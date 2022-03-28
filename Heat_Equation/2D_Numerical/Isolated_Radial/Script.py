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
IntegralList = []

PAPSystemMatrix = []
EndOfYearSystemMatrix = []

for eta in etalist:
    print("Eta:",eta)
    dt = dr**2 / (3*eta)#min(dr**2 / (3*eta),1e-4)
    timelist = np.arange(0,1,dt)
    
    Dimension = int(1/dt)#100
    OSRHalfWidth = int(0.5/dr)
    
    print("dt:",dt)
    print("Dimension Size:",Dimension)
    
    
    
    #######################################
    ###System Setup########################
    #######################################
    
    System = np.ones(int(Dimension)) * (1-Phi)#/dr))
    
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
        System[1:L-1] += (eta*dt/dr**2) * (
            (np.roll(System,-1) - System)[1:L-1]/jlist[1:L-1] + 
            (np.roll(System,1) - 2*System + np.roll(System,-1))[1:L-1])
        
        System[0] += (eta*dt/dr**2) * (2*System[1] - 2*System[0])
        
        if timelist[t] < PAP:
            System *= CSystem
        
        if timelist[t] < PAP and (timelist[t+1]>=PAP):
            PAPSystemMatrix.append(copy.deepcopy(System))
            
            """
            plt.semilogx(np.arange(int(Dimension))*dr,System,label='PAP')
            plt.show()
            """

            """
            Integrand = []
            for i in range(len(System)):
                Integrand.append(2*np.pi * (i*dr) * System[i])
            IntegralList.append(simps(Integrand,np.arange(int(Dimension))*dr))
            """
    """
    plt.semilogx(np.arange(int(Dimension))*dr,System,label='End of Year')
    plt.legend(loc='lower right')
    plt.grid()
    plt.show()
    """
    EndOfYearSystemMatrix.append(copy.deepcopy(System))
    """
    Integrand = []
    for i in range(len(System)):
        Integrand.append(2*np.pi * (i*dr) * System[i])
    IntegralList.append(simps(Integrand,np.arange(int(Dimension))*dr))
        
    print(IntegralList)
    """
    
    #Next year ratio
    Integrand = []
    RatioList = Phi/(Phi+System) - Phi
    for i in range(len(System)):
        Integrand.append(2*np.pi * (i*dr) * RatioList[i])
        
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
    Dimension=Dimension,
    etalist=etalist,
    IntegralList=IntegralList,
    PAPSystemMatrix=PAPSystemMatrix,
    EndOfYearSystemMatrix=EndOfYearSystemMatrix,
    timetaken=timetaken)

"""
plt.figure()
plt.loglog(etalist,NextYearList)
plt.grid()
plt.xlabel("ETA")
plt.ylabel("New Resistant")
plt.show()
"""

