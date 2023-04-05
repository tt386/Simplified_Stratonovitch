from Params import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import scipy.special as special

import copy
import sys

import time

starttime = time.time()


#################################
###Argparse
#################################
from argparse import ArgumentParser
parser = ArgumentParser(description='Dimension Specifier')
parser.add_argument('-N','--Dim',type=float,required=True,help='Dimensions')
args = parser.parse_args()

N = int(args.Dim)

NumSaveDirName = (SaveDirName +
    "/Dimension_%d"%(N))

if not os.path.isdir(NumSaveDirName):
    os.mkdir(NumSaveDirName)
    print("Created Directory for Dimension",N)



#########################################
###Main Process##########################
#########################################

#Data Storage
DimensionList = []


IntegralList = []

PAPSystemMatrix = []
EndOfYearSystemMatrix = []

#Adjust scale to make sure I encompass linearly growing minimum (in space)
etalist = etalist#/(N**2)

for eta in etalist:
    print("Dimension:",N,"Eta:",eta)
    dt = min(dr**2 / ((N+1)*eta*FACTOR),1e-4)
    timelist = np.arange(0,1,dt)
    
    Dimension = int(1/dt)#100
    OSRHalfWidth = int(0.5/dr)
    
    
    
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


        System[1:L-1] += (eta*dt/dr**2) * (
            (N-1)*(np.roll(System,-1) - System)[1:L-1]/jlist[1:L-1] +
            (np.roll(System,-1) - 2*System + np.roll(System,1))[1:L-1])


        """
        System[1:L-1] += (eta*dt/dr**2) * (
            np.roll(System,-1) * ((N-1)/jlist + 1) +
            System * (-(N-1)/jlist - 2) +
            np.roll(System,1)
            )[1:L-1]
        """
        System[0] += (eta*dt/dr**2) * (2*System[1] - 2*System[0])

        """
        prod = System * jlist

        System[1:L-1] += (eta*dt/dr**2) * (
            np.roll(prod,-1) - 2*prod + np.roll(prod,1))[1:L-1]/jlist[1:L-1]

        System[0] += (eta*dt/dr**2) * (2*prod[1] - 2*prod[0])
        """
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
        Integrand.append(2*(np.pi**(N/2)) / special.gamma(N/2) *
        (i*dr)**(N-1) * RatioList[i])
        #Integrand.append(4*np.pi * (i*dr)**2 * RatioList[i])
        
    print("Eta:",eta,"NextYearRatio:",simps(Integrand,np.arange(int(Dimension))*dr))
    IntegralList.append(simps(Integrand,np.arange(int(Dimension))*dr))



endtime = time.time()
timetaken = endtime-starttime

print("Time:",endtime-starttime)

OutputDatafilename = NumSaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    N=N,
    PAP=PAP,
    Phi=Phi,
    dr=dr,
    DimensionList=DimensionList,
    etalist=etalist,
    IntegralList=IntegralList,
    PAPSystemMatrix=PAPSystemMatrix,
    EndOfYearSystemMatrix=EndOfYearSystemMatrix,
    timetaken=timetaken)

