from Params import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

import copy

import time

starttime = time.time()

#################################
###Argparse
#################################
from argparse import ArgumentParser
parser = ArgumentParser(description='Different OSR Separation')
parser.add_argument('-r','--rho',type=float,required=True,help='Sep')
args = parser.parse_args()

rho = args.rho


RhoSaveDirName = (SaveDirName +
    "/OSRrho_" + '{:.3E}'.format(rho))

if not os.path.isdir(RhoSaveDirName):
    os.mkdir(RhoSaveDirName)
    print("Created Directory for Sep",rho)


OSRSize = int(Dimension*(1.-rho))
OSRHalfWidth = int(OSRSize/2)

#################################
###Constructed Arrays############
#################################
OriginalSystem = (1-Phi) * np.ones((Dimension,Dimension),dtype=float)
R_System = Phi * np.ones((Dimension,Dimension),dtype=float)

#Complementary system: multiply during pesticide
CSystem = np.ones((Dimension,Dimension),dtype=float)
for i in range(-OSRHalfWidth,OSRHalfWidth+1):
    for j in range(-OSRHalfWidth,OSRHalfWidth+1):
        if OSR_Shape == "Square":
            CSystem[int(Dimension/2) + j][int(Dimension/2) + i] = 0

        elif OSR_Shape == "Circle":
            if i**2 + j**2 <= OSRHalfWidth**2:
                CSystem[int(Dimension/2) + j][int(Dimension/2) + i] = 0

        else:
            raise ValueError('Specify a Correct OSR geometry') 
#print(CSystem)





#########################################
###Main Process##########################
#########################################

#Data Storage
IntegralList = []

PAPSystemMatrix = []
EndOfYearSystemMatrix = []

for eta in etalist:
    System = copy.deepcopy(OriginalSystem)
    #Constructed Params
    dt = min(h**2 / (4*eta),1e-4)
    timelist = np.arange(0,1,dt)

    for t in range(len(timelist)):
        #print("time:",timelist[t],"TotalPestNum:",np.sum(System))
        System = ((1-4*eta*dt/h**2) * System + 
                  eta*dt/h**2 * 
                  (np.roll(System,1,0) + np.roll(System,-1,0) +
                   np.roll(System,1,1) + np.roll(System,-1,1)))
        
        if timelist[t] < PAP:
            System *= CSystem
        
        if (timelist[t] <= PAP) and (timelist[t+1]>PAP):
            PAPSystemMatrix.append(copy.deepcopy(System))

        """
        if t%(500) == 0:
            plt.figure(t)
            plt.imshow(System, cmap='hot_r', vmin=0, vmax=1,interpolation='nearest')
            plt.colorbar()
            plt.show()
        """
    EndOfYearSystemMatrix.append(copy.deepcopy(System))

    x = np.arange(0,Dimension*h,h)
    y = np.arange(0,Dimension*h,h)

    Integrand = R_System/(R_System+System) - R_System
    Integral = simps(simps(Integrand, y), x)
    IntegralList.append(Integral)
    print("rho:",rho,"eta:",eta,"Integral:",Integral)

    """
    plt.figure(-1)
    plt.imshow(Integrand, cmap='hot_r',interpolation='nearest')
    plt.colorbar()
    plt.show()
    """


endtime = time.time()

timetaken = endtime-starttime

print("Time Taken:",timetaken)

OutputDatafilename = RhoSaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    PAP=PAP,
    Phi=Phi,
    h=h,
    rho=rho,
    Dimension=Dimension,
    etalist=etalist,
    IntegralList=IntegralList,
    PAPSystemMatrix=PAPSystemMatrix,
    EndOfYearSystemMatrix=EndOfYearSystemMatrix,
    timetaken=timetaken)





plt.figure()
plt.loglog(etalist,IntegralList)
plt.grid()
plt.xlabel("Eta")
plt.ylabel("Integral")
plt.savefig(SaveDirName+"/Integral.png")
