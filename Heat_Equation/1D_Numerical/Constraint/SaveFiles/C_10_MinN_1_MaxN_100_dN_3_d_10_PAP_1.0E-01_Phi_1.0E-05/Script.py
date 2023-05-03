from Params import *

import numpy as np
from scipy import integrate

import copy
import sys

import time

from matplotlib import pyplot as plt


starttime = time.time()



#################################
###Argparse
#################################
from argparse import ArgumentParser
parser = ArgumentParser(description='Number of Regions')
parser.add_argument('-N','--Num',type=float,required=True,help='Number of L')
args = parser.parse_args()

N = int(args.Num)

NumSaveDirName = (SaveDirName +
    "/N_%d"%(N))

if not os.path.isdir(NumSaveDirName):
    os.mkdir(NumSaveDirName)
    print("Created Directory for Number",N)


#########################################
###Main Process##########################
#########################################

#Length of the Sub-units
L = C/N


#SusceptiblePests
SNum = np.ones(len(xlist)) * (1-Phi)

#Pesticide Mask
PMask = np.ones(len(SNum))

#Width of periodic subunit
W = (L+d)/dx

PMask[xlist%W < L/dx] = 0
PMask[xlist >= N*W] = 1



PAPSNum = copy.copy(SNum)

t = 0
while t < 1:
    SNum += (dt/dx**2) * (np.roll(SNum,1) + np.roll(SNum,-1) - 2*SNum)

    if t < PAP:
        SNum *= PMask
        PAPSNum = copy.copy(SNum)

    t += dt


ENDSNum = copy.copy(SNum)

"""
plt.figure()
plt.plot(xlist*dx,PAPSNum)
plt.savefig("PAPstate.png")

plt.figure()
plt.plot(xlist*dx,ENDSNum)
plt.savefig("Endstate.png")
"""

#Breeding
RNum = Phi/(SNum + Phi)

#Change
dR = RNum - Phi

ApproxIntegral = integrate.simps(dR,xlist)

print("N:",N,", Integral", ApproxIntegral)

endtime = time.time()
timetaken = endtime-starttime
print("Time Taken:",timetaken)


OutputDatafilename = NumSaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    N=N,
    PAP=PAP,
    Phi=Phi,
    C=C,
    d=d,
    dx=dx,
    dt=dt,
    PAPSNum=PAPSNum,
    ENDSNum=ENDSNum,
    xlist=xlist,
    dR=dR,
    timetaken=timetaken)

