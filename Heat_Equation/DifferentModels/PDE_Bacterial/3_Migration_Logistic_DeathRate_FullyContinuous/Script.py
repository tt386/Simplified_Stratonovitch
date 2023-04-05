import numpy as np
import copy
import matplotlib.pyplot as plt
import time
from Params import *


def StepFunc(x,k):
    """
    Approximation for the Heaviside Step Fnction.
    Used to model continuous cutoffs in time and space
    """
    return 1/(1+np.exp(-2*k*x))

from argparse import ArgumentParser
parser = ArgumentParser(description='L Specifier')
parser.add_argument('-L','--Length',type=float,required=True,help='Length of Cropping region')
args = parser.parse_args()

L = float(args.Length)

LSaveDirName = (SaveDirName +
    "/L_%0.3f"%(L))

if not os.path.isdir(LSaveDirName):
    os.mkdir(LSaveDirName)
    print("Created Directory for Length",L)

print("Beginning L:",L)

###################################

starttime = time.time()

#Initial Conditions
RList = np.ones(len(XList)) * (phi)
SList = 1-RList

"""
#Initial Pesticide:
SList[abs(XList)<L] = 0
"""

SPAPDist = []
RPAPDist = []


#This is time-consuming to create, so make it here
XHeaviside = 1-StepFunc(abs(XList)-L/2,kx)

#Loop
for t in TList:
    tempS = copy.copy(SList)
    tempR = copy.copy(RList)

    #Diffusion:
    tempS += dt/(dx**2) * (np.roll(tempS,1) + np.roll(tempS,-1) - 2*tempS)
    tempR += dt/(dx**2) * (np.roll(tempR,1) + np.roll(tempR,-1) - 2*tempR)    

    
    #Breeding
    tot = SList + RList
    tempS += r*dt*SList*(1-tot)
    tempR += r*dt*RList*(1-tot)
    
    """
    #Pesticide
    if t < P:
        tempS[abs(XList) < L] = 0
    """

    
    #Pesticide Application
    if t < P:
        #Ignore dt: assume it's incorporated into dS
        #tempS[abs(XList) < L] -= SList[abs(XList) < L]*dS
        #tempR[abs(XList) < L] -= RList[abs(XList) < L]*dR
        SPAPDist = copy.copy(tempS)
        RPAPDist = copy.copy(tempR)
    

    #Ignore dt: assume it's incorporated into dS
    tempS -= dS*SList*XHeaviside*(1-StepFunc(t-P,kt))
    tempR -= dR*RList*XHeaviside*(1-StepFunc(t-P,kt))

    SList = tempS
    RList = tempR
    
    """
    plt.figure()
    plt.semilogy(XList,SList)
    plt.semilogy(XList,RList)
    plt.xlim(-10,10)
    plt.savefig("Breeding_t_%0.5f.png"%(t))
    plt.close()
    """

DeltaRList = RList - phi

DeltaR = np.sum(DeltaRList) * dx

print("L value: %0.3f, DeltaR: %0.3f, DeltaR/L: %0.5f"%(L,DeltaR,DeltaR/L))

endtime = time.time()
timetaken = endtime-starttime

OutputDatafilename = LSaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    P=P,
    phi=phi,
    r=r,
    dR=dR,
    dS=dS,
    kx=kx,
    kt=kt,
    L=L,
    DeltaR = DeltaR,
    RList =RList,
    SList = SList,
    XList = XList,
    SPAPDist=SPAPDist,
    RPAPDist=RPAPDist,
    timetaken=timetaken)
