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
def A(n,System,xlist):
    return integrate.simps(System * np.sin(n*np.pi*xlist/rho),xlist)



def PAPDist(N,rho,PAP,eta,System,xlist):
    NewSystem = np.zeros(len(System))
    
    for n in range(N):
        NewSystem += (A(n,System,xlist) * 
                    np.sin(n*np.pi*xlist/rho) * 
                    np.exp(-eta*PAP*(n*np.pi/rho)**2))

    NewSystem[xlist>rho] = 0

    return (2/rho) * NewSystem

##############################################################################
##############################################################################
##############################################################################
#Function for End Of Year 
def A0(System,xlist):
    return integrate.simps(System,xlist)

def An(n,System,xlist):
    return 2*integrate.simps(System * np.cos(2*n*np.pi*xlist),xlist)

def Bn(n,System,xlist):
    return 2*integrate.simps(System * np.sin(2*n*np.pi*xlist),xlist)

def EndOfYearDist(N,rho,t,eta,System,xlist):
    NewSystem = np.zeros(len(System))

    NewSystem += A0(System,xlist)

    for n in range(1,N):
        NewSystem += ((An(n,System,xlist) * np.cos(2*n*np.pi*xlist) +
                    Bn(n,System,xlist) * np.sin(2*n*np.pi*xlist)) *
                    np.exp(-eta*t*(2*n*np.pi)**2))

    return NewSystem


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
 



##############################################################################
##############################################################################
##############################################################################
RTotMatrix = []

SystemMatrix = []
PAPSystemMatrix = []
RSystemMatrix = []

for eta in etaList:
    SystemSize = 1.
    print("(rho,eta) =",rho,eta)

    #Create initial system
    xlist = np.arange(0,1,dx)
    System = np.ones(len(xlist)) * (1-Phi)      #Number of Susceptible
    RSystem = 1-System                          #Number of Resistant

    RTotList = []

    RTotList.append(integrate.simps(RSystem,xlist))

    SystemList = []
    PAPSystemList = []
    RSystemList = []
    for Y in range(Years):
        StartOfYearSystem = copy.copy(System)

        #Pesticide Applied:
        System[xlist>rho] = 0

        #Step 1
        #Solve the analytical result for the system of Refuge with Dirichlet
        #Boundary Conditions
        System = PAPDist(N,rho,PAP,eta,System,xlist)

        PAPSystemList.append(copy.deepcopy(System))

        #Steo 2
        #Solve the analytical result for the system fo Refuge, OSR and 
        #Periodic Boundary conditions
        System = EndOfYearDist(N, rho, 1-PAP, eta, System, xlist)

        RSystem = EndOfYearDist(N, rho, 1, eta, RSystem, xlist)

        #Ensure no negative values
        System[System<0] = 0

        #Storing End of Year Distributions
        SystemList.append(copy.deepcopy(System))
        RSystemList.append(copy.deepcopy(RSystem))


        #Breeding
        System = System/(System + RSystem)
        RSystem = 1-System

        #End of year counting
        RTotList.append(integrate.simps(RSystem,xlist)) 

    RTotMatrix.append(RTotList)

    SystemMatrix.append(SystemList)
    PAPSystemMatrix.append(PAPSystemList)
    RSystemMatrix.append(RSystemList)

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
    RTotMatrix=RTotMatrix,
    SystemMatrix=SystemMatrix,
    PAPSystemMatrix=PAPSystemMatrix,
    RSystemMatrix=RSystemMatrix,
    rhoList=rhoList,
    etaList=etaList,
    timetaken=timetaken,
    dx=dx)
print("rho",rho,"Time taken:",endtime-starttime)
