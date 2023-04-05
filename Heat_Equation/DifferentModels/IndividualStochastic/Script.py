import numpy as np
import copy
import random
import time
import matplotlib.pyplot as plt
from Params import *


#############################
##L Value####################
#############################
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


starttime = time.time()

DeltaRList = []

MeanDeltaRList = []

#for L in LList:
#XList
XList= np.arange(0,alpha*2*M*T+L).astype(int)

for R in range(Repeats):
    print("L, R",L,R)
    #Final Dists
    FinalR = np.zeros(len(XList))
    FinalS = np.zeros(len(XList))
    
    RNum = 0
    for x in XList:
        #print("x,",x)
        for k in range(K):
            xpos = copy.copy(x)
            
            Type = 0
            D = dS
            if random.uniform(0,1) < phi:
                Type = 1
                RNum += 1
                D= dR
                
            Die = False
            for t in TimeList:
                dx = random.randint(-M,M)
                
                xpos += dx
                
                if xpos < 0:
                    xpos += len(XList)
                elif xpos >= len(XList):
                    xpos -= len(XList)
                    
                if t < P:
                    if xpos < L:
                        if random.uniform(0,1)<D:
                            Die = True
                            break
    
            if not Die:
                if Type == 0:
                    FinalS[xpos] += 1
                else:
                    FinalR[xpos] += 1
    endtime = time.time()
    
    timetaken = endtime - starttime
    
    #print("TimeTaken:",timetaken)
           
    """
    #Easy Breeding Algorithm
    denominator = FinalR+FinalS
    denominator[denominator == 0] = 1
    PostBreedR = np.rint(K*FinalR/(denominator)).astype(int)
    """      
    ##############################
    #Stochastic Breeding Algorithm
    ##############################
    PostBreedR = np.zeros(len(XList))
    #print("Start Breeding")
    for x in XList:
        EndSNum = FinalS[x]
        EndRNum = FinalR[x]
        #For each to be generated:
        for k in range(K):
            if random.randint(0,EndSNum+EndRNum-1) < EndRNum:
                PostBreedR[x] += 1
                
    endtime = time.time()
    timetaken = endtime - starttime
    print("TimeTaken:",timetaken)
    print("L",L,"R",R,"RNum",RNum,"FinalR",np.sum(PostBreedR))    

    DeltaR = np.sum(PostBreedR) - RNum
    DeltaRList.append(DeltaR)
    DeltaRL = DeltaR/L

    #MeanDeltaRList.append(np.mean(DeltaRList))

    #print(MeanDeltaRList)

endtime = time.time()
timetaken = endtime-starttime

OutputDatafilename = LSaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    P=P,
    M=M,
    phi=phi,
    dR=dR,
    dS=dS,
    L=L,
    Repeats=Repeats,
    DeltaRList = DeltaRList,
    timetaken=timetaken)
"""
MeanDeltaRList = np.asarray(MeanDeltaRList)
LList = np.asarray(LList)

plt.figure()
plt.plot(LList,MeanDeltaRList/LList)
plt.show()
"""
