import numpy as np
import matplotlib.pyplot as plt

from scipy import integrate


import os
import time

import copy

starttime = time.time()
##############################################################################
##############################################################################
##############################################################################
#For the very first step
def A(n,SRatio):
    if n%2 == 0:
        return 0
    else:
        return (1-SRatio)*4/(np.pi * n)

def u1(x,t,l,k,Phi):

    if  (x < l):
        tot = 0
        for n in range(1,100):
            tot += A(n,Phi) * np.sin(n*np.pi*x/l) * np.exp(-k*t*(n*np.pi/l)**2)
        return tot
    else:
        return 0
##############################################################################
##############################################################################
##############################################################################
#For the second step
def A0(L,init,xlist):
    return (1/(2*L))*integrate.simps(init, xlist)

def An(n,init,xlist,L):
    tempinit = np.asarray(copy.copy(init))

    for x in range(len(xlist)):
        tempinit[x] *= np.cos(n*np.pi*xlist[x]/L)

    return (1/L)*integrate.simps(tempinit, xlist)

def Bn(n,init,xlist,L):
    tempinit = np.asarray(copy.copy(init))

    for x in range(len(xlist)):
        tempinit[x] *= np.sin(n*np.pi*xlist[x]/L)

    return (1/L)*integrate.simps(tempinit, xlist)

def u2(t,L,k,init,xlist):

    A0Val = A0(L,init,xlist)
    AnList = []
    BnList = []

    print("A0",A0Val)

    for n in range(1,200):
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
k = 0.1#1
Year = 0.1#200
PAP = (Year/2)
T = Year - PAP

Phi = 1e-5


k = 2e-5#1
dx = 1e-4#1
dt = 1e-4#1

SystemSizeList = np.arange(0.01,0.5,0.1)#np.arange(1,100,2)
OSR_ProportionList = np.arange(0.05,1,0.05)




#SaveFile
SaveDirName = ("Saved_Plots/" +
    "PAP_%0.3f_Year_%3f_Characteristic_%0.3f_dx_"%(PAP,Year,k) +
    '{:.1E}'.format(dx)+ "_dt_"+'{:.1E}'.format(dt) + "_InitialRRatio_" +
    '{:.1E}'.format(Phi))

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

##############################################################################
##############################################################################
##############################################################################
MinimumList = []
for OSR_Proportion in OSR_ProportionList:
    print("OSR Proportion:",OSR_Proportion)

    Refuge_Proportion = 1- OSR_Proportion
    R_Ratio_List = []

    for SystemSize in SystemSizeList:
        print("System Size:",SystemSize)
        OSRWidth = (1-Refuge_Proportion)*SystemSize
        RefugeWidth = Refuge_Proportion*SystemSize



        #Step 1
        #Solve the analytical result for the system of Refuge with Dirichlet
        #Boundary Conditions
        xlist = np.arange(int(SystemSize/dx))*dx
        ylist = []
        for x in range(len(xlist)):
            ylist.append(u1(xlist[x],PAP,RefugeWidth,k,Phi))

        plt.figure()
        plt.plot(xlist,ylist,label='1st Step')



        L = SystemSize/2
        print("L",L)
        yprimelist = u2(T,L,k,ylist,xlist)

        plt.plot(xlist,yprimelist,label='2nd Step')

        plt.xlabel("Xpos")
        plt.ylabel("S Population")
        plt.legend(loc='upper right')
        plt.grid()
        plt.savefig(SaveDirName + ("/OSR_Proportion_%0.3f_SystemSize_%0.3f.png"%
            (OSR_Proportion,SystemSize)))
        plt.close()




