import numpy as np
import matplotlib.pyplot as plt

from scipy import integrate


import os
import time

import copy
import sys

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

    #print("A0",A0Val)

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
Year = 0.1#200#0.1#200
PAP = (Year/2)
T = Year - PAP

Phi = 1e-5 #Initial Conditions


OSR_Proportion = 0.5

#k = 2e-5#1#2e-5#1
dx =1e-4#1 #1e-4#1
dt =1e-4#1 #1e-4#1

SystemSizeList = np.arange(0.001,0.01,0.0001)#np.arange(1,100,2)#np.arange(0.01,0.5,0.1)#np.arange(1,100,2)

CharacteristicList = np.arange(1e-6,2e-5,1e-6)

OSR_ProportionList = np.arange(0.1,0.9,0.05)




#SaveFile
SaveDirName = ("Saved_Plots/" +
    "PAP_%0.3f_Year_%3f_OSR_Proportion_%0.3f"%(PAP,Year,OSR_Proportion) + 
    "_MaxSize_%0.5f_dx_"%(SystemSizeList[-1]) +
    '{:.1E}'.format(dx)+ "_dt_"+'{:.1E}'.format(dt) + "_InitialRRatio_" +
    '{:.1E}'.format(Phi))

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

##############################################################################
##############################################################################
##############################################################################
MinimumList = []
for k in CharacteristicList:
    #print("OSR Proportion:",OSR_Proportion)

    Refuge_Proportion = 1- OSR_Proportion
    R_Ratio_List = []

    slopelist = []

    for SystemSize in SystemSizeList:
        print("(Characteristic,System Size) =",k,SystemSize)
        OSRWidth = (1-Refuge_Proportion)*SystemSize
        RefugeWidth = Refuge_Proportion*SystemSize



        #Step 1
        #Solve the analytical result for the system of Refuge with Dirichlet
        #Boundary Conditions
        xlist = np.arange(int(SystemSize/dx))*dx
        ylist = []
        for x in range(len(xlist)):
            ylist.append(u1(xlist[x],PAP,RefugeWidth,k,Phi))

        step1area = integrate.simps(ylist,xlist)


        plt.figure()
        plt.plot(xlist,ylist,label='1st Step')



        L = SystemSize/2
        #print("L",L)
        yprimelist = np.asarray(u2(T,L,k,ylist,xlist))

        yprimelist[yprimelist<0] = 0

        step2area = integrate.simps(yprimelist,xlist)
        
        print("1st area / 2nd area:",step1area/step2area)


        plt.plot(xlist,yprimelist,label='2nd Step')

        plt.xlabel("Xpos")
        plt.ylabel("S Population")
        plt.legend(loc='upper right')
        plt.grid()
        plt.savefig(SaveDirName + ("/Characteristic_" + '{:.1E}'.format(k) + 
            "_SystemSize_%0.3f.png"%(SystemSize)))
        plt.close()

        ##########################################################################
        ##########################################################################
        InitialCondition = np.ones(int(SystemSize/dx))*Phi
        InitialTot = np.sum(InitialCondition)

        RatioAfterBreeding = InitialCondition/(InitialCondition+yprimelist)
        np.set_printoptions(threshold=sys.maxsize)
        print(RatioAfterBreeding)

        FinalTot =np.sum(RatioAfterBreeding)

        print(InitialTot)
        print(FinalTot)


        slopelist.append(np.log(FinalTot)- np.log(InitialTot))        

        """   
        InitialCondition = np.ones(int(SystemSize/dx))*Phi

        InitialIntegral = integrate.simps(InitialCondition,xlist)

        InitialRatio = InitialIntegral/integrate.simps(np.ones(int(SystemSize/dx)),xlist)

        FinalRatio = InitialIntegral/(InitialIntegral + step2area)

        
        slopelist.append(np.log(FinalRatio) - np.log(InitialRatio))
        """
        
        """
        RRatioList = []
        SRatioList = []
        for i in range(len(ylist)):
            RRatioList.append(Phi/(yprimelist[i]+Phi))
            SRatioList.append(1-RRatioList[-1])
        
        
        TotalRRatio = sum(RRatioList)/SystemSize    
        #print("SystemSize, TotalRRatio",SystemSize,TotalRRatio)
        
        R_Ratio_List.append(TotalRRatio)
        
        logslope = np.log(TotalRRatio)-np.log(Phi)
        #print("Log of slope",logslope)
        slopelist.append(logslope)
        """
        ##########################################################################
        ##########################################################################
        ##########################################################################
            
            
    plt.figure()
    plt.semilogy(SystemSizeList,slopelist)

    plt.xlabel("System Size")
    plt.ylabel("Log of initial yearly R Allele Slope")

    plt.title("Characteristic " + '{:.1E}'.format(k))
    plt.grid()
    plt.savefig(SaveDirName + 
        "/Minima_Curve_Characteristic_" + '{:.1E}'.format(k) + ".png")
    plt.close()
                   


    #########################################################################
    #########################################################################
    #########################################################################
    #Find the minumum of the curve
    minRAllele = 10000000
    minSystemSize = -1
    for i in range(len(slopelist)):
        if slopelist[i] < minRAllele:
            minRAllele = slopelist[i]
            minSystemSize = SystemSizeList[i]

    MinimumList.append(minSystemSize)
    print("minRAllele:",minRAllele)
    print("minSystemSize",minSystemSize)
    """

        if ((R_Ratio_List[i-1] > R_Ratio_List[i]) and 
                (R_Ratio_List[i+1] > R_Ratio_List[i])):
            MinimumList.append(SystemSizeList[i])
            break
    """
    
plt.figure()
plt.plot(CharacteristicList,MinimumList)

plt.xlabel("Characteristic Migration Scale")
plt.ylabel("System Size at Minimum")

plt.title("PAP: %d, Characteristic "%(PAP) + '{:.1E}'.format(k))
plt.grid()
plt.savefig(SaveDirName + "/ComparingMinima.png")



endtime = time.time()

print("Time taken:",endtime-starttime)
