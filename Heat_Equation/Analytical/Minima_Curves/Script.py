import numpy as np
import matplotlib.pyplot as plt

from scipy import integrate


import os
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
    return (1/2*L)*integrate.simps(init, xlist)

def An(n,init,xlist,L):
    for i in range(len(xlist)):
        init[i] *= np.cos(n*np.pi*x/L)
        
    return (1/L)*integrate.simps(init, xlist)

def Bn(n,init,xlist,L):
    for i in range(len(xlist)):
        init[i] *= np.sin(n*np.pi*x/L)
        
    return (1/L)*integrate.simps(init, xlist)

def u2(x,t,L,k,init,xlist):
    tot = A0(L,init,xlist)
    for n in range(1,100):
        tot += (np.exp(-k*t*(np.pi*n/L)**2) * 
                (An(n,init,xlist,L) * np.cos(n*np.pi*x/L) + 
                Bn(n,init,xlist,L) * np.sin(n*np.pi*x/L)))
    
    return tot
##############################################################################
##############################################################################
##############################################################################
k = 1

SystemSizeList = np.arange(1,30)
OSR_Proportion = 0.25
Refuge_Proportion = 1-OSR_Proportion

PAP = 5
Year = 200
T = Year - PAP

Phi = 1e-5


#SaveFile
SaveDirName = ("Saved_Plots/" + 
    "OSR_Ratio_%0.3f_PAP_%d_Year_%d_Characteristic_%0.3f_InitialRRatio_"%
    (OSR_Proportion,PAP,Year,k) + 
    '{:.1E}'.format(Phi))

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

##############################################################################
##############################################################################
##############################################################################

R_Ratio_List = []
for SystemSize in SystemSizeList:
    OSRWidth = (1-Refuge_Proportion)*SystemSize
    RefugeWidth = Refuge_Proportion*SystemSize
    
    
    
    
    
    #Step 1
    xlist = np.arange(SystemSize)
    ylist = []
    for x in range(len(xlist)):
        ylist.append(u1(xlist[x],PAP,RefugeWidth,k,Phi))
    
    
    #Step 2
    L = SystemSize/2
        
    yprimelist = []
    for x in range(len(xlist)):
        yprimelist.append(u2(x,T,L,k,ylist,xlist))
    ##########################################################################
    ##########################################################################
    ##########################################################################
        
    
    
        
        
    
    
    RRatioList = []
    SRatioList = []
    for i in range(len(ylist)):
        RRatioList.append(Phi/(yprimelist[i]+Phi))
        SRatioList.append(1-RRatioList[-1])
    
    
    TotalRRatio = sum(RRatioList)/SystemSize    
    print(SystemSize,TotalRRatio)
    
    R_Ratio_List.append(TotalRRatio)
    
    logslope = np.log(TotalRRatio)-np.log(Phi)
    print(logslope)
    
    
    ##########################################################################
    ##########################################################################
    ##########################################################################
        
        
plt.figure()
plt.plot(SystemSizeList,R_Ratio_List)

plt.xlabel("System Size")
plt.ylabel("Second Year R Ratio")

plt.title("OSR Proportion %0.3f"%(OSR_Proportion))
plt.grid()
plt.savefig(SaveDirName + "/Minima_Curve.png")
plt.close()
               
        
