
import numpy as np

import matplotlib.pyplot as plt

from scipy import integrate

from copy import copy

import os

import time

starttime = time.time()
##############################################################################
##############################################################################
##############################################################################
#Variables

YearLength = 0.1
SystemLengths = np.arange(0.005,0.2,0.002)

OSR_RatioList = np.arange(0.1,0.9,0.05)

Characteristic = 2e-5
dt = 1e-4
dx = 1e-4

MigrationsPerDay = 1


R_Ratio = 1e-5

#Derived Variables
days = int(YearLength/dt)
PAP = int(days/2)

F = Characteristic * (dt)/(dx)**2


#SaveFile
SaveFileName = ("Saved_Plots/InitialRRatio_"
                + '{:.1E}'.format(R_Ratio) + "_MigrationsPerDay_%d"%
                (MigrationsPerDay) + "_dx_" + '{:.1E}'.format(dx)+
                "_dt_"+'{:.1E}'.format(dt))

if not os.path.isdir(SaveFileName):
    os.mkdir(SaveFileName)




MinList = []
for OSR_Ratio in OSR_RatioList:
    SlopeList = []
    for SystemLength in SystemLengths:
        
        ##############################################################################
        ##############################################################################
        ##############################################################################
        
        print(SystemLength)
        systemsites = int(SystemLength/dx)
        
        OSRWidth = int(systemsites * OSR_Ratio)
        RefugeSize = int(systemsites * (1-OSR_Ratio))

        systemsites = OSRWidth + RefugeSize
        
        
        
        
        S_Domain = np.ones(systemsites)*(1-R_Ratio)
        R_Domain = np.ones(systemsites)*(R_Ratio)
        
        
        PesticideDomain = np.ones(OSRWidth)
        PesticideDomain = np.pad(PesticideDomain,(0,RefugeSize),mode='constant')
        
        #print(S_Domain)
        
        
        ##############################################################################
        ##############################################################################
        ##############################################################################
        for t in range(days):
            #print(t)
            if t < PAP:
                S_Domain[PesticideDomain == 1] = 0
            for tm in range(MigrationsPerDay):
                
                S_Domain = (S_Domain + 
                            F*(np.roll(S_Domain,-1) - 2 * S_Domain + 
                               np.roll(S_Domain,1)))
        
            
        
        R_Breeding = R_Domain/(R_Domain+S_Domain)
        S_Breeding = S_Domain/(R_Domain+S_Domain)
        
        
        TotalR = np.sum(R_Domain/(R_Domain+S_Domain))
        
        
        
        SlopeList.append(np.log(TotalR) - np.log(R_Ratio*systemsites))
        
        print("OSR Ratio,",OSR_Ratio," Total  R",TotalR)
        print("Initial R",R_Ratio*systemsites)
        print("Corresponding slope:",TotalR/np.sum(R_Ratio*systemsites))


    plt.figure()
    plt.plot(SystemLengths,SlopeList)

    plt.xlabel("System Length")
    plt.ylabel("R Proportion")
    plt.title("OSR Proportion: %0.3f"%(OSR_Ratio))
    plt.grid()
    plt.savefig(SaveFileName + "/OSR_Ratio%0.3f.png"%(OSR_Ratio))
    plt.close()



    minval = 1000000
    minsize = 0
    for i in range(len(SlopeList)):
        if SlopeList[i] < minval:
            minval = SlopeList[i]
            minsize = SystemLengths[i]
    MinList.append(minsize)    




plt.figure()
plt.plot(OSR_RatioList,MinList)

plt.xlabel("OSR Ratio")
plt.ylabel("System Lenth at Min Value")
#plt.title("OSR Proportion: %0.3f"%(OSR_Ratio))
plt.grid()
plt.savefig(SaveFileName + "/MinimaOfMinima.png")
plt.close()


endtime = time.time()

print("Time taken:",endtime-starttime)
