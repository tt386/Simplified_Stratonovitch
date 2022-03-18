
import numpy as np

import matplotlib.pyplot as plt

from scipy import integrate

from copy import copy

import os
##############################################################################
##############################################################################
##############################################################################
#Variables

YearLength = 1
SystemLengths = [0.1]#np.arange(0.1,5,0.1)

OSR_Ratio = 0.1

Characteristic = 2e-4
dt = 1e-3
dx = 1e-3

MigrationsPerDay = 50


R_Ratio = 1e-5

#Derived Variables
days = int(YearLength/dt)
F = Characteristic * (dt)/(dx)**2


SlopeList = []
for SystemLength in SystemLengths:
    
    ##############################################################################
    ##############################################################################
    ##############################################################################
    #SaveFileCreation
    
    SaveFileName = ("Saved_Plots/OSRProportion_" + '{:.1E}'.format(OSR_Ratio)
                    + "_SystemSize_%0.3f_InitialRRatio_"%(SystemLength)
                    + '{:.1E}'.format(R_Ratio) + 
                    "_Characteristic_%0.5f_MigrationsPerDay_%d"%(Characteristic,MigrationsPerDay))
    if not os.path.isdir(SaveFileName):
        os.mkdir(SaveFileName)
    
    
    print(SystemLength)
    systemsites = int(SystemLength/dx)
    
    OSRWidth = int(systemsites * OSR_Ratio)
    RefugeSize = int(systemsites * (1-OSR_Ratio))

    systemsites = OSRWidth + RefugeSize
    
    
    
    
    S_Domain = np.ones(systemsites)*(1-R_Ratio)
    R_Domain = np.ones(systemsites)*(R_Ratio)
    
    
    PesticideDomain = np.ones(OSRWidth)
    PesticideDomain = np.pad(PesticideDomain,(0,RefugeSize))
    
    #print(S_Domain)
    
    
    ##############################################################################
    ##############################################################################
    ##############################################################################
    fignum = 0
    for t in range(days):
        #print(t)
        for tm in range(MigrationsPerDay):
            plt.figure()
            plt.plot(np.arange(len(S_Domain)),S_Domain)
            plt.grid()
            plt.savefig(SaveFileName + "/Fig_"+str(fignum).zfill(6))
            plt.close()
            fignum +=1 
            
            #Migration
            S_Domain = (S_Domain + 
                        F*(np.roll(S_Domain,-1) - 2 * S_Domain + 
                           np.roll(S_Domain,1)))
    
        #Death From Pesticides
        S_Domain[PesticideDomain == 1] = 0
    
    TotalR = np.sum(R_Domain/(R_Domain+S_Domain))
    
    SlopeList.append(TotalR/np.sum(R_Ratio*systemsites))
    
    print("OSR Ratio,",OSR_Ratio," Total  R",TotalR)
    print("Corresponding slope:",TotalR/np.sum(R_Ratio*systemsites))
