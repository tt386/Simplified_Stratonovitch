import numpy as np

import shutil, os

Years = 10

PAP = 0.1
T = 1. - PAP

Phi = 1e-5 #Initial Conditions

dx =1e-4#1e-3#1 #1e-4#1
dt =1e-4#1 #1e-4#1

N = 100#50 #Number of terms in Fourier Series


logMinEta = -3
logMaxEta = 1
EtaNum = 100
etaList = np.logspace(logMinEta,logMaxEta,EtaNum)#np.arange(1e-1,1e2,1e-1)#1e-4,1e-0,1e-4)#np.arange(1e-4,5e-2,1e-4)

Maxrho = 1
Minrho = 0.05
drho = 0.05
rhoList = np.arange(Minrho,Maxrho,drho)#np.arange(0,1.0,0.1)#np.arange(0,1.0,0.01)

#SaveFile
SaveDirName = ("Saved_Plots/"+
                "Maxrho_"+'{:.1E}'.format(Maxrho) +
                "_Minrho_"+'{:.1E}'.format(Minrho) +
                "_drho_"+'{:.1E}'.format(drho) +
                "_dx_"+'{:.1E}'.format(dx) +
                "_Years_%d_"%(Years) +
                "_logMaxEta_"+'{:.1E}'.format(logMaxEta)+
                "_logMinEta_"+'{:.1E}'.format(logMinEta)+
                "_EtaNum_%d"%(EtaNum) + "_PAP_" + '{:.1E}'.format(PAP) +
                "_Phi_"+'{:.1E}'.format(Phi) + 
                "_N_%d"%(N))

