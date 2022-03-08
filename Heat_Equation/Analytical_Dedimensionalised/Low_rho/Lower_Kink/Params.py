import numpy as np

import shutil, os

PAP = 0.1
T = 1. - PAP

Phi = 1e-5 #Initial Conditions

dx =1e-4#1e-3#1 #1e-4#1
dt =1e-4#1 #1e-4#1

N = 100#50 #Number of terms in Fourier Series


UpperEta = 6e-2
LowerEta = 1e-3
EtaNum = 1000
etaList = np.linspace(LowerEta,UpperEta,EtaNum)#np.arange(1e-1,1e2,1e-1)#1e-4,1e-0,1e-4)#np.arange(1e-4,5e-2,1e-4)


UpperRho = 6e-2
LowerRho = 4e-2
RhoNum = 20
rhoList = np.linspace(LowerRho,UpperRho,RhoNum)#np.arange(0,1.0,0.1)#np.arange(0,1.0,0.01)

#SaveFile
SaveDirName = ("Saved_Plots/" + 
    "Rho_" + '{:.1E}'.format(LowerRho) + "-" + '{:.1E}'.format(UpperRho) + "_RhoNum_%d_"%(RhoNum) +
    "Eta_" + '{:.1E}'.format(LowerEta) + "-" + '{:.1E}'.format(UpperEta) + "_EtaNum_%d_"%(EtaNum) +
    "PAP_%0.3f_Phi_"%(PAP) + '{:.1E}'.format(Phi) + "_N_%d"%(N))

"""
SaveDirName = ("Saved_Plots/" +
    "PAP_%0.3f_"%(PAP) + "_Maxeta_%0.5f_Maxrho_%0.5f_dx_"%(etaList[-1],
    rhoList[-1]) + '{:.1E}'.format(dx)+ "_dt_"+'{:.1E}'.format(dt) + 
    "_InitialRRatio_" +'{:.1E}'.format(Phi) + "_N_%d"%(N))
"""
