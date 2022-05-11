import numpy as np

import shutil, os
################################
################################
################################
LogMinOSRNum = 0
LogMaxOSRNum = 1.5
OSRNum_Num = 20
OSRNum = np.logspace(LogMinOSRNum,LogMaxOSRNum,OSRNum_Num)

logOSRSepMin = -1
logOSRSepMax = 2
SepNum = 10
OSRSeparationList = np.logspace(logOSRSepMin,logOSRSepMax,SepNum)


Min_x = -10
Max_x = 5000
xNum = 100000#400#200
xlist = np.linspace(Min_x,Max_x,xNum)


"""
LogMaxeta = -2#1e2
LogMineta = 2
etaNum = 100#1000
etalist = np.logspace(LogMineta,LogMaxeta,etaNum)
"""

L = 1e1 #Total area of the OSR
eta = 1/L**2

PAP = 0.1
Phi = 1e-5

N = 1000

ParamDict = {
    "OSRNum" : OSRNum,
    "OSRSeparationList" : OSRSeparationList,
    "xlist" : xlist,
    "L" : L,
    "eta" : eta,
    "PAP" : PAP,
    "Phi" : Phi,
    "N" : N
    }




SaveDirName = ("Saved_Plots/"
    "L_" + '{:.1E}'.format(L) + "_eta_" + '{:.1E}'.format(eta) +
    "_LogMaxSep_"+ '{:.1E}'.format(logOSRSepMax) + "_LogMinSep_" + '{:.1E}'.format(logOSRSepMin) + 
    "_SepNum_" + '{:.1E}'.format(SepNum) +
    "_LogMinNum_%d_LogMaxNum_%d_dNum_%d"
    "_MinX_%d_MaxX_%d_xNum_%d"
    %(LogMinOSRNum,LogMaxOSRNum,OSRNum_Num,Min_x,Max_x,xNum) + 
    "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi) +
    "_N_%d"%(N))


