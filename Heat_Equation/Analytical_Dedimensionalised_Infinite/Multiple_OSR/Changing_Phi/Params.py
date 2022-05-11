import numpy as np

import shutil, os
################################
################################
################################
OSRNum = 4#4

logMinPhi = -6
logMaxPhi = -1
PhiNum = 30
Phi = np.logspace(logMinPhi,logMaxPhi,PhiNum)

"""
MinOSRSep = 0
MaxOSRSep = 4
dOSRSep = 0.2
OSRSeparation = np.arange(MinOSRSep,MaxOSRSep+dOSRSep,dOSRSep)
"""

SingleOSR = True
if SingleOSR:
    OSRNum = 1
    OSRSeparation = np.asarray([0])

Max_x = 2#10
xNum = 100#200
xlist = np.linspace(-Max_x,Max_x,xNum)


LogMaxeta = 0#1e2
LogMineta = -4
etaNum = 1000#1000
etalist = np.logspace(LogMineta,LogMaxeta,etaNum)

#PAP = 0.1
PAP = 0.1



ParamDict = {
    "OSRNum" : OSRNum,
    "OSRSeparation" : OSRSeparation,
    "xlist" : xlist,
    "etalist" : etalist,
    "PAP" : PAP,
    "Phi" : Phi,
    "xNum" : xNum
    }



if not SingleOSR:
    SaveDirName = ("Saved_Plots/"
        "NEWXBOUNDS_OSRNum_%d_"
        "MinSep_%d_MaxSep_%d_dSep_%0.3f_"
        "MaxX_%d_xNum_%d_"
        %(OSRNum,MinOSRSep,MaxOSRSep,dOSRSep,Max_x,xNum) + 
        "_LogMaxeta_"+ '{:.1E}'.format(LogMaxeta) + "_LogMineta_" + '{:.1E}'.format(LogMineta)
        + "_etaNum_" + '{:.1E}'.format(etaNum) +
        "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi))

else:
    SaveDirName = ("Saved_Plots/"
        "OSRNum_%d"
        "_MaxX_%d_xNum_%d"
        %(OSRNum,Max_x,xNum) +
        "_LogMaxeta_"+ '{:.1E}'.format(LogMaxeta) + "_LogMineta_" + '{:.1E}'.format(LogMineta)
        + "_etaNum_" + '{:.1E}'.format(etaNum) +
        "_LogMaxPhi_" + '{:.1E}'.format(logMaxPhi) + "_LogMinPhi_" + '{:.1E}'.format(logMinPhi) +
        "_PhiNum_" + '{:.1E}'.format(PhiNum) + 
        "_PAP_" + '{:.1E}'.format(PAP))
