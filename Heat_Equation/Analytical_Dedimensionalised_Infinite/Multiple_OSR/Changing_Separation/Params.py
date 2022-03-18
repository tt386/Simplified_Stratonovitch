import numpy as np

import shutil, os
################################
################################
################################
OSRNum = 4

MinOSRSep = 0
MaxOSRSep = 4
dOSRSep = 0.2
OSRSeparation = np.arange(MinOSRSep,MaxOSRSep+dOSRSep,dOSRSep)


Max_x = 50
xNum = 400#200
xlist = np.linspace(-Max_x,Max_x,xNum)


Maxeta = 5e-1#1e2
Mineta = 1e-1
etaNum = 100#1000
etalist = np.linspace(Mineta,Maxeta,etaNum)

PAP = 0.1
Phi = 1e-5



ParamDict = {
    "OSRNum" : OSRNum,
    "OSRSeparation" : OSRSeparation,
    "xlist" : xlist,
    "etalist" : etalist,
    "PAP" : PAP,
    "Phi" : Phi
    }




SaveDirName = ("Saved_Plots/"
    "OSRNum_%d_"
    "MinSep_%d_MaxSep_%d_dSep_%0.3f_"
    "MaxX_%d_xNum_%d_"
    %(OSRNum,MinOSRSep,MaxOSRSep,dOSRSep,Max_x,xNum) + 
    "_Maxeta_"+ '{:.1E}'.format(Maxeta) + "_Mineta_" + '{:.1E}'.format(Mineta)
    + "_etaNum_" + '{:.1E}'.format(etaNum) +
    "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi))


