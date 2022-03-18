import numpy as np

import shutil, os
################################
################################
################################
OSRNum = 2

MinOSRSep = 0
MaxOSRSep = 0
dOSRSep = 0.2
OSRSeparation = np.arange(MinOSRSep,MaxOSRSep+dOSRSep,dOSRSep)


Max_x = 20
xNum = 100
xlist = np.linspace(-Max_x,Max_x,xNum)


Maxeta = 5e1
Mineta = 1e-1
etaNum = 100
etalist = np.linspace(Mineta,Maxeta,etaNum)

PAP = 0.1
Phi = 1e-5



ParamDict = {
    "OSRNum" = OSRNum,
    "OSRSeparation" = OSRSeparation,
    "xlist" = xlist,
    "etalist" = etalist,
    "PAP" = PAP,
    "Phi" = Phi
    }




SaveDirName = ("Saved_Plots/
    OSRNum_%d_
    MinSep_%d_MaxSep_%d_dSep_%0.3f_
    MaxX_%d_xNum_%d_"
    %(OSRNum,MinOSRSep,MaxOSRSep,dOSRSep,Max_x,xNum) + "_
    Maxeta_"+ '{0:.1f}'.format(Maxeta) + "_Mineta_" + '{0:.1f}'.format(Mineta)
    + "_etaNum_" + '{0:.1f}'.format(etaNum) + "_
    PAP_" + '{0:.1f}'.format(PAP) + "_Phi_" + '{0:.1f}'.format(Phi))


