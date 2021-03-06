import numpy as np

import shutil, os
################################
################################
################################
MinOSRNum = 2
MaxOSRNum = 20
dOSRNum = 1
OSRNum = np.arange(MinOSRNum,MaxOSRNum,dOSRNum)

OSRSeparation = 4


Min_x = -10
Max_x = 300
xNum = 1000#400#200
xlist = np.linspace(Min_x,Max_x,xNum)


LogMaxeta = -2#1e2
LogMineta = 2
etaNum = 100#1000
etalist = np.logspace(LogMineta,LogMaxeta,etaNum)

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
    "N_OSRSep_%0.3f_"
    "MinNum_%d_MaxNum_%d_dNum_%d"
    "MinX_%d_MaxX_%d_xNum_%d_"
    %(OSRSeparation,MinOSRNum,MaxOSRNum,dOSRNum,Min_x,Max_x,xNum) + 
    "_LogMaxeta_"+ '{:.1E}'.format(LogMaxeta) + "_LogMineta_" + '{:.1E}'.format(LogMineta)
    + "_etaNum_" + '{:.1E}'.format(etaNum) +
    "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi))


