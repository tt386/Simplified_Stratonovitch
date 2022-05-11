import numpy as np

import shutil, os
################################
################################
################################
OSRNum = 1
OSRSeparation = np.asarray([0])

Max_x = 2
xNum =  100#200
xlist = np.linspace(-Max_x,Max_x,xNum)


LogMaxeta = 0#1e2
LogMineta = -5
etaNum = 1000#1000
etalist = np.logspace(LogMineta,LogMaxeta,etaNum)

PAP = 0.1
Phi = 1e-5



ParamDict = {
    "OSRNum" : OSRNum,
    "OSRSeparation" : OSRSeparation,
    "xlist" : xlist,
    "etalist" : etalist,
    "PAP" : PAP,
    "Phi" : Phi,
    "xNum" : xNum
    }



SaveDirName = ("Saved_Plots/"
    "VariableBounds2_OSRNum_%d"
    "_MaxX_%d_xNum_%d"
    %(OSRNum,Max_x,xNum) +
    "_LogMaxeta_"+ '{:.1E}'.format(LogMaxeta) + "_LogMineta_" + '{:.1E}'.format(LogMineta)
    + "_etaNum_" + '{:.1E}'.format(etaNum) +
    "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi))
