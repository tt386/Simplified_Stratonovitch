import numpy as np

import shutil, os
################################
################################
################################
OSRSeparation = 1#4

MinOSRNum = 1
MaxOSRNum = 42
dNum = 2
OSRNum=np.arange(MinOSRNum,MaxOSRNum,dNum)


SingleOSR = True
if SingleOSR:
    OSRNum = np.asarray([1])
    OSRSeparation = 0

x_buffer = 100
dx = 1e-1

x_buffer2 = 100#30
x_buffer2Factor = 1.25#2

"""
x_buffer = 200
dx = 1e-1
"""

LogMaxeta = np.log10(1/((1.1)**2))
LogMineta = np.log10(1/((1.2)**2))
etaNum = 1000#1000
etalist = np.logspace(LogMineta,LogMaxeta,etaNum)

PAP = 0.1
Phi = 1e-5

N = 50#100

ParamDict = {
    "OSRNum" : OSRNum,
    "OSRSeparation" : OSRSeparation,
    #"xlist" : xlist,
    "etalist" : etalist,
    "PAP" : PAP,
    "Phi" : Phi,
    #"xNum" : xNum,
    "x_buffer" : x_buffer,
    "x_buffer2": x_buffer2,
    "dx" : dx,
    "N" : N
    }


"""
if not SingleOSR:
    SaveDirName = ("Saved_Plots/"
        "NEWXBOUNDS_OSRNum_%d_"
        "MinSep_%d_MaxSep_%d_dSep_%0.3f_"
        "MaxX_%d_xNum_%d_"
        %(OSRNum,MinOSRSep,MaxOSRSep,dOSRSep,Max_x,xNum) + 
        "_LogMaxeta_"+ '{:.1E}'.format(LogMaxeta) + "_LogMineta_" + '{:.1E}'.format(LogMineta)
        + "_etaNum_" + '{:.1E}'.format(etaNum) +
        "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi) +
        "_N_%d"%(N))

else:
    SaveDirName = ("Saved_Plots/"
        "NEWXBOUNDS_OSRNum_%d"
        "_MaxX_%d_xNum_%d"
        %(OSRNum,Max_x,xNum) +
        "_LogMaxeta_"+ '{:.1E}'.format(LogMaxeta) + "_LogMineta_" + '{:.1E}'.format(LogMineta)
        + "_etaNum_" + '{:.1E}'.format(etaNum) +
        "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi) + 
        "_N_%d"%(N))
"""
"""
SaveDirName = ("Saved_Plots/"
    "LOG_OSRNum_%d_"
    "LogMinSep_%d_LogMaxSep_%d_SepNum_%d_"
    "MaxX_%d_xNum_%d_"
    %(OSRNum,LogMinOSRSep,LogMaxOSRSep,SepNum,Max_x,xNum) +
    "_LogMaxeta_"+ '{:.1E}'.format(LogMaxeta) + "_LogMineta_" + '{:.1E}'.format(LogMineta)
    + "_etaNum_" + '{:.1E}'.format(etaNum) +
    "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi) +
    "_N_%d"%(N))
"""

SaveDirName = ("Saved_Plots/"
    "DynamicIntegration_ExponentialBuffer_DYNAMICXBOUNDS_LOG_"
    "XBuffer_%d_dx_%0.3f_XBuffer2_%d_XBuffer2Factor_%0.3f"
    %(x_buffer,dx,x_buffer2,x_buffer2Factor) +
    "_LogMaxeta_"+ '{:.1E}'.format(LogMaxeta) + "_LogMineta_" + '{:.1E}'.format(LogMineta)
    + "_etaNum_" + '{:.1E}'.format(etaNum) +
    "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi) +
    "_N_%d"%(N))
