import numpy as np

import shutil, os
################################
################################
################################
OSRSeparation = 1#4

OSRNum = 1


SingleOSR = False
if SingleOSR:
    OSRNum = 1
    OSRSeparation = np.asarray([0])

x_buffer = 10
dx = 1e-3

x_buffer2 = 10#30
x_buffer2Factor = 1.25#2


"""
LogMaxeta = 6#1e2
LogMineta = -2
etaNum = 1000
etalist = np.logspace(LogMineta,LogMaxeta,etaNum)
"""
etalist = [1/(20**2)]

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


SaveDirName = ("Saved_Plots/"
    "DynamicIntegration_ExponentialBuffer_DYNAMICXBOUNDS_LOG_"
    "XBuffer_%d_dx_%0.3f_XBuffer2_%d_XBuffer2Factor_%0.3f"
    %(x_buffer,dx,x_buffer2,x_buffer2Factor) +
    "_eta_"+ '{:.1E}'.format(etalist[0]) +
    "_PAP_" + '{:.1E}'.format(PAP) + "_Phi_" + '{:.1E}'.format(Phi) +
    "_N_%d"%(N))




if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copy("Params.py",SaveDirName)
