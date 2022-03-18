import numpy as np

import shutil, os

Year = 0.1#200#0.1#200
PAP = (Year/2)
T = Year - PAP

Phi = 1e-5 #Initial Conditions


k = 1e-5

dx =1e-4#1 #1e-4#1
dt =1e-4#1 #1e-4#1


SystemSizeList = np.arange(4e-3,6e-2,1e-4)#np.arange(1,100,2)#np.arange(0.01,0.5,0.1)#np.arange(1,100,2)

OSRRatioList = np.arange(0,1.0,0.01)

#OSR_ProportionList = np.arange(0.1,0.9,0.05)




#SaveFile
SaveDirName = ("Saved_Plots/" +
    "PAP_%0.3f_Year_%3f_Chara_"%(PAP,Year) + '{:.1E}'.format(k) +
    "_MaxSize_%0.5f_MaxOSRProp_%0.5f_dx_"%(SystemSizeList[-1],OSRRatioList[-1]) +
    '{:.1E}'.format(dx)+ "_dt_"+'{:.1E}'.format(dt) + "_InitialRRatio_" +
    '{:.1E}'.format(Phi))

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copy("Params.py",SaveDirName)
