import numpy as np

import shutil, os

Year = 0.1#200#0.1#200
PAP = (Year/2)
T = Year - PAP

Phi = 1e-5 #Initial Conditions


OSR_Proportion = 0.5

dx =1e-4#1 #1e-4#1
dt =1e-4#1 #1e-4#1



SystemSizeList = np.arange(1e-3,2e-2,1e-4)#np.arange(1,100,2)#np.arange(0.01,0.5,0.1)#np.arange(1,100,2)
CharacteristicList = np.arange(1e-6,1e-4,1e-6)

#OSR_ProportionList = np.arange(0.1,0.9,0.05)




#SaveFile
SaveDirName = ("Saved_Plots/" +
    "PAP_%0.3f_Year_%3f_OSR_Proportion_%0.3f"%(PAP,Year,OSR_Proportion) +
    "_MaxSize_%0.5f_MaxChar_%0.5f_dx_"%(SystemSizeList[-1],CharacteristicList[-1]) +
    '{:.1E}'.format(dx)+ "_dt_"+'{:.1E}'.format(dt) + "_InitialRRatio_" +
    '{:.1E}'.format(Phi))

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copy("Params.py",SaveDirName)
