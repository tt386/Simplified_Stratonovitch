import numpy as np

import shutil, os

Year = 0.1#200#0.1#200
PAP = 0.05
T = Year*(1-PAP)
Phi = 1e-5 #Initial Conditions
k = 1e-5

OSR_Proportion = 0.8

dx =1e-4#1 #1e-4#1
dt =1e-4#1 #1e-4#1


SystemSizeList = np.arange(4e-3,3e-2,1e-4)

SumTerms1 = 100
SumTerms2 = 100

#SaveFile
SaveDirName = ("Saved_Plots/" +
    "PAP_%0.3f_OSR_Proportion_%0.3f_Year_%3f_Chara_"%(PAP,OSR_Proportion,Year)
     + '{:.1E}'.format(k) + "_MaxSize_%0.5f__dx_"%
    (SystemSizeList[-1]) + '{:.1E}'.format(dx)+ "_dt_"+
    '{:.1E}'.format(dt) + "_InitialRRatio_" +'{:.1E}'.format(Phi) + 
    "_SumTerms1_%d_SumTerms2_%d"%(SumTerms1,SumTerms2))

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)

