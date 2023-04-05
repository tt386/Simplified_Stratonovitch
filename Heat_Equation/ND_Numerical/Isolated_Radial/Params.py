import numpy as np
import os
import shutil

Min_N = 1
Max_N = 31
NList = np.arange(Min_N,Max_N)

print(NList)

dr = 0.01

LogMaxEta = 2
LogMinEta = -4
EtaNum = 200
etalist = np.logspace(LogMinEta,LogMaxEta,EtaNum)#100

Phi = 1e-5
PAP = 0.1


SaveDirName = ("Saved_Plots/"+
                "MinN_%d_MaxN_%d"%(Min_N,Max_N) + 
                "_dr_"+'{:.1E}'.format(dr) +
                "_LogMaxEta_"+'{:.1E}'.format(LogMaxEta)+
                "_LogMinEta_"+'{:.1E}'.format(LogMinEta)+
                "_EtaNum_%d"%(EtaNum) + 
                "_PAP_" + '{:.1E}'.format(PAP) +
                "_Phi_"+'{:.1E}'.format(Phi))


if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copy("Params.py",SaveDirName)
