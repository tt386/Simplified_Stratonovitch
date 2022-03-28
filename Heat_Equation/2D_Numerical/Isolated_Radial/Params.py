import numpy as np
import os
import shutil


dr = 0.005

LogMaxEta = 0
LogMinEta = -2
EtaNum = 20
etalist = np.logspace(LogMinEta,LogMaxEta,EtaNum)#100

Phi = 1e-5
PAP = 0.1


SaveDirName = ("Saved_Plots/"+
                "dr_"+'{:.1E}'.format(dr) +
                "_LogMaxEta_"+'{:.1E}'.format(LogMaxEta)+
                "_LogMinEta_"+'{:.1E}'.format(LogMinEta)+
                "_EtaNum_%d"%(EtaNum) + 
                "_PAP_" + '{:.1E}'.format(PAP) +
                "_Phi_"+'{:.1E}'.format(Phi))


if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copy("Params.py",SaveDirName)
