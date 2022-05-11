import numpy as np
import os
import shutil

h = 0.1

Dimension = 100 #Number of points

Maxrho = 1.0
Minrho = 0.9#0.05
drho = 0.005
rhoList = np.arange(Minrho,Maxrho,drho)

print(rhoList)

OSR_Shape = "Circle"
if OSR_Shape not in ["Circle","Rect","Square"]:
    raise ValueError('Specify a Correct OSR geometry')

#OSRSize = int(Dimension*rho)
#OSRHalfWidth = int(OSRSize/2)

"""
MaxEta = 1e-1
MinEta = 1e-5
EtaNum = 10000
etalist = np.linspace(MinEta,MaxEta,EtaNum)
"""
logMaxEta = 4
logMinEta = -2
EtaNum = 1000
etalist = np.logspace(logMinEta,logMaxEta,EtaNum)

PAP = 0.1

Phi = 1e-5


SaveDirName = ("Saved_Plots/"+
                "OSR_Shape_"+OSR_Shape+
                "_Maxrho_"+'{:.1E}'.format(Maxrho) +
                "_Minrho_"+'{:.1E}'.format(Minrho) +
                "_drho_"+'{:.1E}'.format(drho) + 
                "_h_"+'{:.1E}'.format(h)+"_Dimension_%d"%(Dimension)+
                "_logMaxEta_"+'{:.1E}'.format(logMaxEta)+
                "_logMinEta_"+'{:.1E}'.format(logMinEta)+
                "_EtaNum_%d"%(EtaNum) + "_PAP_" + '{:.1E}'.format(PAP) +
                "_Phi_"+'{:.1E}'.format(Phi))

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copy("Params.py",SaveDirName)
