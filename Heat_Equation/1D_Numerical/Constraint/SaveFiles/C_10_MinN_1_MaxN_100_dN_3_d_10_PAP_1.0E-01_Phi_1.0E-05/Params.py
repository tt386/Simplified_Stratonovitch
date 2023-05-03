import numpy as np
import os
import shutil



#Total width of Selection Region
C = 10

#Number of sub Selection regions
NList = np.arange(1,101,3)

#Separations
d = 10#dList = np.arange(1,11)

#Make dx 10 times smaller than the smallest gap: C/N=L
dx = min(C/NList)/10

#Corresponding dt for VOn Neumann stability analysis
dt = dx**2 /10


#XList
xlist = np.arange(2*(1/dt + C*max(NList)*d))



#Initial proportion of R
Phi = 1e-5

#Time for which the pesticide is applied
PAP = 0.1


SaveDirName = ("SaveFiles/"+
            "C_%d"%(C) +
            "_MinN_%d_MaxN_%d_dN_%d"%(min(NList),max(NList),NList[1]-NList[0]) +
            "_d_%d"%(d) +
            "_PAP_" + '{:.1E}'.format(PAP) +
            "_Phi_"+'{:.1E}'.format(Phi))


if not os.path.isdir("SaveFiles"):
    os.mkdir("SaveFiles")


if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copyfile("Params.py", SaveDirName+'/Params.py')
shutil.copyfile("Script.py", SaveDirName+"/Script.py")

