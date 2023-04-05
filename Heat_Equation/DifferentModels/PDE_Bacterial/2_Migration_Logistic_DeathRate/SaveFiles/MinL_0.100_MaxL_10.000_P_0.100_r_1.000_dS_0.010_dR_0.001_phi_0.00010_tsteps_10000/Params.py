import numpy as np
import os,sys,shutil


r = 1
phi = 1e-4
P = 0.1

LList = np.logspace(-1,1,30)#np.arange(0.1,100,10)

tsteps = 10000
dt = 1/tsteps
TList = np.linspace(0,1,tsteps)

dS = 0.01
dR = 0.001



dx = 2 * np.sqrt(2*dt)
MaxX = 1000

print("dx,",dx)
print("dt,",dt)

XList = np.arange(-MaxX,MaxX,dx)


SaveDirName = "SaveFiles/MinL_%0.3f_MaxL_%0.3f_P_%0.3f_r_%0.3f_dS_%0.3f_dR_%0.3f_phi_%0.5f_tsteps_%d"%(min(LList),max(LList),P,r,dS,dR,phi,tsteps)

if not os.path.isdir(SaveDirName):
    os.mkdir(SaveDirName)
    print("Created Directory")

shutil.copy("Params.py",SaveDirName)
