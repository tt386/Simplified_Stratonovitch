import numpy as np
import os,sys,shutil


T = 100
TimeList = np.linspace(0,1,100)

LList = [10,20,30,40,50,60,70,80]
P = 0.1

#Max Jump
M = 3



#Carrying Capacity
K = 100


#Mutant proportion
phi = 10/K

#Death Probabilities
dS = 0.9
dR = 0.1

Repeats = 5



SaveDirName = "SaveFiles/MinL_%d_MaxL_%d_MaxJump_%d_CarryingCapacity_%d_phi_%0.5f_P_%0.3f_dS_%0.3f_dR_%0.3f_Repeats_%d"%(LList[0],LList[-1],M,K,phi,P,dS,dR,Repeats)
