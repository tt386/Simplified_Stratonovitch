import os
import time

import copy
import sys
sys.path.insert(0,'../CoreFunctions')

from Core import Core
from Params import *


starttime = time.time()

#################################
###Argparse
#################################
from argparse import ArgumentParser
parser = ArgumentParser(description='Different OSR Separation')
parser.add_argument('-P','--Phi',type=float,required=True,help='Phi')
args = parser.parse_args()

Phi = args.Phi

print("PAP:",PAP)

PhiSaveDirName = (SaveDirName + 
    "/Phi_" + '{:.1E}'.format(Phi))

if not os.path.isdir(PhiSaveDirName):
    os.mkdir(PhiSaveDirName)
    print("Created Directory for Phi",Phi)


ParamDict["Phi"] = Phi


#################################
###Main Process##################
#################################
ResultsDict = Core(ParamDict)

endtime = time.time()

timetaken = endtime-starttime
#################################
###Unpack and Save###############
#################################
PAPMatrix = ResultsDict["PAPMatrix"]
YearMatrix = ResultsDict["YearMatrix"]
RDifferenceList = ResultsDict["RDifference"]

OutputDatafilename = PhiSaveDirName + '/datafile.npz'
np.savez(OutputDatafilename,
    PAP=PAP,
    Phi=Phi,
    OSRSeparation=OSRSeparation,
    OSRNum=OSRNum,
    etalist=etalist,
    PAPMatrix=PAPMatrix,
    YearMatrix=YearMatrix,
    RDifferenceList=RDifferenceList,
    xlist=xlist,
    timetaken=timetaken)

