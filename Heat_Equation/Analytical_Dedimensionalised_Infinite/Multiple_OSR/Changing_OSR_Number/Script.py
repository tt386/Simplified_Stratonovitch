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
parser.add_argument('-N','--Num',type=float,required=True,help='Sep')
args = parser.parse_args()

OSRNum = args.Num


NumSaveDirName = (SaveDirName + 
    "/OSRNum_" + '{:.1E}'.format(OSRNum))

if not os.path.isdir(NumSaveDirName):
    os.mkdir(NumSaveDirName)
    print("Created Directory for Sep",OSRSeparation)


ParamDict["OSRNum"] = OSRNum


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

OutputDatafilename = NumSaveDirName + '/datafile.npz'
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

