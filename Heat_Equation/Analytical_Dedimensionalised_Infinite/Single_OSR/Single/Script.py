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
parser.add_argument('-S','--Sep',type=float,required=True,help='Sep')
args = parser.parse_args()

OSRSeparation = args.Sep


SepSaveDirName = (SaveDirName + 
    "/OSRSeparation_" + '{:.1E}'.format(OSRSeparation))

if not os.path.isdir(SepSaveDirName):
    os.mkdir(SepSaveDirName)
    print("Created Directory for Sep",OSRSeparation)


ParamDict["OSRSeparation"] = OSRSeparation


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

OutputDatafilename = SepSaveDirName + '/datafile.npz'
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

