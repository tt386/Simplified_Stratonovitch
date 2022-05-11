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
parser.add_argument('-P','--PAP',type=float,required=True,help='PAP')
args = parser.parse_args()

PAP = args.PAP

print("PAP:",PAP)

PAPSaveDirName = (SaveDirName + 
    "/PAP_" + '{:.1E}'.format(PAP))

if not os.path.isdir(PAPSaveDirName):
    os.mkdir(PAPSaveDirName)
    print("Created Directory for PAP",PAP)


ParamDict["PAP"] = PAP


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

OutputDatafilename = PAPSaveDirName + '/datafile.npz'
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

