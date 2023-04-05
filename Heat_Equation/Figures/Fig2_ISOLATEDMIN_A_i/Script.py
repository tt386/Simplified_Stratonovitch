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
parser.add_argument('-S','--Num',type=float,required=True,help='Sep')
args = parser.parse_args()

OSRNum = args.Num


SepSaveDirName = (SaveDirName + 
    "/OSRNum_" + '{:.1E}'.format(OSRNum))

if not os.path.isdir(SepSaveDirName):
    os.mkdir(SepSaveDirName)
    print("Created Directory for Num",OSRNum)


ParamDict["OSRNum"] = OSRNum

"""
#Define dynamic x-list:
xlist = np.arange(-x_buffer,OSRNum + (OSRNum-1)*OSRSeparation + x_buffer,dx)

#Add on the exponentially distributed buffer:
nlist = np.arange(x_buffer2)

ConcatList = dx*x_buffer2Factor**(nlist)

PreList = xlist[0] - np.flip(ConcatList)

PostList = xlist[-1] + ConcatList

xlist = np.concatenate((PreList,xlist,PostList))
"""
#Define dynamic x-list:
LowerLimit = -x_buffer
UpperLimit = OSRNum + (OSRNum-1)*OSRSeparation + x_buffer
xlist = np.arange((LowerLimit+UpperLimit)/2,UpperLimit,dx)


#Add on the exponentially distributed buffer:
nlist = np.arange(x_buffer2)

ConcatList = dx*x_buffer2Factor**(nlist)

PostList = xlist[-1] + ConcatList

xlist = np.concatenate((xlist,PostList))

ParamDict["xlist"] = xlist

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

