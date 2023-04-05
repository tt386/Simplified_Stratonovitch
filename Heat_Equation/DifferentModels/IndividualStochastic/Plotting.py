import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

plt.rcParams['text.usetex'] = True

import numpy as np
import time

from scipy.stats import linregress

from scipy.optimize import curve_fit

import sys

import scipy.integrate as integrate
import scipy.special as special

from scipy.signal import argrelextrema

import os

starttime = time.time()
################################
##ArgParse######################
################################
import os.path

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The Directory %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

from argparse import ArgumentParser

parser = ArgumentParser(description='Plotting')
parser.add_argument('-d','--directory',help='The directory of the data')
args = parser.parse_args()

###############################
##Extract Data#################
###############################
filename = 'datafile.npz'


#Find list of all the datafiles
tempdirlist = os.listdir(args.directory)
dirlist = []
for i in tempdirlist:
    if os.path.isdir(os.path.join(args.directory,i)):
        dirlist.append(os.path.join(args.directory,i))

print("Dirlist:",dirlist)


LList = []
MeanDeltaRList = []
MedianDeltaRList = []
for i in dirlist:
    try:
        with np.load(os.path.join(i,filename)) as data:
            phi = data["phi"]
            L = data["L"]
            DeltaRList = data["DeltaRList"]
            print(DeltaRList)
            LList.append(L)

            MeanDeltaRList.append(np.mean(DeltaRList))
            MedianDeltaRList.append(np.median(DeltaRList))
        print("L value:",L)

    except:
        print("Whoops")
    

#Sort these lists
zipped_lists = zip(LList, MeanDeltaRList,MedianDeltaRList)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
LList, MeanDeltaRList, MedianDeltaRList= [ list(tuple) for tuple in  tuples]

LList = np.asarray(LList)
MeanDeltaRList = np.asarray(MeanDeltaRList) 
MedianDeltaRList = np.asarray(MedianDeltaRList)


plt.figure() 
plt.plot(LList,MeanDeltaRList/LList,label='Mean')
plt.plot(LList,MedianDeltaRList/LList,label='Median')
plt.legend(loc='lower right')
plt.savefig(str(args.directory)+"/MinPlot.png")
plt.close()

"""
min_index = np.argmin(abs(DeltaRPropList/LList))

print("Min L value is",LList[min_index])

plt.figure()
plt.loglog(LList,abs(DeltaRPropList/LList))
plt.savefig(str(args.directory)+"/MinPropPlot.png")
plt.close()


plt.figure()
plt.loglog(LList,abs(DeltaRList/DeltaSList)/LList)
plt.savefig(str(args.directory)+"/ProportionOfChangesPlot.png")
plt.close()

plt.figure()
plt.loglog(LList,abs(PopDiffList/LList))
plt.savefig(str(args.directory)+"/PopDiff.png")
plt.close()
"""
