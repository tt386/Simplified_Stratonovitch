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
DeltaRList = []
DeltaSList = []
DeltaRPropList = []
PopDiffList = []
for i in dirlist:
    #try:
    with np.load(os.path.join(i,filename)) as data:
        phi = data["phi"]
        r = data["r"]
        L = data["L"]
        DeltaR = data["DeltaR"]

        RList = data["RList"]
        SList = data["SList"]

        SPAPDist = data["SPAPDist"]
        RPAPDist = data["RPAPDist"]

        XList = data["XList"]

        LList.append(L)
        DeltaRList.append(DeltaR)

        dx = XList[1] - XList[0]

        DeltaSList.append(np.sum((1-phi) - SList)* dx)


        PopDiffList.append(np.sum(1-RList-SList))


        DeltaRProp = np.sum(RList/(SList+RList)-phi) * dx

        DeltaRPropList.append(DeltaRProp)

    print("L value:",L)

    
    print(len(XList))
    print(len(SPAPDist))
    plt.figure()
    plt.semilogy(XList,SPAPDist,linewidth=5)
    plt.semilogy(XList,RPAPDist,linewidth=5)
    plt.xlim(-L*5,L*5)
    plt.ylim(1e-5,1.1)
    plt.savefig(str(i) + "/PAPDist.png")

    plt.close()
    


    plt.figure()
    plt.semilogy(XList,SList,linewidth=5)
    plt.semilogy(XList,RList,linewidth=5)
    plt.xlim(-L*5,L*5)
    plt.ylim(1e-5,1.1)
    plt.savefig(str(i) + "/EndState.png")

    plt.close()


    plt.figure()
    plt.semilogy(XList,RList-phi,linewidth=5)
    plt.xlim(-L*5,L*5)
    plt.ylim(1e-5,1.1)
    plt.savefig(str(i) + "/ChangeInR.png")
    plt.close()

    plt.figure()
    plt.semilogy(XList,RList/(RList+SList)-phi,linewidth=5)
    plt.xlim(-L*5,L*5)
    plt.ylim(1e-5,1.1)
    plt.savefig(str(i) + "/ChangeInRProp.png")
    plt.close()

#Sort these lists
zipped_lists = zip(LList, DeltaRList,DeltaSList, PopDiffList, DeltaRPropList)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
LList, DeltaRList,DeltaSList, PopDiffList, DeltaRPropList= [ list(tuple) for tuple in  tuples]

LList = np.asarray(LList)
DeltaRList = np.asarray(DeltaRList) 
PopDiffList = np.asarray(PopDiffList)
DeltaRPropList = np.asarray(DeltaRPropList)
DeltaSList = np.asarray(DeltaSList)

plt.figure() 
plt.loglog(LList,abs(DeltaRList/LList))
plt.savefig(str(args.directory)+"/MinPlot.png")
plt.close()


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
