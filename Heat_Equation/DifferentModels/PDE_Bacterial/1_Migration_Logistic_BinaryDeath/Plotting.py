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

        XList = data["XList"]

        LList.append(L)
        DeltaRList.append(DeltaR)

        PopDiffList.append(np.sum(1-RList-SList))

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

#Sort these lists
zipped_lists = zip(LList, DeltaRList, PopDiffList)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
LList, DeltaRList, PopDiffList= [ list(tuple) for tuple in  tuples]

LList = np.asarray(LList)
DeltaRList = np.asarray(DeltaRList) 
PopDiffList = np.asarray(PopDiffList)


plt.figure() 
plt.loglog(LList,DeltaRList/LList)
plt.savefig(str(args.directory)+"/MinPlot.png")
plt.close()

plt.figure()
plt.loglog(LList,PopDiffList/LList)
plt.savefig(str(args.directory)+"/PopDiff.png")
plt.close()
