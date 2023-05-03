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

NList = []
dRList = []

for i in dirlist:
    try:
        with np.load(os.path.join(i,filename)) as data:
            N = data["N"]
            print("Number of regions being plotted:",N)
            PAP = data["PAP"]
            Phi = data["Phi"]
            C = data["C"]
            d = data["d"]
            dx = data["dx"]
            dt = data["dt"]
            PAPSNum = data["PAPSNum"]
            ENDSNum = data["ENDSNum"]
            xlist = data["xlist"]
            dR = data["dR"]
            timetaken = data["timetaken"]


    #Collect data
    NList.append(N)
    dRList.append(dR)

    #Plot the EndStatesetc

    plt.figure()
    plt.plot(xlist*dx,PAPSNum)
    plt.xlabel("Space")
    plt.ylabel("Number")
    plt.title("PAP Dist")
    plt.savefig(str(i) + "/PAP.png")
    plt.close()

    plt.figure()
    plt.plot(xlist*dx,ENDSNum)
    plt.xlabel("Space")
    plt.ylabel("Number")
    plt.title("End Dist")
    plt.savefig(str(i) + "/End.png")
    plt.close()

    plt.figure()
    plt.plot(xlist*dx,Phi/(Phi+ENDSNum))
    plt.xlabel("Space")
    plt.ylabel("Number")
    plt.title("Changing R")
    plt.savefig(str(i) + "/DR.png")
    plt.close()


NList,dRList = zip(*sorted(zip(CList, MeanList)))


plt.figure()
plt.semilogy(NList,dRList)
plt.xlabel("Number of Patches")
plt.ylabel("Change in R")
plt.savefig(str(args.directory) + /"dR.png")
plt.close()

