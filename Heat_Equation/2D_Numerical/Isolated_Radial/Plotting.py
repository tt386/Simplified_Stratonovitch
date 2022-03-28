import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
import time

from scipy.stats import linregress

from scipy.optimize import curve_fit

import sys

import scipy.integrate as integrate
import scipy.special as special

from scipy.signal import argrelextrema

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

with np.load(os.path.join(str(args.directory),filename)) as data:
    PAP = data["PAP"]
    Phi = data["Phi"]
    dr = data["dr"]
    Dimension = data["Dimension"]
    etalist = data["etalist"]
    IntegralList = data["IntegralList"]
    PAPSystemMatrix = data["PAPSystemMatrix"]
    EndOfYearSystemMatrix = data["EndOfYearSystemMatrix"]
    timetaken = data["timetaken"]


    plt.figure()
    plt.loglog(
        etalist,IntegralList,
        label='Min: %0.3f'%(etalist[IntegralList.argmin()]))
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("ETA")
    plt.ylabel("New Resistant")
    plt.savefig(str(args.directory) + "/MinimumPlot.png")
