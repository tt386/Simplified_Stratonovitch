import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
import time

from scipy.stats import linregress

from scipy.optimize import curve_fit

import sys

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
with np.load(os.path.join(args.directory, filename)) as data:
    SlopeList = data['SlopeList']
    #RawMatrix = data['RawMatrix']
    Year = data["Year"]
    k = data["k"]
    OSR_Proportion =data["OSR_Proportion"]
    SystemSizeList = data["SystemSizeList"]
    timetaken = data["timetaken"]
    PAP = data["PAP"]
    AnalyticalList = data["AnalyticalList"]
    NumericalList = data["NumericalList"]
    dx = data["dx"]
    dt = data["dt"]


for i in range(len(SystemSizeList)):
    print("Plotting SystemSize:",SystemSizeList[i])
    ylist = AnalyticalList[i][0]
    yprimelist = AnalyticalList[i][1]

    EndOfPAPDomain = NumericalList[i][0]
    S_Domain = NumericalList[i][1]


    xlist = np.arange(len(ylist))*dx 

    
    fig, axs = plt.subplots(2)
    fig.suptitle("System Size %0.5f"%(SystemSizeList[i]))

    axs[0].semilogy(xlist,ylist,'--b',label='After PAP: Analytical')
    #axs[0].semilogy(np.arange(len(S_Domain))*dx,EndOfPAPDomain,'--r',label='After PAP: Numerical')

    axs[0].semilogy(xlist,yprimelist,'b',label='End of Year: Analytical')
    #axs[0].semilogy(np.arange(len(S_Domain))*dx,S_Domain,'r',label='End of Year: Numerical')
    #axs[0].legend(loc='upper center')
    axs[0].grid()
    axs[0].set_ylim(bottom=1e-5,top=1)

    axs[1].plot(SystemSizeList,SlopeList)
    axs[1].plot(SystemSizeList[i],SlopeList[i],'ro')
    axs[1].grid()

    plt.savefig(args.directory + "/SystemSizeNum_" + str(i).zfill(5))






