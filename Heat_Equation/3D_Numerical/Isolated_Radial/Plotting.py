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
    #DimensionList = data["DimensionList"]
    etalist = data["etalist"]
    IntegralList = data["IntegralList"]
    PAPSystemMatrix = data["PAPSystemMatrix"]
    EndOfYearSystemMatrix = data["EndOfYearSystemMatrix"]
    timetaken = data["timetaken"]
    print("TIMETAKEN:",timetaken)
    #########################################################################

    plt.figure()
    plt.loglog(
        etalist,IntegralList,
        label='Min: %0.3f'%(etalist[IntegralList.argmin()]))
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("ETA")
    plt.ylabel("New Resistant")
    plt.savefig(str(args.directory) + "/MinimumPlot.png")
    plt.close()

    #########################################################################

    plt.figure()
    plt.loglog(
        1/np.sqrt(etalist),IntegralList,
        label='Min: %0.3f'%((1/np.sqrt(etalist))[IntegralList.argmin()]))
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("L")
    plt.ylabel("New Resistant")
    plt.savefig(str(args.directory) + "/MinimumPlot_OSRLength.png")
    plt.close()

    print("L:",list(1/np.sqrt(etalist)))
    print("Slopelist:",list(IntegralList))

    #########################################################################
    
    plt.figure()
    plt.loglog(
        np.pi/etalist,IntegralList,
        label='Min: %0.3f'%((np.pi/(etalist))[IntegralList.argmin()]))
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("L")
    plt.ylabel("New Resistant")
    plt.savefig(str(args.directory) + "/MinimumPlot_OSRArea.png")
    plt.close()


    #########################################################################


    ThirdOrder = abs(np.gradient(np.gradient(np.gradient(IntegralList,etalist),etalist),etalist))
    plt.figure()
    plt.loglog(
        etalist,ThirdOrder)
    #plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel("ETA")
    plt.ylabel("New Resistant Third Order")
    plt.savefig(str(args.directory) + "/ThirdOrderPlot.png")
    plt.close()


    for i in range(len(PAPSystemMatrix)):
        print("Image:",i)
        eta = etalist[i]
        PAPSystem = PAPSystemMatrix[i]
        EndOfYearSystem = EndOfYearSystemMatrix[i]
        Integrand = Phi/(Phi+EndOfYearSystem)

        Dimension = DimensionList[i]

        print(PAPSystem)

        plt.figure()  
        
        xlist = np.arange(int(Dimension))*dr
 
        plt.semilogy(xlist[xlist<20],PAPSystem[xlist<20],label='PAP')
        plt.semilogy(xlist[xlist<20],EndOfYearSystem[xlist<20],label="EOY")
        plt.semilogy(xlist[xlist<20],Integrand[xlist<20],label='Integrand')
    

        plt.xlim(0,10)
        plt.ylim(1e-6,10)   
        plt.grid() 

        plt.ylabel("Susceptible Pest Dist")
        plt.xlabel("Radial Space")

        plt.legend(loc='lower right')

        plt.title("Eta: " + '{:.1E}'.format(eta) )
        plt.savefig(
            str(args.directory) + "/Frame_" + str(i).zfill(5) + '.png')
        plt.close()        


        """
        plt.figure()
        plt.imshow(EndOfYearSystem, cmap='hot_r', vmin=0, vmax=1,interpolation='nearest')
        plt.colorbar()
        plt.title("EOYDist, Eta: " + '{:.1E}'.format(eta) )
        plt.savefig(
            str(args.directory) + "/EOYSystem_" + str(i).zfill(5) + '.png')
        plt.close()
        """
