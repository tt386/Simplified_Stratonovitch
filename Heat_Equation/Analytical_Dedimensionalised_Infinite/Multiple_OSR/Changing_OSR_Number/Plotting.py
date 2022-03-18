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
SlopeMatrix = []
OSRRatioList = []
minmatrix = []

#Find list of all the datafiles
tempdirlist = os.listdir(args.directory)
dirlist = []
for i in tempdirlist:
    if os.path.isdir(os.path.join(args.directory,i)):
        dirlist.append(os.path.join(args.directory,i))



etaList = 0
timetaken = 0
PAP = 0
OSRNumList = []

Phi = 0
dx=0

OSRSeparation = 0
LeftMostMinima = []
RightMostMinima = []

for i in dirlist:
    ######################
    ###Downloading Data###
    ######################
    with np.load(os.path.join(i, filename)) as data:
        SlopeMatrix.append(data['RDifferenceList'])
        OSRNumList.append(data['OSRNum'])
        etalist = data["etalist"]
        timetaken = data["timetaken"]
        PAP = data["PAP"]
        Phi = data["Phi"]
        OSRSeparation = data["OSRSeparation"]
        
        xlist = data['xlist']
        PAPMatrix = data["PAPMatrix"]
        YearMatrix = data["YearMatrix"]
        
        timetaken = data['timetaken']


    ################
    ###Find Minima##
    ################
    """    
    minlistindices = argrelextrema(SlopeMatrix[-1],np.less)
    minlist = etalist[minlistindices[0]]
    """
    minlist = []
    
    for j in range(1,len(SlopeMatrix[-1])-1):
        if  (SlopeMatrix[-1][j-1] > SlopeMatrix[-1][j] and 
            SlopeMatrix[-1][j] < SlopeMatrix[-1][j+1]):
            minlist.append(etalist[j])
    print(minlist) 

    if len(minlist)==0: 
        minlist.append(0)

    LeftMostMinima.append(minlist[0])
    RightMostMinima.append(minlist[-1])
    
    rlist = np.ones(len(xlist))*(Phi)
    print("Plotting Sep:",OSRNumList[-1]," Time:",timetaken)
    Plotting = [0,5,-1]
    for p in Plotting:
        eta = etalist[p]
        PAPylist = PAPMatrix[p]
        Yearylist = YearMatrix[p]

        plt.figure()
        plt.semilogy(xlist,PAPylist,label='PAP')
        plt.semilogy(xlist,Yearylist,label='Year')
        plt.semilogy(xlist,rlist/(rlist+Yearylist)-rlist,label='Integrand')
        plt.legend(loc='upper left')
        plt.title("Diffusion Coeff: "+ '{0:.5f}'.format(eta))
        
        plt.xlabel("Space")
        plt.grid()
        plt.savefig(i + "/Eta_" + '{0:.5f}'.format(eta) + ".png")
        plt.close()    


#Sort these lists
zipped_lists = zip(OSRNumList, SlopeMatrix, LeftMostMinima, RightMostMinima)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
OSRNumList, SlopeMatrix, LeftMostMinima, RightMostMinima= [ list(tuple) for tuple in  tuples]

OSRNumList = np.asarray(OSRNumList)
SlopeMatrix = np.asarray(SlopeMatrix)



#Minima Slopes
for c in range(len(SlopeMatrix)):

    slopelist = SlopeMatrix[c]
    OSRNum = OSRNumList[c]

    plt.figure()
    plt.loglog(etalist,slopelist)


    #Add markers at minima points
    leftindex = np.where(etalist==LeftMostMinima[c])[0]
    rightindex = np.where(etalist==RightMostMinima[c])[0]
    print("leftindex",leftindex)
    if len(leftindex) > 0:
        plt.loglog(LeftMostMinima[c],slopelist[leftindex],"ko")
        plt.loglog(RightMostMinima[c],slopelist[rightindex],"ko")

    plt.xlabel("ETA")
    plt.ylabel("# New Resistant")

    plt.grid()

    plt.title("OSRNum: %d, OSRSep: %0.3f"%(OSRNum,OSRSeparation))

    plt.savefig(str(args.directory) + 
        "/Minima_Curve_Num_" + '{:.1E}'.format(OSRNum) + ".png")
    plt.savefig(str(args.directory) + "/Frame_" + str(c).zfill(5) + ".png")
    plt.close()



####################################
###Plotting Minima##################
####################################
RightMostMinimaTheory = ((RightMostMinima[0]/(OSRNumList)**2) 
    * (OSRNumList+(OSRNumList-1)*OSRSeparation)**2)

print("Leftmost: ",LeftMostMinima)
plt.figure()

plt.plot(OSRNumList,LeftMostMinima,label='Leftmost')
plt.plot(OSRNumList,RightMostMinima,label='Rightmost')

if OSRNumList[0] == 0:
    plt.plot(OSRNumList,RightMostMinimaTheory,"--",label='RightMost Theory')



"""
Rho     Prediction
0.5     0.164
0.2     0.0468
0.8     0.6
0.6     0.24375
"""

PREDICTED = 0.6
plt.plot(OSRNumList,np.ones(len(OSRNumList))*PREDICTED,'--',label='Infinite Theory')


plt.legend(loc='upper center')

plt.grid()

plt.xlabel("OSR Number")
plt.ylabel("ETA value of Minima")

plt.title("OSR Separation: %0.3f"%(OSRSeparation))

plt.savefig(str(args.directory) + "/Locations_Of_Minima.png")

plt.close()

