import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
import time

import scipy

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




timetaken = 0
PAP = 0
rhoList = []

Phi = 0
dx=0

OSRSeparation = 0
LeftMostMinima = []
RightMostMinima = []

PAPSystemMatrixList = []
EndOfYearSystemMatrixList = []

for i in dirlist:
    ######################
    ###Downloading Data###
    ######################
    with np.load(os.path.join(i, filename)) as data:
        SlopeMatrix.append(data['IntegralList'])
        etalist = data["etalist"]
        timetaken = data["timetaken"]
        PAP = data["PAP"]
        Phi = data["Phi"]
        rhoList.append(data["rho"])
        PAPSystemMatrix = data["PAPSystemMatrix"]
        EndOfYearSystemMatrix = data["EndOfYearSystemMatrix"]

        print("Rho:",data["rho"],", Time:",timetaken)

    """
    rho = rhoList[-1]
    print("rho: ",rho)
    
    for j in range(len(EndOfYearSystemMatrix)):
        eta = etalist[j]
        print("eta",eta)
        #PAPSystem = PAPSystemMatrix[j]
        EndOfYearSystem = np.asarray(EndOfYearSystemMatrix[j])
        print(EndOfYearSystem)


        plt.figure()
        plt.imshow(EndOfYearSystem, cmap='hot_r', vmin=0, vmax=1,interpolation='nearest')
        plt.colorbar()
        plt.title("EndOfYearDist, Rho: %0.3f, Eta: "%(rho) + '{:.1E}'.format(eta) )
        plt.savefig(i + "/EndOfYearDist_" + str(j).zfill(6) + ".png")
        plt.close()    
    """
#Sort these lists
zipped_lists = zip(rhoList, SlopeMatrix)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
rhoList, SlopeMatrix = [ list(tuple) for tuple in  tuples]

rhoList = np.asarray(rhoList)
SlopeMatrix = np.asarray(SlopeMatrix)



"""
#############
###Video Images
#############
for c in range(len(PAPSystemMatrixList)):
    PAPSystemMatrix = PAPSystemMatrixList[c]
    EndOfYearSystemMatrix = EndOfYearSystemMatrixList[c]

    rho = rhoList[c]

    for i in range(len(PAPSystemMatrix)):
        eta = etalist[i]
        PAPSystem = PAPSystemMatrix[i]
        EndOfYearSystem = EndOfYearSystemMatrix[i]

        plt.figure()
        plt.imshow(PAPSystem, cmap='hot_r', vmin=0, vmax=1,interpolation='nearest')
        plt.colorbar()
        plt.title("PAPDist, Rho: %0.3f, Eta: "%(rho) + '{:.1E}'.format(eta) )
        plt.close()
"""



MinimaList = []

SecondOrderDerivList = []
SecondOrderDerivList_rho = []

SecondOrderDerivList_min = []
SecondOrderDerivList_rho_min = []

LowerKinkEta = []
UpperKinkEta = []

LowerKinkEta_rho = []
UpperKinkEta_rho = []

#Minima Slopes
for c in range(len(SlopeMatrix)):

    slopelist = SlopeMatrix[c]

    MinimaList.append(etalist[np.argmin(slopelist)])

    rho = rhoList[c]

    print("rho",rho)

    plt.figure()
    plt.loglog(etalist,slopelist/etalist)

    plt.xlabel("ETA")
    plt.ylabel("# New Resistant")

    plt.grid()

    plt.title("rho: %0.3f"%(rho))

    plt.savefig(str(args.directory) +
        "/Minima_Curve_Rho_" + '{:.3E}'.format(rho) + ".png")
    plt.close()

    #####################################################################
    #####################################################################
    #####################################################################

    plt.figure()
    plt.loglog(etalist,abs(np.gradient(slopelist,etalist)))
    plt.legend(loc='lower left')

    plt.title("Refuge Proportion " + '{:.1E}'.format(rho))

    plt.xlabel("eta")
    plt.ylabel("First Derivative")

    plt.grid()

    plt.savefig(str(args.directory) +
        "/FirstOrderDeriv_Curve_rho_" + '{:.3E}'.format(rho) + ".png")
    plt.close()

    #####################################################################
    #####################################################################
    #####################################################################

    #Second order derivative
    plt.figure()
    #plt.loglog(etalist,abs(np.gradient(np.gradient(slopelist))),label="Data")

    """
    Spline = scipy.interpolate.splrep(etalist,slopelist,k=3)
    SecondOrder = scipy.interpolate.splev(etalist,Spline,der=2)
    
    plt.loglog(etalist,SecondOrder,label="Data")
    """
    plt.loglog(etalist,abs(np.gradient(np.gradient(slopelist,etalist),etalist)),label="Data")


    plt.legend(loc='lower left')

    plt.title("Refuge Proportion " + '{:.1E}'.format(rho))

    plt.xlabel("eta")
    plt.ylabel("Second Derivative")

    plt.grid()

    plt.savefig(str(args.directory) +
        "/SecondOrderDeriv_Curve_rho_" + '{:.3E}'.format(rho) + ".png")
    plt.close()

    #####################################################################
    #####################################################################
    #####################################################################


    plt.figure()
    #plt.loglog(etalist,abs(np.gradient(np.gradient(slopelist))),label="Data")

    """
    Spline = scipy.interpolate.splrep(etalist,slopelist,k=3)
    SecondOrder = scipy.interpolate.splev(etalist,Spline,der=2)

    plt.loglog(etalist,SecondOrder,label="Data")
    """
    plt.loglog(etalist,abs(np.gradient(np.gradient(np.gradient(slopelist,etalist),etalist),etalist)),label="Data")


    plt.legend(loc='lower left')

    plt.title("Refuge Proportion " + '{:.1E}'.format(rho))

    plt.xlabel("eta")
    plt.ylabel("Third Derivative")

    plt.grid()

    plt.savefig(str(args.directory) +
        "/ThirdOrderDeriv_Curve_rho_" + '{:.3E}'.format(rho) + ".png")
    plt.close()

    #####################################################################
    #####################################################################
    #####################################################################

    plt.figure()
    plt.loglog(etalist,abs(np.gradient(np.gradient(np.gradient(np.gradient(slopelist,etalist),etalist),etalist),etalist)),label="Data")

    plt.legend(loc='lower left')

    plt.title("Refuge Proportion " + '{:.1E}'.format(rho))

    plt.xlabel("eta")
    plt.ylabel("Fourth Derivative")

    plt.grid()

    plt.savefig(str(args.directory) +
        "/FourthOrderDeriv_Curve_rho_" + '{:.3E}'.format(rho) + ".png")
    plt.close()




    index_min = np.argmin(slopelist)    #index of global minimum

    #Upper kink: look at miximum in 1st order
    FirstOrder = abs(np.gradient(slopelist,etalist))
    
    UpperIndex = 0
    CurrentMax = 0
    for i in range(index_min,len(slopelist)):
        if FirstOrder[i] > CurrentMax:
            UpperIndex = i
            CurrentMax = FirstOrder[i]
    UpperKinkEta.append(etalist[UpperIndex])
    UpperKinkEta_rho.append(rho)

    #Lower kink: look at maximum in 3rd order
    ThirdOrder = abs(np.gradient(np.gradient(np.gradient(slopelist,etalist),etalist),etalist))
    LowerIndex = 0
    CurrentMax = 0
    for i in range(index_min,1,-1):
        if (ThirdOrder[i+1] < ThirdOrder[i]) and (ThirdOrder[i-1]<ThirdOrder[i]):
            LowerKinkEta.append(etalist[i])
            LowerKinkEta_rho.append(rho)
            break

    """
     #finding max points
    if rho > 0:
        SecondOrderDerivMaxList = [argrelextrema(abs(np.gradient(np.gradient(slopelist,etalist),etalist)), np.greater)]
        #print("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
        #print(SecondOrderDerivMaxList)
        #print(SecondOrderDerivMaxList[0])
        print(SecondOrderDerivMaxList[0][0])
        #print("BBBBBBBBBBBBBBBBBBBBBBBBBBB")
        try:
            #print(SecondOrderDerivMaxList)
            #print(SecondOrderDerivMaxList[0][0][1])
            SecondOrderDerivList.append(etalist[SecondOrderDerivMaxList[0][0][1]])
            SecondOrderDerivList_rho.append(rho)
        except:
            print("no maxima found")

        #print("SecondOrderDerivList:",SecondOrderDerivList)


    #finding min points
    if rho > 0:
        secondderiv = abs(np.gradient(np.gradient(slopelist,etalist),etalist))
        SecondOrderDerivMinList = [argrelextrema(secondderiv[secondderiv>10e-10], np.less)]
        #print("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
        #print("Second order min:",SecondOrderDerivMinList)
        #print("Second order min:",SecondOrderDerivMinList[0])
        print("Second order min:",SecondOrderDerivMinList[0][0])
        #print("BBBBBBBBBBBBBBBBBBBBBBBBBBB")
        try:
            #print(SecondOrderDerivMinList)
            #print(SecondOrderDerivMinList[0][0][0])
            SecondOrderDerivList_min.append(etalist[SecondOrderDerivMinList[0][0][-1]])
            SecondOrderDerivList_rho_min.append(rho)
        except:
            print("no minima found")

    """


plt.figure()
plt.semilogy(rhoList,MinimaList)

plt.grid()
plt.xlabel("RHO")
plt.ylabel("ETA at minimum")

plt.title("How Minima Change with RHO")

plt.savefig(str(args.directory) + "/ChangingMinima.png")
plt.close()

#########################################################################
#########################################################################
#########################################################################
if len(LowerKinkEta_rho) > 0:
    x,y = 1-np.asarray(LowerKinkEta_rho),np.asarray(LowerKinkEta)

    theory = x**2 /(16*(1-PAP)* np.log(1/Phi))

    theory =  x**2 /(16*(1-PAP)) * 10

    slope, intercept, r,p,se = scipy.stats.linregress(np.log(x),np.log(y))

    plt.figure()
    plt.loglog(x,y,label="Data")
    plt.loglog(x,np.exp(intercept) * x**slope,label='Fit: slope: %0.5f, Int: %0.5f'%(slope,np.exp(intercept)))
    plt.loglog(x,theory,"--",label='Theory')
    plt.legend(loc='upper left')
    plt.grid()
    plt.xlabel("1-rho")
    plt.ylabel("eta at Lower Kink")
    plt.title("How Lower Kink Changed with rho")

    #plt.xlim(0.3,1)
    plt.savefig(str(args.directory) + "/LowerKink_ChangingRho.png")
    plt.close()
else:
    print("No Minima Found")
#########################################################################
#########################################################################
#########################################################################

x,y = np.asarray(UpperKinkEta_rho),np.asarray(UpperKinkEta)

slope, intercept, r,p,se = scipy.stats.linregress(x,np.log(y))

print("Slope:",slope)
print("Intercept:",intercept)

plt.figure()
plt.semilogy(x,y,label="Data")
plt.semilogy(x,np.exp(intercept + slope*x),label='Fit slope: %0.5f'%(slope))
#plt.loglog(x,theorykinks,"--",label='Theory')
plt.legend(loc='upper left')
plt.grid()
plt.xlabel("rho")
plt.ylabel("eta at Upper Kink")
plt.title("How Upper Kink changes with rho")

plt.savefig(str(args.directory) + "/UpperKink_ChangingRho.png")
plt.close()




"""
############################################################################
############################################################################
############################################################################
x = 1.-np.asarray(SecondOrderDerivList_rho)
y = np.asarray(SecondOrderDerivList)

x_slope = np.log(x[x>0.3])
y_slope = np.log(y[x>0.3])

slope, intercept, r, p, se = linregress(x_slope, y_slope)


theorykinks = (x)**2/(16*(1-PAP) * special.erfinv(1-Phi)**2)

plt.figure()
plt.loglog(x,y,label="Data")
plt.loglog(x,theorykinks,"--",label='Theory')
plt.legend(loc='upper left')
plt.grid()
plt.xlabel("1-rho")
plt.ylabel("eta at kink")
plt.title("How initial kink changes with rho")

plt.xlim(0.3,1)
plt.savefig(str(args.directory) + "/LowerKink_ChangingRho.png")
plt.close()

############################################################################
############################################################################
############################################################################


x = np.asarray(SecondOrderDerivList_rho_min)
y = np.asarray(SecondOrderDerivList_min)

x_slope = np.log(x)
y_slope = np.log(y)

slope, intercept, r, p, se = linregress(x_slope, y_slope)

theorykinks = (x**2/(np.pi**2 * PAP)) * np.log(8*x/(Phi*np.pi**2))

plt.figure()
plt.loglog(SecondOrderDerivList_rho_min,SecondOrderDerivList_min,label="Data")
plt.loglog(x,theorykinks,"--",label='Theory')
plt.legend(loc='upper left')
plt.grid()
plt.xlabel("rho")
plt.ylabel("eta at Second order min")
plt.title("How Second order min changes with rho")

plt.savefig(str(args.directory) + "/UpperKink_ChangingRho.png")
plt.close()
"""
