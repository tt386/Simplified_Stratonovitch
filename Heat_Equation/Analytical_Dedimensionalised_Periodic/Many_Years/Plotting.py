import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

from matplotlib.pyplot import cm

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
RTotMatrixMatrix = []

#Find list of all the datafiles
tempdirlist = os.listdir(args.directory)
dirlist = []
for i in tempdirlist:
    if os.path.isdir(os.path.join(args.directory,i)):
        dirlist.append(os.path.join(args.directory,i))

etaList = 0
timetaken = 0
PAP = 0
rhoList = []

Phi = 0
dx=0

SystemMatrixMatrix = []
PAPSystemMatrixMatrix = []
RSystemMatrixMatrix = []
for i in dirlist:
    with np.load(os.path.join(i, filename)) as data:
        RTotMatrixMatrix.append(data['RTotMatrix'])

        SystemMatrixMatrix.append(data["SystemMatrix"])
        PAPSystemMatrixMatrix.append(data["PAPSystemMatrix"])
        RSystemMatrixMatrix.append(data["RSystemMatrix"])

        rhoList.append(float(data['rho']))
        etaList = data["etaList"]
        timetaken = data["timetaken"]
        PAP = data["PAP"]
        Phi = data["Phi"]
        dx = data["dx"]

    

#Sort Lists correctly
zipped_lists = zip(rhoList, RTotMatrixMatrix, SystemMatrixMatrix, PAPSystemMatrixMatrix,  RSystemMatrixMatrix)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
rhoList, RTotMatrixMatrix, SystemMatrixMatrix, PAPSystemMatrixMatrix, RSystemMatrixMatrix = [ list(tuple) for tuple in  tuples]

rhoList = np.asarray(rhoList)
RTotMatrixMatrix = np.asarray(RTotMatrixMatrix)



for i in range(len(RTotMatrixMatrix)):  #For each rho
    rho = rhoList[i]
    print("Plotting Rho:",rho)
    color = cm.coolwarm(np.linspace(0,1,len(RTotMatrixMatrix[i])))

    NumYearLists = []
    for j in range(len(RTotMatrixMatrix[i][0])):
        NumYearLists.append([])

    fig, ax = plt.subplots()

    for j in range(len(RTotMatrixMatrix[i])): #For each eta
        TotList = RTotMatrixMatrix[i][j]
        YearList = np.arange(len(TotList))

        for Y in range(len(RTotMatrixMatrix[i][j])):
            NumYearLists[Y].append(TotList[Y])

        #FirstYearList.append(TotList[1])
        #LastYearList.append(TotList[-1])

        c = color[j]    

        plt.semilogy(YearList,TotList,c=c,linewidth=0.5)

    
    #fig.colorbar(color,orientation='vertical')
    plt.grid()
    plt.xlabel("Year")
    plt.ylabel("Integrated nummber of Resistant Pests")
    plt.title("RHO = %0.3f"%(rho))

    plt.savefig(str(args.directory) + "/ResistancePerYear_RHO_"+
        '{:.3E}'.format(rho) + ".png")
    plt.close()
        
    
    plt.figure()

    NumYearLists = np.asarray(NumYearLists)

    color = cm.Reds(np.linspace(0,1,len(NumYearLists)))

    MinEtaList = []
    MinResistantList = []

    MinSystem = 0
    MinPAPSystem = 0
    MinRSystem = 0

    for Y in range(len(NumYearLists)):
        c = color[Y]

        plt.loglog(etaList,NumYearLists[Y],c=c)

        min_index = np.argmin(NumYearLists[Y])

        if Y == 1:  #Remember 0 is the initial system
            MinSystem = SystemMatrixMatrix[i][min_index][1]
            MinPAPSystem = PAPSystemMatrixMatrix[i][min_index][1]
            print(MinSystem)
            print(len(MinSystem))
            MinRSystem = RSystemMatrixMatrix[i][min_index][1]

        MinEtaList.append(etaList[min_index])
        MinResistantList.append(NumYearLists[Y][min_index])

        plt.scatter(etaList[min_index],NumYearLists[Y][min_index],c='k')

    #plt.loglog(etaList,FirstYearList,label='First Year')
    #plt.loglog(etaList,LastYearList,label='Last Year')
    
    #plt.legend(loc='upper left')

    plt.grid()
    plt.xlabel("ETA")
    plt.ylabel("Number of Resistant")

    plt.title("RHO = %0.3f"%(rho))

    plt.savefig(str(args.directory) + "/MinimaCurves_RHO_"+
        '{:.3E}'.format(rho) + ".png")
    plt.close()


    ###########################################################
    #System At Minimised Eta
    plt.figure()
    plt.plot(np.linspace(0,1,len(MinSystem)),MinSystem,label='Susceptible Year')
    plt.plot(np.linspace(0,1,len(MinSystem)),MinPAPSystem,label='Susceptible PAP')
    plt.plot(np.linspace(0,1,len(MinRSystem)),MinRSystem,label='Resistant')

    plt.grid()

    plt.legend(loc='upper right')

    plt.xlabel("Space")
    plt.ylabel("Number")

    plt.title("RHO: %0.3f, ETA: "%(rho) + '{:.3E}'.format(MinEtaList[1]))

    plt.savefig(str(args.directory) + "SystemAtMinEta_RHO_"+
        '{:.3E}'.format(rho) + ".png")
    plt.close()
    #####################
    plt.figure()
    plt.semilogy(np.arange(len(MinEtaList))[1:],MinEtaList[1:])

    plt.grid()

    plt.xlabel("Year")
    plt.ylabel("Minimum ETA")

    plt.title("RHO = %0.3f"%(rho))

    plt.savefig(str(args.directory) + "/YearlyEtaMinimaCurves_RHO_"+
        '{:.3E}'.format(rho) + ".png")
    plt.close()

    ######################
    plt.figure()
    plt.semilogy(np.arange(len(MinEtaList)),MinResistantList)

    plt.grid()

    plt.xlabel("Year")
    plt.ylabel("Minimum Resistant Num")

    plt.title("RHO = %0.3f"%(rho))

    plt.savefig(str(args.directory) + "/YearlyResistantMinimaCurves_RHO_"+
        '{:.3E}'.format(rho) + ".png")
    plt.close()




for Y in range(len(RTotMatrixMatrix[0][0])):    #For each year
    print("Year:",Y)
    plt.figure()

    color = cm.Reds(np.linspace(0,1,len(RTotMatrixMatrix)))

    for i in range(len(RTotMatrixMatrix)):
        rho = rhoList[i]
        TotList = []

        c = color[i]

        for j in range(len(RTotMatrixMatrix[i])):
            TotList.append(RTotMatrixMatrix[i][j][Y])

        plt.loglog(1/np.sqrt(etaList),TotList/np.sqrt(etaList),label="RHO %0.2f"%(rho),c=c)


    plt.legend(bbox_to_anchor=(1.04,1),loc='upper left')
    plt.grid()
    plt.xlabel("SYSTEM SIZE")   #Used to eb ETA
    plt.ylabel("Resistant Number")

    plt.title("YEAR = %d"%(Y))

    plt.savefig(str(args.directory) + "/Year_" +
        str(Y).zfill(3) + ".png",bbox_inches="tight")
    plt.close()
