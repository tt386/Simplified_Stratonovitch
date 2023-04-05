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


HeatMap_SepList = []
HeatMap_NumList = []
HeatMap_SizeList = []
HeatMap_SlopeList = []

slopelist = []
for i in dirlist:
    ######################
    ###Downloading Data###
    ######################
    with np.load(os.path.join(i, filename)) as data:
        SlopeMatrix.append(data['RDifferenceList'])     #A list for each different separation lenght
        OSRNumList.append(data['OSRNum'])
        L = data["L"]
        eta = data["eta"]
        timetaken = data["timetaken"]
        print("TIME TAKEN:",timetaken)
        PAP = data["PAP"]
        Phi = data["Phi"]
        OSRSeparationList = data["OSRSeparationList"]
        
        xlist = data['xlist']
        PAPMatrix = data["PAPMatrix"]
        YearMatrix = data["YearMatrix"]
        
        timetaken = data['timetaken']
        
        for j in range(len(OSRSeparationList)):    
            HeatMap_SepList.append(OSRSeparationList[j])
            HeatMap_NumList.append(OSRNumList[-1])
            HeatMap_SizeList.append(L/OSRNumList[-1])
            HeatMap_SlopeList.append(SlopeMatrix[-1][j])

    slopelist.append(SlopeMatrix[-1][0])

    for j in range(len(OSRSeparationList)):
        print("OSR Separation:",OSRSeparationList[j])
        Sep = OSRSeparationList[j]
        
        PAPList = PAPMatrix[j]
        YearList= YearMatrix[j]

        ##########################
        plt.figure()
        plt.semilogy(xlist,PAPList,label='PAP')
        plt.semilogy(xlist,YearList,label='Year')
        plt.semilogy(xlist,Phi/(Phi+YearList)-Phi,label='Integrand: New R: ' + '{:.3E}'.format(SlopeMatrix[-1][j]))

        plt.title("Sep Length %0.3f"%(Sep))
        plt.xlim(-1,2)


        plt.legend(loc='lower right')    

        #plt.xlim(0,(L+Sep*OSRNumList[-1])*2)

        plt.grid()
        plt.xlabel("Space")
        plt.ylabel("Number")

        plt.savefig(i + "/System_Sep_%0.3f.png"%(Sep))
        plt.close()
    ############################


#Sort these lists
zipped_lists = zip(OSRNumList, SlopeMatrix)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
OSRNumList, SlopeMatrix= [ list(tuple) for tuple in  tuples]

OSRNumList = np.asarray(OSRNumList)
SlopeMatrix = np.asarray(SlopeMatrix)


for i in range(len(SlopeMatrix)):   #For each OSR Num
    OSRNum = OSRNumList[i]

    print("OSRNum,",OSRNum)

    slopelist = SlopeMatrix[i]

    plt.figure()
    plt.loglog(OSRSeparationList,slopelist)
    plt.grid()

    plt.title("OSR Number: %0.3f"%(OSRNum))

    plt.xlabel("OSR Separation")
    plt.ylabel("Number new R")

    plt.savefig(str(args.directory) + "/MinPlot_OSRNum_" + '{:.3E}'.format(OSRNum) + ".png")
    plt.close()


#############################################################################
#############################################################################
#############################################################################
SepValues = list(set(HeatMap_SepList))
SepValues.sort()
Sep_List = []
for i in SepValues:
    Sep_List.append([i,[]])

for i in range(len(HeatMap_SizeList)):
    for j in range(len(Sep_List)):
        if HeatMap_SepList[i] == Sep_List[j][0]:
            Sep_List[j][1].append([HeatMap_SizeList[i],HeatMap_SlopeList[i]])
            break

for i in range(len(Sep_List)):
    Sep = Sep_List[i][0]
    print("SEP:",Sep)

    unsortedarray = np.asarray(Sep_List[i][1])
    print(unsortedarray)


    sortedarray = unsortedarray[unsortedarray[:, 0].argsort()]

    plt.figure()
    plt.loglog(sortedarray[:,0],sortedarray[:,1])

    plt.grid()
    plt.xlabel("Size of SubFields")
    plt.ylabel("Number of New R Pests")

    plt.title("Separation: %0.5f"%(Sep))

    plt.savefig(str(args.directory) + '/MinPlot_Sep_' + '{:.3E}'.format(Sep) + ".png")
    plt.close()




logHeatMap_SepList = np.log10(np.asarray(HeatMap_SepList))
logHeatMap_NumList = np.log10(np.asarray(HeatMap_NumList))
logHeatMap_SizeList = np.log10(np.asarray(HeatMap_SizeList))
logHeatMap_SlopeList = np.log10(np.asarray(HeatMap_SlopeList))


fig = plt.figure(2)
ax = fig.add_subplot(111)

cp = ax.tricontourf(logHeatMap_SepList.ravel(), logHeatMap_SizeList, logHeatMap_SlopeList.ravel(),10,cmap='coolwarm')

#fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.scatter(logHeatMap_SepList,logHeatMap_SizeList,c=logHeatMap_SlopeList,s=150,cmap='coolwarm')


C = 10*L
ConstrainSepList = np.linspace(min(HeatMap_SepList),max(HeatMap_SepList),100)
ConstraintSize = ConstrainSepList*L/(C+ConstrainSepList+L)
 

plt.plot(np.log10(ConstrainSepList),np.log10(ConstraintSize),'--k',label="C: %0.3f"%(C),linewidth=3)
plt.ylim(min(logHeatMap_SizeList),max(logHeatMap_SizeList))

#plt.legend(loc='lower left',fontsize=20)

ax.set_ylabel(r'$\log_{10}(\frac{L}{\sqrt{DY}}$)',fontsize=30)
ax.set_xlabel(r'$\log_{10}(\frac{\delta}{\sqrt{DY}}$)',fontsize=30)


cbar = plt.colorbar(cp)
cbar.set_label(label=r'$\log_{10}(\frac{\sqrt{DY}}{L}\Delta R$)',size=30)

cbar.ax.tick_params(labelsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)

plt.xticks(fontsize=20,rotation=45)
plt.yticks(fontsize=20,rotation=45)


plt.tight_layout()
plt.grid(True)

plt.savefig(str(args.directory) + '/Scatter.png')
plt.close()



minindex = logHeatMap_SlopeList.argmin()

print("Min Sep:",logHeatMap_SepList[minindex])
print("Min Size:",logHeatMap_SizeList[minindex])
print("Min Diff:",logHeatMap_SlopeList[minindex])
