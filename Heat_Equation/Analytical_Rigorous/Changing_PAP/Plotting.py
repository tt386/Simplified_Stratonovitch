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
SlopeMatrix = []
PAPList = []

#Find list of all the datafiles
tempdirlist = os.listdir(args.directory)
dirlist = []
for i in tempdirlist:
    if os.path.isdir(os.path.join(args.directory,i)):
        dirlist.append(os.path.join(args.directory,i))



k = 0
SystemSizeList = 0
timetaken = 0
OSR_Proportion = 0
Year = 0
for i in dirlist:
    try:
        with np.load(os.path.join(i, filename)) as data:
            SlopeMatrix.append(data['SlopeList'])
            #RawMatrix = data['RawMatrix']
            PAPList.append(float(data['PAP']))
            k = data["k"]
            OSR_Proportion =data["OSR_Proportion"]
            SystemSizeList = data["SystemSizeList"]
            timetaken = data["timetaken"]
            PAP = data["PAP"]
            Year = data["Year"]
    except:
        pass
#SlopeMatrix = np.zeros((len(CharacteristicList),len(SystemSizeList)))
print("PAPList",PAPList)
print("SystemSizeList",SystemSizeList)
#Correctly order SlopeMatrxi and OSRProportion
zipped_lists = zip(PAPList, SlopeMatrix)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
PAPList, SlopeMatrix = [ list(tuple) for tuple in  tuples]

PAPList = np.asarray(PAPList)
SlopeMatrix = np.asarray(SlopeMatrix)

print(PAPList)

#sys.exit()


MinimumList = []

#Minima Slopes:
for c in range(len(SlopeMatrix)):
    PAP = PAPList[c]
    print("Starting minima slopes for PAP=",PAP)
    slopelist = SlopeMatrix[c]


    plt.figure()
    plt.semilogy(SystemSizeList,slopelist)

    plt.xlabel("System Size")
    plt.ylabel("Log of initial yearly R Allele Slope")

    plt.title("PAP " + '{:.1E}'.format(PAP))
    plt.grid()
    plt.savefig(str(args.directory) +
        "/Minima_Curve_PAP_" + '{:.1E}'.format(PAP) + ".png")
    plt.close()

    #Find the minumum of the curve
    minRAllele = 10000000
    minSystemSize = -1
    for i in range(len(slopelist)):
        if slopelist[i] < minRAllele:
            minRAllele = slopelist[i]
            minSystemSize = SystemSizeList[i]

    MinimumList.append(minSystemSize)
    print("minRAllele:",minRAllele)
    print("minSystemSize",minSystemSize)



#Fitting
MAGNTIDUEINCREASE = 10**10
def fit_func(x,a):
    return  a * (x*(MAGNTIDUEINCREASE-x))**0.25

firstplateau = 0
lastplateau = 0

PAPList/=Year

for i in range(len(PAPList)):
    if PAPList[i] > 0.0:
        firstplateau = i
        break

for i in range(len(PAPList)):
    if PAPList[i] > 0.4:
        lastplateau = i
        break
    

fitMinimumList = np.asarray(MinimumList[firstplateau:lastplateau])*MAGNTIDUEINCREASE
fitPAPList = PAPList[firstplateau:lastplateau]*MAGNTIDUEINCREASE

sigmaweights = np.ones(len(fitPAPList))
#sigmaweights[:20] = 0.01
#sigmaweights[20:] = 0.01
sigmaweights*=0.0001

print("fitMinimumList",fitMinimumList)
print("fitPAPList",fitPAPList)
popt,pcov = curve_fit(fit_func,fitPAPList,fitMinimumList,maxfev=100000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
a = popt


fitPAPList = fitPAPList/MAGNTIDUEINCREASE
fittedmin1 =  a * (fitPAPList*(1-fitPAPList))**0.25 * MAGNTIDUEINCREASE**(0.5-1)

print("DifferenceList: ",fittedmin1 - fitMinimumList/MAGNTIDUEINCREASE)
"""
#Second fit
fitPAPList*=MAGNTIDUEINCREASE
popt,pcov = curve_fit(fit_func2,fitPAPList,fitMinimumList,maxfev=100000000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
a,b,c,d,e = popt

fitPAPList = fitPAPList/MAGNTIDUEINCREASE
fittedmin2 = a * MAGNTIDUEINCREASE**(b-1) * (1-fitPAPList)**b + c*MAGNTIDUEINCREASE**(d-1) * fitPAPList**d + e*MAGNTIDUEINCREASE**-1
"""

plt.figure()
plt.semilogy(PAPList,MinimumList,label='data')#,label="Slope = %0.3f"%(result.slope))
plt.loglog(fitPAPList,fittedmin1,label='fit: theory')
#plt.loglog(fitPAPList,fittedmin2,label='fit: sum')
plt.legend(loc='lower left')

plt.xlabel("PAP")
plt.ylabel("System Size at Minimum")

plt.title("Characteristic: " + '{0:.1E}'.format(k)+ ", OSR Proportion: %0.5f"%(OSR_Proportion))
plt.grid()
plt.savefig(str(args.directory) + "/ComparingMinima.png")
plt.close()

endtime = time.time()

print("Time taken,",endtime-starttime)
