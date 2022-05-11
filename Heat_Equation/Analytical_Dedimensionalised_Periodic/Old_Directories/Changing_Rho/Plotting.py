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
rhoList = []

for i in dirlist:
    with np.load(os.path.join(i, filename)) as data:
        SlopeMatrix.append(data['SlopeList'])
        minmatrix.append(data['minlist'])
        rhoList.append(float(data['rho']))
        etaList = data["etaList"]
        timetaken = data["timetaken"]
        PAP = data["PAP"]
#SlopeMatrix = np.zeros((len(CharacteristicList),len(SystemSizeList)))

#Correctly order SlopeMatrxi and OSRProportion
"""
rhoList,SlopeMatrix,minmatrix = zip(*sorted(zip(rhoList,SlopeMatrix,minmatrix)))
"""
zipped_lists = zip(rhoList, SlopeMatrix,minmatrix)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
rhoList, SlopeMatrix, minmatrix= [ list(tuple) for tuple in  tuples]

rhoList = np.asarray(rhoList)
SlopeMatrix = np.asarray(SlopeMatrix)
minmatrix = np.asarray(minmatrix)
print(rhoList)

#sys.exit()


MinimumList = []

max_minpointlist = []

#Minima Slopes:
for c in range(len(SlopeMatrix)):
    Refuge_Proportion = rhoList[c]
    print("Starting minima slopes for k=",Refuge_Proportion)
    slopelist = SlopeMatrix[c]


    plt.figure()
    plt.semilogy(etaList,slopelist,label="Slope")
    plt.semilogy(etaList,minmatrix[c],label="Min Point")

    #Plot minima and maxima
    val, idx = max((val, idx) for (idx, val) in enumerate(minmatrix[c]))
    plt.semilogy(etaList[idx],val,"ro")

    max_minpointlist.append(etaList[idx])

    val, idx = min((val, idx) for (idx, val) in enumerate(slopelist))
    plt.semilogy(etaList[idx],val,"bo")

    plt.legend(loc='upper right')

    plt.xlabel("eta")
    plt.ylabel("Log of initial yearly R Allele Slope")

    plt.title("Refuge Proportion " + '{:.1E}'.format(Refuge_Proportion))
    plt.grid()
    plt.savefig(str(args.directory) +
        "/Minima_Curve_rho_" + '{:.1E}'.format(Refuge_Proportion) + ".png")
    plt.close()

    #Find the minumum of the curve
    minRAllele = 10000000
    minSystemSize = -1
    for i in range(len(slopelist)):
        if slopelist[i] < minRAllele:
            minRAllele = slopelist[i]
            minSystemSize = etaList[i]

    MinimumList.append(minSystemSize)
    print("minRAllele:",minRAllele)
    print("minSystemSize",minSystemSize)



plt.figure()
plt.plot(MinimumList,max_minpointlist,"ro")
plt.xlabel("eta at min slope")
plt.ylabel("eta at max minpoint")
plt.grid()
plt.savefig(str(args.directory) + "/Correlation.png")

plt.close()

"""
#Fitting
MAGNTIDUEINCREASE = 10**10
def fit_func(x,b,c,d):
    return  b *( (1*MAGNTIDUEINCREASE-x)**c * x**d)

def fit_func2(x,a,b,c,d,e):
    return a*(MAGNTIDUEINCREASE-x)**b + c*x**d + e


def fit_func3(x,a,b,c):
    return a * ()


firstplateau = 0
lastplateau = 0

for i in range(len(rhoList)):
    if rhoList[i] > 0.01:
        firstplateau = i
        break

for i in range(len(rhoList)):
    if rhoList[i] > 0.98:
        lastplateau = i
        break
    

fitMinimumList = np.asarray(MinimumList[firstplateau:lastplateau])*MAGNTIDUEINCREASE
fitrhoList = rhoList[firstplateau:lastplateau]*MAGNTIDUEINCREASE

sigmaweights = np.ones(len(fitrhoList))
#sigmaweights[:20] = 0.01
#sigmaweights[20:] = 0.01
sigmaweights*=0.01
popt,pcov = curve_fit(fit_func,fitrhoList,fitMinimumList,maxfev=100000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
b,c,d = popt


fitrhoList = fitrhoList/MAGNTIDUEINCREASE
fittedmin1 =  b*( (1 - fitrhoList)**c * fitrhoList**d) * MAGNTIDUEINCREASE**(c+d-1)

print("DifferenceList: ",fittedmin1 - fitMinimumList/MAGNTIDUEINCREASE)

#Second fit
fitrhoList*=MAGNTIDUEINCREASE
popt,pcov = curve_fit(fit_func2,fitrhoList,fitMinimumList,maxfev=100000000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
a,b,c,d,e = popt

fitrhoList = fitrhoList/MAGNTIDUEINCREASE
fittedmin2 = a * MAGNTIDUEINCREASE**(b-1) * (1-fitrhoList)**b + c*MAGNTIDUEINCREASE**(d-1) * fitrhoList**d + e*MAGNTIDUEINCREASE**-1
"""

t = 1-PAP

testtheorylist = (1-np.sqrt(1+(np.pi/rhoList)**2 *PAP* (1-rhoList)**2 / t))/(-4*(np.pi/rhoList)**2 * PAP)

plt.figure()
plt.loglog(rhoList,MinimumList,label='data')#,label="Slope = %0.3f"%(result.slope))
plt.loglog(rhoList,testtheorylist,label='Theory')
#plt.loglog(fitrhoList,fittedmin1,label='fit: product')
#plt.loglog(fitOSRRatioList,fittedmin2,label='fit: sum')
plt.legend(loc='lower left')

plt.xlabel("Refuge Proportion")
plt.ylabel("Eta at Minimum")

plt.title("PAP: %0.5f"%(PAP))
plt.grid()
plt.savefig(str(args.directory) + "/ComparingMinima.png")
plt.close()

endtime = time.time()

print("Time taken,",endtime-starttime)
