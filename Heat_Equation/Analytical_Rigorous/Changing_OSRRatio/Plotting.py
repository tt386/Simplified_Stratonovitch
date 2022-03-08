import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
import time

from scipy.stats import linregress

from scipy.optimize import curve_fit

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
    SlopeMatrix = data['SlopeMatrix']
    #RawMatrix = data['RawMatrix']
    OSRRatioList = data['OSRRatioList']
    k = data["k"]
    SystemSizeList = data["SystemSizeList"]
    timetaken = data["timetaken"]
    PAP = data["PAP"]
#SlopeMatrix = np.zeros((len(CharacteristicList),len(SystemSizeList)))

MinimumList = []

#Minima Slopes:
for c in range(len(SlopeMatrix)):
    OSR_Proportion = OSRRatioList[c]
    print("Starting minima slopes for k=",OSR_Proportion)
    slopelist = SlopeMatrix[c]


    plt.figure()
    plt.semilogy(SystemSizeList,slopelist)

    plt.xlabel("System Size")
    plt.ylabel("Log of initial yearly R Allele Slope")

    plt.title("OSR Proportion " + '{:.1E}'.format(OSR_Proportion))
    plt.grid()
    plt.savefig(str(args.directory) +
        "/Minima_Curve_OSRProportion_" + '{:.1E}'.format(OSR_Proportion) + ".png")
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
def fit_func(x,b,c,d):
    return  b *( (1*MAGNTIDUEINCREASE-x)**c * x**d)

def fit_func2(x,a,b,c,d,e):
    return a*(MAGNTIDUEINCREASE-x)**b + c*x**d + e

"""
def fit_func3(x,a,b,c):
    return a * ()
"""

def fit_func4(x,a):
    return a*( (1*MAGNTIDUEINCREASE-x)**-0.5 * x**-0.5)

firstplateau = 0
lastplateau = 0

for i in range(len(OSRRatioList)):
    if OSRRatioList[i] > 0.01:
        firstplateau = i
        break

for i in range(len(OSRRatioList)):
    if OSRRatioList[i] > 0.98:
        lastplateau = i
        break
    

fitMinimumList = np.asarray(MinimumList[firstplateau:lastplateau])*MAGNTIDUEINCREASE
fitOSRRatioList = OSRRatioList[firstplateau:lastplateau]*MAGNTIDUEINCREASE

sigmaweights = np.ones(len(fitOSRRatioList))
sigmaweights[:20] = 0.01
sigmaweights[20:] = 0.01
popt,pcov = curve_fit(fit_func,fitOSRRatioList,fitMinimumList,maxfev=100000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
b,c,d = popt


fitOSRRatioList = fitOSRRatioList/MAGNTIDUEINCREASE
fittedmin1 =  b*( (1 - fitOSRRatioList)**c * fitOSRRatioList**d) * MAGNTIDUEINCREASE**(c+d-1)

print("DifferenceList: ",fittedmin1 - fitMinimumList/MAGNTIDUEINCREASE)

#Second fit
fitOSRRatioList*=MAGNTIDUEINCREASE
popt,pcov = curve_fit(fit_func2,fitOSRRatioList,fitMinimumList,maxfev=100000000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
a,b,c,d,e = popt

fitOSRRatioList = fitOSRRatioList/MAGNTIDUEINCREASE
fittedmin2 = a * MAGNTIDUEINCREASE**(b-1) * (1-fitOSRRatioList)**b + c*MAGNTIDUEINCREASE**(d-1) * fitOSRRatioList**d + e*MAGNTIDUEINCREASE**-1



#4th fit
print("Starting 4th fit")
fitOSRRatioList*=MAGNTIDUEINCREASE
popt,pcov = curve_fit(fit_func4,fitOSRRatioList,fitMinimumList,maxfev=100000000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
a = popt

fitOSRRatioList = fitOSRRatioList/MAGNTIDUEINCREASE
fittedmin4 = a*( (1 - fitOSRRatioList)**-0.5 * fitOSRRatioList**-0.5) * MAGNTIDUEINCREASE**(-2)






plt.figure()
plt.loglog(OSRRatioList,MinimumList,label='data')#,label="Slope = %0.3f"%(result.slope))
plt.loglog(fitOSRRatioList,fittedmin1,label='fit: product')
#plt.loglog(fitOSRRatioList,fittedmin2,label='fit: sum')
plt.loglog(fitOSRRatioList,fittedmin4,label='fit: Theory')
plt.legend(loc='lower left')

plt.xlabel("OSR Proportion")
plt.ylabel("System Size at Minimum")

plt.title("Characteristic: " + '{0:.1f}'.format(k)+ ", PAP: %0.5f"%(PAP))
plt.grid()
plt.savefig(str(args.directory) + "/ComparingMinima.png")
plt.close()

endtime = time.time()

print("Time taken,",endtime-starttime)
