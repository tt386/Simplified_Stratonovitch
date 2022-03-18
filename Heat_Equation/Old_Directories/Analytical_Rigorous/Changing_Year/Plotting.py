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
YearList = []

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
for i in dirlist:
    try:
        with np.load(os.path.join(i, filename)) as data:
            SlopeMatrix.append(data['SlopeList'])
            #RawMatrix = data['RawMatrix']
            YearList.append(float(data['Year']))
            k = data["k"]
            OSR_Proportion =data["OSR_Proportion"]
            SystemSizeList = data["SystemSizeList"]
            timetaken = data["timetaken"]
            PAP = data["PAP"]
    except:
        pass
#SlopeMatrix = np.zeros((len(CharacteristicList),len(SystemSizeList)))
print("YearList",YearList)
print("SystemSizeList",SystemSizeList)
#Correctly order SlopeMatrxi and OSRProportion
zipped_lists = zip(YearList, SlopeMatrix)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
YearList, SlopeMatrix = [ list(tuple) for tuple in  tuples]

YearList = np.asarray(YearList)
SlopeMatrix = np.asarray(SlopeMatrix)

print(YearList)

#sys.exit()


MinimumList = []

#Minima Slopes:
for c in range(len(SlopeMatrix)):
    Year = YearList[c]
    print("Starting minima slopes for Year=",Year)
    slopelist = SlopeMatrix[c]


    plt.figure()
    plt.semilogy(SystemSizeList,slopelist)

    plt.xlabel("System Size")
    plt.ylabel("Log of initial yearly R Allele Slope")

    plt.title("Year " + '{:.1E}'.format(Year))
    plt.grid()
    plt.savefig(str(args.directory) +
        "/Minima_Curve_Year_" + '{:.1E}'.format(Year) + ".png")
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

"""
firstplateau = 0
lastplateau = 0

for i in range(len(YearList)):
    if YearList[i] > 0.01:
        firstplateau = i
        break

for i in range(len(YearList)):
    if YearList[i] > 0.98:
        lastplateau = i
        break
    

fitMinimumList = np.asarray(MinimumList[firstplateau:lastplateau])*MAGNTIDUEINCREASE
fitYearList = YearList[firstplateau:lastplateau]*MAGNTIDUEINCREASE

sigmaweights = np.ones(len(fitYearList))
#sigmaweights[:20] = 0.01
#sigmaweights[20:] = 0.01
sigmaweights*=0.01
popt,pcov = curve_fit(fit_func,fitYearList,fitMinimumList,maxfev=100000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
b,c,d = popt


fitYearList = fitYearList/MAGNTIDUEINCREASE
fittedmin1 =  b*( (1 - fitYearList)**c * fitYearList**d) * MAGNTIDUEINCREASE**(c+d-1)

print("DifferenceList: ",fittedmin1 - fitMinimumList/MAGNTIDUEINCREASE)

#Second fit
fitYearList*=MAGNTIDUEINCREASE
popt,pcov = curve_fit(fit_func2,fitYearList,fitMinimumList,maxfev=100000000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
a,b,c,d,e = popt

fitYearList = fitYearList/MAGNTIDUEINCREASE
fittedmin2 = a * MAGNTIDUEINCREASE**(b-1) * (1-fitYearList)**b + c*MAGNTIDUEINCREASE**(d-1) * fitYearList**d + e*MAGNTIDUEINCREASE**-1


"""

result = linregress(np.log(YearList[-20:]), np.log(MinimumList[-20:]))

plt.figure()
plt.loglog(YearList,MinimumList,label="Slope = %0.3f"%(result.slope))
#plt.loglog(fitYearList,fittedmin1,label='fit: product')
#plt.loglog(fitYearList,fittedmin2,label='fit: sum')
plt.legend(loc='upper left')

plt.xlabel("Year")
plt.ylabel("System Size at Minimum")

plt.title("Characteristic: " + '{0:.1E}'.format(k)+ ", OSR Proportion: %0.3f"%(OSR_Proportion) + " PAP: %0.3f"%(PAP))
plt.grid()
plt.savefig(str(args.directory) + "/ComparingMinima.png")
plt.close()

endtime = time.time()

print("Time taken,",endtime-starttime)
