import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
import time

from scipy.stats import linregress

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
    CharacteristicList = data['CharacteristicList']
    SystemSizeList = data["SystemSizeList"]
    timetaken = data["timetaken"]
    PAP = data["PAP"]
#SlopeMatrix = np.zeros((len(CharacteristicList),len(SystemSizeList)))

MinimumList = []

#Minima Slopes:
for c in range(len(SlopeMatrix)):
    k = CharacteristicList[c]
    print("Starting minima slopes for k=",k)
    slopelist = SlopeMatrix[c]


    plt.figure()
    plt.semilogy(SystemSizeList,slopelist)

    plt.xlabel("System Size")
    plt.ylabel("Log of initial yearly R Allele Slope")

    plt.title("Characteristic " + '{:.1E}'.format(k))
    plt.grid()
    plt.savefig(str(args.directory) +
        "/Minima_Curve_Characteristic_" + '{:.1E}'.format(k) + ".png")
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



result = linregress(np.log(CharacteristicList), np.log(MinimumList))

plt.figure()
plt.loglog(CharacteristicList,MinimumList,label="Slope = %0.3f"%(result.slope))

plt.legend(loc='upper left')

plt.xlabel("Characteristic Migration Scale")
plt.ylabel("System Size at Minimum")

plt.title("PAP: %0.5f"%(PAP))
plt.grid()
plt.savefig(str(args.directory) + "/ComparingMinima.png")
plt.close()

endtime = time.time()

print("Time taken,",endtime-starttime)
