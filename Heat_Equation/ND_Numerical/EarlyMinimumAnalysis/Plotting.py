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

import os

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


#Find list of all the datafiles
tempdirlist = os.listdir(args.directory)
dirlist = []
for i in tempdirlist:
    if os.path.isdir(os.path.join(args.directory,i)):
        dirlist.append(os.path.join(args.directory,i))

print("Dirlist:",dirlist)


DimList = []
LAtMinList = []

for i in dirlist:
    try:
        with np.load(os.path.join(i,filename)) as data:
            N = data["N"]
            print("Dimension Being Plotted:",N)
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
                1/np.sqrt(etalist),IntegralList,
                label='Min: %0.3f'%((1/np.sqrt(etalist))[IntegralList.argmin()]))
            plt.legend(loc='upper right')
            plt.grid()
            plt.xlabel("L")
            plt.ylabel("New Resistant")
            plt.title("Dimension: %d"%(N))
            plt.savefig(str(args.directory) + "/MinimumPlot_Dimension_%d.png"%(N))
            plt.close()

            #print("L:",list(1/np.sqrt(etalist)))
            #print("Slopelist:",list(IntegralList))

            #########################################################################
            """
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
            """
            #########################################################################

            #Locate the Minimum:
            DimList.append(N)
            LAtMinList.append((1/np.sqrt(etalist))[IntegralList.argmin()])
    except:
        print(i,"Unavailable")


#Sort these lists
zipped_lists = zip(DimList, LAtMinList)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
DimList, LAtMinList= [ list(tuple) for tuple in  tuples]

DimList = np.asarray(DimList)
LAtMinList = np.asarray(LAtMinList)


#############################################################################
#############################################################################
#############################################################################

#Curious as to the error on the straight line approx
slope,intercept,r_value,p_value,std_err = linregress(DimList,LAtMinList)

print("Slope:",slope)
print("intercept:",intercept)
print("rvalue:",r_value)
print("pvalue:",p_value)
print("stderr:",std_err)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(DimList,LAtMinList,marker='o')


ax.set_ylabel('Optimal ' + r'$\frac{L}{\sqrt{DY}}$',fontsize=20)
ax.set_xlabel('Dimension',fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=20)

plt.xticks(fontsize=20,rotation=45)
plt.yticks(fontsize=20,rotation=45)


plt.tight_layout()
ax.grid()

plt.savefig(str(args.directory) + "ChangindMinwithDim.png")

plt.close()
#############################################################################
#############################################################################
#############################################################################

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(DimList,LAtMinList/DimList,marker='o')


ax.set_ylabel('Optimal ' + r'$\frac{L}{\sqrt{DY}} (/D)$',fontsize=20)
ax.set_xlabel('Dimension ' + r"$D$",fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=20)

plt.xticks(fontsize=20,rotation=45)
plt.yticks(fontsize=20,rotation=45)


plt.tight_layout()
plt.grid()


ax.grid()

plt.savefig(str(args.directory) + "ChangindMinwithDim_RescaleDim.png")

plt.close()

#############################################################################
#############################################################################
#############################################################################

fig = plt.figure()
ax = fig.add_subplot(111)

plt.scatter(DimList,LAtMinList/LAtMinList[0],color='blue')

Rescale = LAtMinList/LAtMinList[0]

slope,intercept,r_value,p_value,std_err = linregress(DimList,Rescale)

plt.plot(DimList,slope*DimList+intercept,'--',color='orange')

ax.set_ylabel('Optimal ' + r'$\frac{L}{\sqrt{DY}}$' + ' (/1D Optimum)',fontsize=20)
ax.set_xlabel('Dimension ' + r'$D$',fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=20)

plt.xticks(fontsize=20,rotation=45)
plt.yticks(fontsize=20,rotation=45)


plt.tight_layout()
plt.grid()


ax.grid()

plt.savefig(str(args.directory) + "ChangindMinwithDim_RescaleFirst.png")

plt.close()



