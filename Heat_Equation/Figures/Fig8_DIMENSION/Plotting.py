import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

#plt.rcParams['text.usetex'] = True

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

rightmin = []
leftminD = []
leftmin = []

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

            min_indices = []
            for i in range(1,len(IntegralList)-1):
                if  ((IntegralList[i-2] > IntegralList[i-1]) and
                    (IntegralList[i-1] > IntegralList[i]) and
                    (IntegralList[i+1] > IntegralList[i]) and
                    (IntegralList[i+2] > IntegralList[i+1])):
                    min_indices.append(i)

            if len(min_indices) > 2:
                print("TOOOOO MANNNYYYYYY")

            if len(min_indices) > 1:
                rightmin.append(1/np.sqrt(etalist)[min_indices[0]])
                leftminD.append(N)
                leftmin.append(1/np.sqrt(etalist)[min_indices[1]])

            else:
                rightmin.append(1/np.sqrt(etalist)[min_indices[0]])
                leftminD.append(N)
                leftmin.append(0)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.plot(
                np.log10(1/np.sqrt(etalist)),np.log10(IntegralList),
                'k', linewidth=5,
                label='Min: %0.3f'%((1/np.sqrt(etalist))[IntegralList.argmin()]))

            if len(min_indices) > 1:
                plt.scatter([np.log10(1/np.sqrt(etalist)[min_indices[1]])],[np.log10(IntegralList[min_indices[1]])],marker='x',s=80,zorder=5)
                plt.scatter([np.log10(1/np.sqrt(etalist)[min_indices[0]])],[np.log10(IntegralList[min_indices[0]])],marker='o',s=80,zorder=5)
            
            else:
                plt.scatter([np.log10(1/np.sqrt(etalist)[min_indices[0]])],[np.log10(IntegralList[min_indices[0]])],marker='x',s=80,zorder=5)

            #plt.legend(loc='upper right')
            #plt.grid()
            #plt.xlabel("L")
            #plt.ylabel("New Resistant")
            #plt.title("Dimension: %d"%(N))

            ax.tick_params(axis='both', which='major', labelsize=50)

            plt.xticks(fontsize=50,fontname='Arial')
            plt.yticks(fontsize=50,fontname='Arial')


            ax.set_xticks([0,1,2])
            ax.set_xticklabels([r'$0$',r'$1$',r'$2$'])

            """
            ax.set_yticks([0,10,20,30])
            ax.set_yticklabels(
                [r'$0$',r'$10$',r'$20$',r'$30$'])
            """
            plt.tick_params(left = False,labelleft=False)
            plt.tight_layout()
            plt.savefig(str(args.directory) + "/MinimumPlot_Dimension_%d.png"%(N))
            plt.savefig(str(args.directory) + "/MinimumPlot_Dimension_%d.eps"%(N))
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
zipped_lists = zip(DimList, LAtMinList,leftmin,leftminD,rightmin)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
DimList, LAtMinList,leftmin,leftminD,rightmin= [ list(tuple) for tuple in  tuples]

DimList = np.asarray(DimList)
LAtMinList = np.asarray(LAtMinList)


#############################################################################
#############################################################################
#############################################################################

#Curious as to the error on the straight line approx
slope,intercept,r_value,p_value,std_err = linregress(DimList,LAtMinList)

print("For main min")
print("Slope:",slope)
print("intercept:",intercept)
print("rvalue:",r_value)
print("pvalue:",p_value)
print("stderr:",std_err)


slope,intercept,r_value,p_value,std_err = linregress(leftminD[-10:],leftmin[-10:])
print("For Second min")
print("Slope:",slope)
print("intercept:",intercept)
print("rvalue:",r_value)
print("pvalue:",p_value)
print("stderr:",std_err)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(DimList,LAtMinList,linewidth=5,marker='o',markersize=10)

plt.plot(leftminD,leftmin,linewidth=5,marker='o',markersize=10)
plt.plot(DimList,rightmin,linewidth=5,marker='o',markersize=10)

for i in range(len(leftmin)):
    print(leftminD[i],leftmin[i])

plt.plot([1,11],[3,13],'k',linewidth=5)

#ax.set_ylabel('Optimal ' + r'$\frac{L}{\sqrt{DY}}$',fontsize=20)
#ax.set_xlabel('Dimension',fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=50)

plt.xticks(fontsize=50,fontname='Arial')
plt.yticks(fontsize=50,fontname='Arial')


ax.set_xticks([0,10,20])
ax.set_xticklabels([r'$0$',r'$10$',r'$20$'])

ax.set_yticks([0,10,20,30])
ax.set_yticklabels(
    [r'$0$',r'$10$',r'$20$',r'$30$'])

plt.tick_params(left = False,labelleft=False)

plt.xlim(0,25)
plt.ylim(0,35)

plt.tight_layout()
#ax.grid()

plt.savefig(str(args.directory) + "ChangindMinwithDim.png")
plt.savefig(str(args.directory) + "ChangindMinwithDim.eps")


plt.close()



fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(DimList[:10],LAtMinList[:10],linewidth=5,marker='o',markersize=10)

#plt.plot(leftminD,leftmin,linewidth=5,marker='o',markersize=10)
#plt.plot(DimList,rightmin,linewidth=5,marker='o',markersize=10)

#plt.plot([1,11],[3,13],'k',linewidth=5)

#ax.set_ylabel('Optimal ' + r'$\frac{L}{\sqrt{DY}}$',fontsize=20)
#ax.set_xlabel('Dimension',fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=50)

plt.xticks(fontsize=50,fontname='Arial')
plt.yticks(fontsize=50,fontname='Arial')


ax.set_xticks([0,5,10])
ax.set_xticklabels([r'$0$',r'$5$',r'$10$'])

ax.set_yticks([0,10,20,30])
ax.set_yticklabels(
    [r'$0$',r'$10$',r'$20$',r'$30$'])

plt.tick_params(left = False,labelleft=False)

plt.xlim(0,10)
plt.ylim(0,LAtMinList[10])

plt.tight_layout()
#ax.grid()

plt.savefig(str(args.directory) + "10DChangindMinwithDim.png")
plt.savefig(str(args.directory) + "10DChangindMinwithDim.eps")


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

plt.scatter(DimList,LAtMinList/LAtMinList[0],color='blue',s=80)

Rescale = LAtMinList/LAtMinList[0]

slope,intercept,r_value,p_value,std_err = linregress(DimList,Rescale)

plt.plot(DimList,slope*DimList+intercept,'--',color='orange',linewidth=3)

ax.set_ylabel('Optimal ' + r'$\frac{\tilde{L}}{\sqrt{DY}}$' + ' (/1D Optimum)',fontsize=20)
ax.set_xlabel('Dimension ' + r'$D$',fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=30)


plt.xticks(fontsize=30,rotation=45)
plt.yticks(fontsize=30,rotation=45)

plt.grid()

plt.tight_layout()

plt.savefig(str(args.directory) + "ChangindMinwithDim_RescaleFirst.png")

plt.close()



