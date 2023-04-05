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

from scipy.signal import savgol_filter

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
PhiList = []
OSRSep = 0

PAP = 0
dx=0

OSRNum = 0

LeftMostMinima = []
RightMostMinima = []
MinList = []

SecondOrderMaxList = []

for i in dirlist:
    ######################
    ###Downloading Data###
    ######################
    with np.load(os.path.join(i, filename)) as data:
        SlopeMatrix.append(data['RDifferenceList'])
        OSRSep = data['OSRSeparation']
        etalist = data["etalist"]
        timetaken = data["timetaken"]
        PhiList.append(data["Phi"])
        PAP = data["PAP"]
        OSRNum = data["OSRNum"]
        
        xlist = data['xlist']
        PAPMatrix = data["PAPMatrix"]
        YearMatrix = data["YearMatrix"]
        
        timetaken = data['timetaken']
        print("TIME TAKEN:",timetaken)

    ################
    ###Find Minima##
    ################
    """    
    minlistindices = argrelextrema(SlopeMatrix[-1],np.less)
    minlist = etalist[minlistindices[0]]
    """

    """
    minlist = []

    #SlopeMatrix[-1] = (np.asarray(SlopeMatrix[-1])/np.sqrt(etalist)).tolist()
 
    for j in range(1,len(SlopeMatrix[-1])-1):
        if  (SlopeMatrix[-1][j-1] > SlopeMatrix[-1][j] and 
            SlopeMatrix[-1][j] < SlopeMatrix[-1][j+1]):
            minlist.append(etalist[j])
    print(minlist) 

    if len(minlist)==0: 
        minlist.append(0)
    
    LeftMostMinima.append(minlist[0])
    RightMostMinima.append(minlist[-1])
    """
    MinList.append(etalist[SlopeMatrix[-1].argmin()])

    """
    SecondOrderMaxList.append(etalist[abs(np.gradient())])
    """



    
    rlist = np.ones(len(xlist))*(PhiList[-1])
    print("Plotting Phi:",PhiList[-1]," Time:",timetaken)
    Plotting = [0,5,-1]

    #Plotting.append(np.where(etalist==1e-1)[0][0])

    print("Plotting List:",Plotting)

    for p in Plotting:
        eta = etalist[p]
        PAPylist = PAPMatrix[p]
        Yearylist = YearMatrix[p]

        plt.figure()
        plt.semilogy(xlist,PAPylist,label='PAP')
        plt.semilogy(xlist,Yearylist,label='Year')
        plt.semilogy(xlist,rlist/(rlist+Yearylist)-rlist,label='Integrand')
        plt.legend(loc='upper left')
        plt.title("Diffusion Coeff: "+ '{0:.5f}'.format(eta))
        
        plt.xlabel("Space")
        plt.grid()
        plt.savefig(i + "/Eta_" + '{0:.5f}'.format(eta) + ".png")
        plt.close()    
    

#Sort these lists
zipped_lists = zip(PhiList, SlopeMatrix, MinList)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
PhiList, SlopeMatrix, MinList= [ list(tuple) for tuple in  tuples]

PhiList = np.asarray(PhiList)
SlopeMatrix = np.asarray(SlopeMatrix)
MinList = np.asarray(MinList)


#Minima Slopes
for c in range(len(SlopeMatrix)):

    slopelist = np.asarray(SlopeMatrix[c])
    Phi = PhiList[c]

    SmoothedData = savgol_filter(slopelist, 51, 3)

    print("L List:",1/np.sqrt(etalist))
    print("Slopelist:",slopelist)
    #smooth = savgol_filter(slopelist,51,10)
    plt.figure()
    plt.loglog(1/np.sqrt(etalist),slopelist,
        label="Min L: %0.3f"%((1/np.sqrt(etalist))[slopelist.argmin()]))

    plt.loglog(1/np.sqrt(etalist),SmoothedData,'k--',label='Smoothed')

    plt.legend(loc='lower right')

    plt.xlabel("Size of Isolated Field (L)")
    plt.ylabel(r"\# New Resistant per unit Field Size (ie per crop)")

    plt.grid()

    plt.title("Phi: %0.6f, OSRNum: %d, OSRSep: %0.3f"%(Phi,OSRNum,OSRSep))

    plt.savefig(str(args.directory) + 
        "/Minima_Curve_Phi_" + '{:.1E}'.format(Phi) + ".png")
    plt.savefig(str(args.directory) + "/Frame_" + str(c).zfill(5) + ".png")
    plt.close()


    #print("L:",list(1/np.sqrt(etalist)))
    #print("S:",list(slopelist))

    ################################
    plt.figure()
    plt.loglog(1/np.sqrt(etalist),slopelist/np.sqrt(etalist))

    """
    #Add markers at minima points
    leftindex = np.where(etalist==LeftMostMinima[c])[0]
    rightindex = np.where(etalist==RightMostMinima[c])[0]
    print("leftindex",leftindex)
    if len(leftindex) > 0:
        plt.loglog(LeftMostMinima[c],slopelist[leftindex],"ko")
        plt.loglog(RightMostMinima[c],slopelist[rightindex],"ko")
    """
    plt.xlabel("Size of Isolated Field (L)")
    plt.ylabel(r"Total \# New Resistant")

    plt.grid()

    plt.title("Phi: %0.6f, OSRNum: %d, OSRSep: %0.3f"%(Phi,OSRNum,OSRSep))

    plt.savefig(str(args.directory) +
        "/Minima_Curve_TOTAL_Phi_" + '{:.1E}'.format(Phi) + ".png")
    plt.close()

    ################################
    #Refer to this as 'constrained' because we constain to wanting a certain
    #amount of crop (=TOT)
    TOT = 1e4 #The total amount of crop desired.

    Num = TOT/(1/np.sqrt(etalist)) #Number of widths L that fit into TOT

    plt.figure()
    plt.loglog(1/np.sqrt(etalist),slopelist/np.sqrt(etalist) * Num)

    """
    #Add markers at minima points
    leftindex = np.where(etalist==LeftMostMinima[c])[0]
    rightindex = np.where(etalist==RightMostMinima[c])[0]
    print("leftindex",leftindex)
    if len(leftindex) > 0:
        plt.loglog(LeftMostMinima[c],slopelist[leftindex],"ko")
        plt.loglog(RightMostMinima[c],slopelist[rightindex],"ko")
    """
    plt.xlabel("Size of Isolated Field (L)")
    plt.ylabel(r"Total \# New Resistant if we split Massive OSR into n regions size L")

    plt.grid()

    plt.title("Phi: %0.6f, OSRNum: %d, OSRSep: %0.3f"%(Phi,OSRNum,OSRSep))

    plt.savefig(str(args.directory) +
        "/Minima_Curve_TOTALCONSTRAINED_Phi_" + '{:.1E}'.format(Phi) + ".png")
    plt.close()



    #smooth = savgol_filter(slopelist,11,3)

    maxval = max(abs(np.gradient(np.gradient(slopelist,1/np.sqrt(etalist)),1/np.sqrt(etalist))))
    maxx = (1/np.sqrt(etalist))[(abs(np.gradient(np.gradient(slopelist,1/np.sqrt(etalist)),1/np.sqrt(etalist)))).argmax()]

    #SecondOrderMaxList.append(maxx)
    print("Phi:",Phi)
    smoothed = []
    if PAP < 1:
        smoothed = abs(np.gradient(np.gradient(SmoothedData,1/np.sqrt(etalist)),1/np.sqrt(etalist)))
        #smoothed = savgol_filter(abs(np.gradient(np.gradient(slopelist,1/np.sqrt(etalist)),1/np.sqrt(etalist))), 51, 3)
    else:
        smoothed = slopelist
    
    maximaindices = argrelextrema(smoothed,np.less)[0]  #LOOK FOR 2ND MIN RATHER THAN FIRST MAX
    maxindex = maximaindices[0]

    dist = 100
    if PAP <1:
        for i in reversed(maximaindices):
            if abs((1/np.sqrt(etalist))[i] - 4*np.sqrt(1-PAP)*special.erfinv(1-Phi/(1-Phi))) < dist:
                maxindex = i
                dist = abs((1/np.sqrt(etalist))[i] - 4*np.sqrt(1-PAP)*special.erfinv(1-Phi/(1-Phi)))



    """
    for i in reversed(maximaindices):
        if smoothed[i] > 1e-2:
            maxindex = i
            break
    """
    print("MAXINDEX:",maxindex)
    SecondOrderMaxList.append((1/np.sqrt(etalist))[maxindex])

    plt.figure()
    plt.loglog(1/np.sqrt(etalist),abs(np.gradient(np.gradient(slopelist,1/np.sqrt(etalist)),1/np.sqrt(etalist))))
    plt.loglog(1/np.sqrt(etalist),smoothed,label='Smoother')


    plt.legend(loc='upper left')

    plt.xlabel("Size of Isolated Field (L)")
    plt.ylabel("Second Order Deriv")

    plt.grid()

    plt.title("Phi: %0.6f, OSRNum: %d, OSRSep: %0.3f"%(Phi,OSRNum,OSRSep))

    plt.savefig(str(args.directory) +
        "/SecondOrder_Phi_" + '{:.1E}'.format(Phi) + ".png")
    plt.close()

    #########################################################################

    plt.figure()
    plt.loglog(etalist,abs(np.gradient(np.gradient(slopelist,etalist),1/etalist)))

    #plt.legend(loc='upper left')

    plt.xlabel("Size of Isolated Field (L)")
    plt.ylabel("Second Order Deriv")

    plt.grid()

    plt.title("Phi: %0.6f, OSRNum: %d, OSRSep: %0.3f"%(Phi,OSRNum,OSRSep))

    plt.savefig(str(args.directory) +
        "/SecondOrderEta_Phi_" + '{:.1E}'.format(Phi) + ".png")
    plt.close()

####################################
###Plotting Minima##################
####################################
print("Leftmost: ",LeftMostMinima)
plt.figure()

#plt.loglog(PAPList,1/np.sqrt(LeftMostMinima),label='Leftmost')
plt.loglog(PhiList[:-1],1/np.sqrt(MinList)[:-1])#,label='Rightmost')

print("PAPLIST:",np.log10(PhiList))
print("Min L:",np.log10(1/np.sqrt(MinList)))

"""
if PAPList[0] == 0:
    plt.plot(PAPList,RightMostMinimaTheory,"--",label='RightMost Theory')
"""
#plt.legend(loc='upper center')

plt.grid()

plt.xlabel("Phi")
plt.ylabel("L of Minima")

plt.title("Number of OSR: %d"%(OSRNum))

plt.savefig(str(args.directory) + "/LOGLocations_Of_Minima.png")

plt.close()

#############################################################################

slope, intercept, r, p, se = linregress(PhiList[:-1], 1/np.sqrt(MinList)[:-1])
fitxlist = np.linspace(min(PhiList[:-1]),max(PhiList[:-1]),100)

print("Slope:",slope)
print("PhiList:",PhiList)
print("fitxlist:",fitxlist)

plt.figure()

plt.plot(PhiList[:-1],1/np.sqrt(MinList)[:-1])#,label='Rightmost')
plt.plot(fitxlist,slope*fitxlist+intercept,'--k',label='Fit Slope: %0.5f'%(slope))

plt.legend(loc='upper left')

plt.grid()

plt.xlabel("Phi")
plt.ylabel("L of Minima")

plt.title("Number of OSR: %d"%(OSRNum))

plt.savefig(str(args.directory) + "/Locations_Of_Minima.png")

plt.close()


#############################################################################



Theoryx = PhiList[:-1]
Theoryy = 4*special.erfinv(1-Theoryx/(1-Theoryx))*np.sqrt(1-PAP)


plt.figure()

#plt.loglog(PAPList,1/np.sqrt(LeftMostMinima),label='Leftmost')
plt.scatter(np.log10(PhiList[:-1]),np.log10(SecondOrderMaxList[:-1]),color='blue')#,label='Rightmost')
plt.plot(np.log10(Theoryx),np.log10(Theoryy),'--',color='orange',linewidth=3,label='Theory')

plt.grid()

#plt.legend(loc='lower left',fontsize=20)

plt.xlabel(r"$\log_{10}\left(\phi\right)$",fontsize=30)
plt.ylabel(r"$\log_{10}\left(\frac{L}{\sqrt{DY}}\right)$",fontsize=30)

plt.xticks(fontsize=30,rotation=45)
plt.yticks(fontsize=30,rotation=45)


plt.tight_layout()

plt.savefig(str(args.directory) + "/LowerKinkWithPhi.png")

plt.close()

