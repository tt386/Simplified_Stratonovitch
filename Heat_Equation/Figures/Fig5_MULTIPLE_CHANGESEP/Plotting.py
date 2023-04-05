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
PAP = 0
OSRSeparationList = []

Phi = 0
dx=0

OSRNum = 0

LeftMostMinima = []
RightMostMinima = []

for i in dirlist:
    try:
        ######################
        ###Downloading Data###
        ######################
        with np.load(os.path.join(i, filename)) as data:
            SlopeMatrix.append(data['RDifferenceList'])
            OSRSeparationList.append(data['OSRSeparation'])
            etalist = data["etalist"]
            timetaken = data["timetaken"]
            PAP = data["PAP"]
            Phi = data["Phi"]
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
        minlist = []

        #SlopeMatrix[-1] = (np.asarray(SlopeMatrix[-1])/np.sqrt(etalist)).tolist()
     
        for j in range(2,len(SlopeMatrix[-1])-2):
            if  (SlopeMatrix[-1][j-1] > SlopeMatrix[-1][j] and 
                SlopeMatrix[-1][j-2] > SlopeMatrix[-1][j] and
                SlopeMatrix[-1][j] < SlopeMatrix[-1][j+1] and
                SlopeMatrix[-1][j] < SlopeMatrix[-1][j+2]):
                minlist.append(etalist[j])
        print(minlist) 

        if len(minlist)==0: 
            minlist.append(0)

        LeftMostMinima.append(minlist[0])
        RightMostMinima.append(minlist[-1])
        
        rlist = np.ones(len(xlist))*(Phi)
        print("Plotting Sep:",OSRSeparationList[-1]," Time:",timetaken)
        Plotting = [0,5,-1]
    
        for p in range(len(etalist)):
            if p%100 == 0:
                eta = etalist[p]
                PAPylist = PAPMatrix[p]
                Yearylist = YearMatrix[p]

                plt.figure()
                plt.loglog(xlist,PAPylist,label='PAP')
                plt.loglog(xlist,Yearylist,label='Year')
                plt.loglog(xlist,rlist/(rlist+Yearylist)-rlist,label='Integrand')
                plt.legend(loc='upper left')
                plt.title("Diffusion Coeff: "+ '{0:.5f}'.format(eta))
                
                plt.xlabel("Space")
                plt.grid()
                plt.savefig(i + "/Eta_" + '{0:.5f}'.format(eta) + ".png")
                plt.close()    
        
    except:
        print(i," not found")

#Sort these lists
zipped_lists = zip(OSRSeparationList, SlopeMatrix, LeftMostMinima, RightMostMinima)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
OSRSeparationList, SlopeMatrix, LeftMostMinima, RightMostMinima= [ list(tuple) for tuple in  tuples]

OSRSeparationList = np.asarray(OSRSeparationList)
SlopeMatrix = np.asarray(SlopeMatrix)



#Minima Slopes
for c in range(len(SlopeMatrix)):

    slopelist = np.asarray(SlopeMatrix[c])
    OSRSep = OSRSeparationList[c]

    """
    leftindex = slopelist.index(LeftMostMinima[c])
    rightindex = slopelist.index(RightMostMinima[c])
    """
    """
    leftindex = np.where(etalist == LeftMostMinima[c])
    rightindex = np.where(etalist == RightMostMinima[c])

    leftL = np.log10(1/np.sqrt(etalist))[leftindex]
    rightL = np.log10(1/np.sqrt(etalist))[rightindex]
    """

    leftL = np.log10(1/np.sqrt(LeftMostMinima[c]))
    rightL = np.log10(1/np.sqrt(RightMostMinima[c]))

    print(OSRSep)
    print("LeftL:",leftL)
    print("RightL:",rightL)

    #smooth = savgol_filter(slopelist,51,10)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(np.log10(1/np.sqrt(etalist)),np.log10(slopelist),color='k',linewidth='5')
    #label="Min L: %0.3f"%((1/np.sqrt(etalist))[slopelist.argmin()]))


    #plt.legend(loc='lower right')

    #plt.xlabel(r"$\log_{10}\frac{L}{\sqrt{DY}}$",fontsize=30)
    #plt.ylabel(r"$\log_{10}\left(\frac{\sqrt{DY}}{L}\Delta R\right)$",fontsize=30)

    ax.tick_params(axis='both', which='major', labelsize=50)


    if rightL != leftL:
        ax.set_xticks([rightL,leftL])
        ax.set_xticklabels(
            [r'$L_1$',r'$L_2$'])

    else:
        ax.set_xticks([rightL])
        ax.set_xticklabels(
            [r'$L_2$'])

    ax.set_yticks([-4,-3,-2])
    ax.set_yticklabels(
        [r'$-4$',r'$-3$',r'$-2$'])


    plt.ylim(-4.5,-1.5)

    plt.xticks(fontsize=50,fontname = "Arial")#rotation=45)
    plt.yticks(fontsize=50,fontname = "Arial")


    #plt.grid()

    #plt.title(r"$\delta = $ %0.3f"%(OSRSep),fontsize=30)


    plt.tight_layout()

    plt.savefig(str(args.directory) + 
        "/Minima_Curve_Sep_" + '{:.1E}'.format(OSRSep) + ".eps")
    plt.savefig(str(args.directory) + "/Frame_" + str(c).zfill(5) + ".png")
    plt.close()


    print("L:",list(1/np.sqrt(etalist)))
    print("S:",list(slopelist))

    ################################
    """
    plt.figure()
    plt.loglog(1/np.sqrt(etalist),slopelist/np.sqrt(etalist))

    plt.xlabel("Size of Isolated Field (L)")
    plt.ylabel("Total # New Resistant")

    plt.grid()

    plt.title("OSRNum: %d, OSRSep: %0.3f"%(OSRNum,OSRSep))

    plt.savefig(str(args.directory) +
        "/Minima_Curve_TOTAL_Sep_" + '{:.1E}'.format(OSRSep) + ".png")
    plt.close()
    """
    ################################
    #Refer to this as 'constrained' because we constain to wanting a certain
    #amount of crop (=TOT)
    TOT = 1e4 #The total amount of crop desired.

    Num = TOT/(1/np.sqrt(etalist)) #Number of widths L that fit into TOT

    """
    plt.figure()
    plt.loglog(1/np.sqrt(etalist),slopelist/np.sqrt(etalist) * Num)

    plt.xlabel("Size of Isolated Field (L)")
    plt.ylabel("Total # New Resistant if we split Massive OSR into n regions size L")

    plt.grid()

    plt.title("OSRNum: %d, OSRSep: %0.3f"%(OSRNum,OSRSep))

    plt.savefig(str(args.directory) +
        "/Minima_Curve_TOTALCONSTRAINED_Sep_" + '{:.1E}'.format(OSRSep) + ".png")
    plt.close()
    """


    #smooth = savgol_filter(slopelist,11,3)

    maxval = max(abs(np.gradient(np.gradient(slopelist,1/np.sqrt(etalist)),1/np.sqrt(etalist))))
    maxx = (1/np.sqrt(etalist))[(abs(np.gradient(np.gradient(slopelist,1/np.sqrt(etalist)),1/np.sqrt(etalist)))).argmax()]

    print("MAX:",maxx,maxval)

    plt.figure()
    plt.loglog(1/np.sqrt(etalist),abs(np.gradient(np.gradient(slopelist,1/np.sqrt(etalist)),1/np.sqrt(etalist))))

    plt.xlabel("Size of Isolated Field (L)")
    plt.ylabel("Second Order Deriv")

    plt.grid()

    plt.title("OSRNum: %d, OSRSep: %0.3f"%(OSRNum,OSRSep))

    plt.savefig(str(args.directory) +
        "/SecondOrder_Sep_" + '{:.1E}'.format(OSRSep) + ".png")
    plt.close()

####################################
###Plotting Minima##################
####################################
L1Min  = 1.15

RightMostMinimaTheory = L1Min/(OSRNum + (OSRNum-1)*OSRSeparationList)

"""
RightMostMinimaTheory = ((RightMostMinima[0]/(OSRNum)**2) 
    * (OSRNum+(OSRNum-1)*OSRSeparationList)**2)
"""

print("Leftmost: ",LeftMostMinima)

fig = plt.figure()
ax = fig.add_subplot(111)



#Get rid of leftmost that occur at same time as rightmost
JUSTRightMostMinima = []
JUSTRightOSRSeparationList = []
for i in range(len(OSRSeparationList)):
    if 1/np.sqrt(LeftMostMinima[i]) != 1/np.sqrt(RightMostMinima[i]):
        JUSTRightMostMinima.append(RightMostMinima[i])
        JUSTRightOSRSeparationList.append(OSRSeparationList[i])


plt.scatter(np.log10(OSRSeparationList),np.log10(1/np.sqrt(LeftMostMinima)),label=r'$L_2$',s=400)

plt.scatter(np.log10(JUSTRightOSRSeparationList),np.log10(1/np.sqrt(JUSTRightMostMinima)),marker='x',label=r'$L_1$',s=400)

plt.plot(np.log10(OSRSeparationList),np.log10(RightMostMinimaTheory),"--",color='#d81b60',label=r'$L_1$ Theory',linewidth=5)

plt.plot(np.log10(OSRSeparationList),np.log10(L1Min*np.ones(len(OSRSeparationList))),'--k',label=r'$L^*$',linewidth=5)

plt.plot([-1,0],[1,0.5],'k',linewidth=10)

#plt.legend(loc='lower left')

handles, labels = plt.gca().get_legend_handles_labels()

#specify order of items in legend
order = [3,2,0,1]


# Put a legend to the right of the current axis
ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': 50},
    markerscale=1)#5)

#lgnd.legendHandles[2]._sizes(150)
#lgnd.legendHandles[3]._sizes(150)

#plt.grid()

#plt.xlabel(r"$\log_{10}\frac{\delta}{\sqrt{DY}}$",fontsize=30)
#plt.ylabel(r"$\log_{10}\frac{L_{min}}{\sqrt{DY}}$",fontsize=30)

#plt.title(r"\#"+" Cropping Regions: %d"%(OSRNum),fontsize=30)

#plt.xticks(fontsize=20,rotation=45)
#plt.yticks(fontsize=20,rotation=45)

ax.tick_params(axis='both', which='major', labelsize=50)


ax.set_xticks([-1,0,1])
ax.set_xticklabels(
    [r'$-1$',r'$0$',r'$1$'])

ax.set_yticks([-1,0,1])
ax.set_yticklabels(
    [r'$-1$',r'$0$',r'$1$'])


plt.ylim(-1,1)
plt.xlim(-1,1)

plt.xticks(fontsize=50,fontname = "Arial")#rotation=45)
plt.yticks(fontsize=50,fontname = "Arial")

fig.set_figwidth(20)

plt.tick_params(left = False,labelleft=False)

#plt.tight_layout()

plt.savefig(str(args.directory) + "/Locations_Of_Minima.png")
plt.savefig(str(args.directory) + "/Locations_Of_Minima.eps",bbox_inches="tight")

plt.close()


slope, intercept, r, p, se = linregress(np.log(OSRSeparationList)[:20],np.log(1/np.sqrt(LeftMostMinima))[:20])

print("SLOPE",slope)
print("stderr",se)

################################
plt.figure()

plt.plot(OSRSeparationList,LeftMostMinima,label='Leftmost')
plt.plot(OSRSeparationList,RightMostMinima,label='Rightmost')

plt.plot(OSRSeparationList,1/RightMostMinimaTheory**2,"--",label='Rightmost Theory')


plt.plot(OSRSeparationList,(1/1.15**2)*np.ones(len(OSRSeparationList)),'--k',label='Isolated Result')

plt.legend(loc='upper right')

plt.grid()

plt.xlabel(r"$\delta$",fontsize=30)
plt.ylabel(r"$\eta$ value of Minima",fontsize=30)

plt.title("Number of Cropping Regions: %d"%(OSRNum))

plt.tight_layout()

plt.savefig(str(args.directory) + "/ETALocations_Of_Minima.png")

plt.close()





print("Final value:",(1/np.sqrt(LeftMostMinima))[-1])
