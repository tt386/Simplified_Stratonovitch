import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

#plt.rcParams['text.usetex'] = True

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


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
OSRNumList = []

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
            OSRNumList.append(data['OSRNum'])
            etalist = data["etalist"]
            timetaken = data["timetaken"]
            PAP = data["PAP"]
            Phi = data["Phi"]
            OSRSeparation = data["OSRSeparation"]
            
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
        print("Plotting Sep:",OSRNumList[-1]," Time:",timetaken)
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
                plt.legend(loc='upper right')
                plt.title("Diffusion Coeff: "+ '{0:.5f}'.format(eta))
                
                plt.xlabel("Space")
                plt.grid()
                plt.savefig(i + "/Eta_" + '{0:.5f}'.format(eta) + ".png")
                plt.close()    
        
    except Exception as e:
        print(e)
        print(i," not found")

#Sort these lists
zipped_lists = zip(OSRNumList, SlopeMatrix, LeftMostMinima, RightMostMinima)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
OSRNumList, SlopeMatrix, LeftMostMinima, RightMostMinima= [ list(tuple) for tuple in  tuples]

OSRNumList = np.asarray(OSRNumList)
SlopeMatrix = np.asarray(SlopeMatrix)



#Minima Slopes
for c in range(len(SlopeMatrix)):

    slopelist = np.asarray(SlopeMatrix[c])
    OSRNum = OSRNumList[c]


    #smooth = savgol_filter(slopelist,51,10)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.plot(np.log10(1/np.sqrt(etalist)),np.log10(slopelist),'k',linewidth=5,zorder=0)
    #label="Min L: %0.3f"%((1/np.sqrt(etalist))[slopelist.argmin()]))

    index_min = min(range(len(slopelist)), key=slopelist.__getitem__)

    minL = 1/np.sqrt(etalist)[index_min]
    minR = min(slopelist)

    print("MinL:",minL)
    print("MinR:",minR)

    #plt.scatter([np.log10(minL)],[np.log10(minR)],color='y',s=60,label=r"$L = %0.3f $"%(minL),zorder=1)


    #plt.legend(loc='upper left',fontsize=30)
    """
    #Add markers at minima points
    leftindex = np.where(etalist==LeftMostMinima[c])[0]
    rightindex = np.where(etalist==RightMostMinima[c])[0]
    print("leftindex",leftindex)
    if len(leftindex) > 0:
        plt.loglog(LeftMostMinima[c],slopelist[leftindex],"ko")
        plt.loglog(RightMostMinima[c],slopelist[rightindex],"ko")
    """

    #plt.legend(loc='lower right')

    #plt.xlabel("width of cropping region, " + r"$L$",fontsize=30)
    #plt.ylabel("change in resistant " + r"$\# / L$",fontsize=30)

    plt.xticks([-2,-1,0,1,2])

    ax.set_xticks([-2,-1,0,1,2])
    ax.set_xticklabels(
        ['$10^{-2}$',r'$10^{-1}$',r'$10^{0}$',r'$10^{1}$',r'$10^2$'])

    ax.set_yticks([-4,-3,-2,-1,0])
    ax.set_yticklabels(
        ['$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',r'$10^0$'])

    #plt.grid()

    #plt.title(r"$Num = $ %0.3f"%(OSRNum),fontsize=30)

    plt.xticks(fontsize=30,fontname='Arial')
    plt.yticks(fontsize=30,fontname='Arial')#,rotation=45)


    plt.xticks(fontname = "Arial")
    plt.yticks(fontname = "Arial")

    plt.tight_layout()

    plt.savefig(str(args.directory) + 
        "/Minima_Curve_Sep_" + '{:.1E}'.format(OSRNum) + ".png")
    plt.savefig(str(args.directory) + "/Frame_" + str(c).zfill(5) + ".png")
    plt.savefig(str(args.directory) + "/Abstract_Minimum.png")
    plt.savefig(str(args.directory) + "/Abstract_Minimum.pdf")
    plt.savefig(str(args.directory) + "/Abstract_Minimum.eps")
    plt.close()



    #######################################################################

    ##########################
    #Find maximum of 2nd order
    SmoothedData = savgol_filter(slopelist, 51, 3,mode='mirror')

    smoothed = []
    if PAP < 1:
        smoothed = abs(np.gradient(np.gradient(SmoothedData,1/np.sqrt(etalist)),1/np.sqrt(etalist)))
        #smoothed = savgol_filter(abs(np.gradient(np.gradient(slopelist,1/np.sqrt(etalist)),1/np.sqrt(etalist))), 51, 3)
    else:
        smoothed = slopelist

    maximaindices = argrelextrema(smoothed,np.less)[0]

    maxindex = 0
    
    dist = 100  #Distance from the kink
    if PAP <1:
        for i in reversed(maximaindices):
            if abs((1/np.sqrt(etalist))[i] - 4*np.sqrt(1-PAP)*special.erfinv(1-Phi/(1-Phi))) < dist:
                maxindex = i
                dist = abs((1/np.sqrt(etalist))[i] - 4*np.sqrt(1-PAP)*special.erfinv(1-Phi/(1-Phi)))

    print("MAXINDEX:",maxindex)
    SecondOrderMax_L = ((1/np.sqrt(etalist))[maxindex])
    SecondOrderMax_val = slopelist[maxindex]
    #########################

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.plot(np.log10(1/np.sqrt(etalist)),np.log10(slopelist),'k',linewidth=5,zorder=0)

    plt.plot([np.log10(SecondOrderMax_L),np.log10(SecondOrderMax_L)],[-10,np.log10(SecondOrderMax_val)],'--k',linewidth=5)

    #plt.xticks([-2,-1,0,1,2])

    plt.plot([np.log10(1/np.sqrt(etalist))[-1],np.log10(minL/20)],[-2.5,-np.log10(minL/20)+np.log10(1/np.sqrt(etalist))[-1]-2.5],'g',linewidth='5',zorder=1)

    #Find indices of the heuristic points L=10,20,40

    #10
    # calculate the difference array
    x = 10
    difference_array = np.absolute(1/np.sqrt(etalist)-x)
 
    # find the index of minimum element from the array
    index = difference_array.argmin()

    plt.scatter([np.log10(1/np.sqrt(etalist))[index]],[np.log10(slopelist)[index]],color='Blue',s=750)

    #20
    # calculate the difference array
    x = 20
    difference_array = np.absolute(1/np.sqrt(etalist)-x)

    # find the index of minimum element from the array
    index = difference_array.argmin()

    plt.scatter([np.log10(1/np.sqrt(etalist))[index]],[np.log10(slopelist)[index]],color='red',s=750)

    #40
    # calculate the difference array
    x = 40
    difference_array = np.absolute(1/np.sqrt(etalist)-x)

    # find the index of minimum element from the array
    index = difference_array.argmin()

    plt.scatter([np.log10(1/np.sqrt(etalist))[index]],[np.log10(slopelist)[index]],color='yellow',s=750)





    ax.set_xticks([np.log10(minL),np.log10(SecondOrderMax_L)])
    ax.set_xticklabels(
        [r'$L^{*}$',r'$L_{U}$'])

    
    ax.set_yticks([-4,-2,0])
    ax.set_yticklabels(
        ['$-4$',r'$-2$',r'$0$'])
    
    ax.tick_params(axis='x', which='major', pad=15)

    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    #plt.grid()

    #plt.title(r"$Num = $ %0.3f"%(OSRNum),fontsize=30)

    plt.xticks(fontsize=50,fontname='Arial')
    plt.yticks(fontsize=50,fontname='Arial')#,rotation=45)

    plt.ylim(min(np.log10(slopelist))-0.1,0.1)

    plt.xticks(fontname = "Arial")
    plt.yticks(fontname = "Arial")

    plt.tight_layout()

    plt.savefig(str(args.directory) + "/Abstract_Minimum_UpperBound.eps")


    fig.set_figwidth(20)

    #fig.set_figheight(7.5)

    ax.tick_params(axis='x', which='major', pad=15)

    fig.tight_layout()

    plt.savefig(str(args.directory) +
        "/STRETCHEDMinima_Curve_Sep_" + '{:.1E}'.format(OSRNum) + ".png")
    #plt.savefig(str(args.directory) + "/Frame_" + str(c).zfill(5) + ".png")
    #plt.savefig(str(args.directory) + "/Abstract_Minimum.png")
    #plt.savefig(str(args.directory) + "/Abstract_Minimum.pdf")
    plt.savefig(str(args.directory) + "/STRETCHED_Abstract_Minimum_UpperBound.eps")

    


    #######################################################################


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

    plt.title("OSRNum: %d, OSRSep: %0.3f"%(OSRNum,OSRSeparation))

    plt.savefig(str(args.directory) +
        "/SecondOrder_Sep_" + '{:.1E}'.format(OSRSeparation) + ".png")
    plt.close()

####################################
###Plotting Minima##################
####################################
L1Min  = 1.15

RightMostMinimaTheory = L1Min/(OSRNumList + (OSRNumList-1)*OSRSeparation)

"""
RightMostMinimaTheory = ((RightMostMinima[0]/(OSRNum)**2) 
    * (OSRNum+(OSRNum-1)*OSRSeparationList)**2)
"""

print("Leftmost: ",LeftMostMinima)
plt.figure()

plt.loglog(OSRNumList,1/np.sqrt(LeftMostMinima),label='Rightmost')
plt.loglog(OSRNumList,1/np.sqrt(RightMostMinima),label='Leftmost')

plt.loglog(OSRNumList,RightMostMinimaTheory,"--",label='Leftmost Theory')


plt.loglog(OSRNumList,1.15*np.ones(len(OSRNumList)),'--k',label='Isolated Result')

#w-dependent value
"""
w   l       L
0   100     100
0.1 7.225   6.529
0.2 5.793   4.634
0.3 5.227   3.659
0.4 4.994   2.996
0.5 4.945   2.472
0.6 5.051   2.020
0.7 5.392   1.617
0.8 6.402   1.280
0.9 11.868  1.185
"""
Periodic = 2.472

plt.loglog(OSRNumList,Periodic*np.ones(len(OSRNumList)),'k',label='Periodic Result')


plt.legend(loc='lower left')

plt.grid()

plt.xlabel(r"Num",fontsize=30)
plt.ylabel(r"$\frac{L}{\sqrt{DY}}$ value of Minima",fontsize=30)

plt.title(r"$\delta$: %0.3f"%(OSRSeparation))

plt.tight_layout()

plt.savefig(str(args.directory) + "/Locations_Of_Minima.png")

plt.close()



################################
plt.figure()

plt.plot(OSRNumList,LeftMostMinima,label='Leftmost')
plt.plot(OSRNumList,RightMostMinima,label='Rightmost')

plt.plot(OSRNumList,1/RightMostMinimaTheory**2,"--",label='Rightmost Theory')


plt.plot(OSRNumList,(1/1.15**2)*np.ones(len(OSRNumList)),'--k',label='Isolated Result')

plt.legend(loc='upper right')

plt.grid()

plt.xlabel(r"Num",fontsize=30)
plt.ylabel(r"$\eta$ value of Minima",fontsize=30)

plt.title(r"$\delta$: %0.3f"%(OSRSeparation))

plt.tight_layout()

plt.savefig(str(args.directory) + "/ETALocations_Of_Minima.png")

plt.close()





print("Final value:",(1/np.sqrt(LeftMostMinima))[-1])
