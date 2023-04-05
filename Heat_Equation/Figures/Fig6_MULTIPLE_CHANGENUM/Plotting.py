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
                plt.legend(loc='upper left')
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

    leftL = np.log10(1/np.sqrt(LeftMostMinima[c]))
    rightL = np.log10(1/np.sqrt(RightMostMinima[c]))

    """
    #smooth = savgol_filter(slopelist,51,10)
    plt.figure()
    plt.loglog(1/np.sqrt(etalist),slopelist)
    #label="Min L: %0.3f"%((1/np.sqrt(etalist))[slopelist.argmin()]))


    #plt.legend(loc='lower right')

    plt.xlabel(r"$\frac{L}{\sqrt{DY}}$",fontsize=30)
    plt.ylabel(r"$\Delta R$",fontsize=30)

    plt.grid()

    plt.title(r"$Num = $ %0.3f"%(OSRNum),fontsize=30)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20,rotation=45)

    plt.tight_layout()

    plt.savefig(str(args.directory) + 
        "/Minima_Curve_Sep_" + '{:.1E}'.format(OSRNum) + ".png")
    plt.savefig(str(args.directory) + "/Frame_" + str(c).zfill(5) + ".png")
    plt.close()
    """

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.plot(np.log10(1/np.sqrt(etalist)),np.log10(slopelist),color='k',linewidth=5)

    #plt.plot(np.log10(1/np.sqrt(etalist)),np.log10(slopelist/np.sqrt(etalist)),color='k',linewidth=5)

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


    #plt.grid()

    plt.xticks(fontsize=50,fontname='Arial')
    plt.yticks(fontsize=50,fontname='Arial')

    #plt.title(r"\#=%d"%(OSRNum),fontsize=30)

    """
    plt.xticks([-3,-1,1])

    ax.set_xticks([-3,-1,1])
    ax.set_xticklabels(
        ['$-3$',r'$-1$',r'$1$'])
    """
    ax.set_yticks([-5,-4,-3,-2,-1])
    ax.set_yticklabels(
        [r'$-5$',r'$-4$','$-3$',r'$-2$',r'$-1$'])

    plt.ylim(-5,-0.5)

    plt.tight_layout()

    plt.savefig(str(args.directory) +
        "/Minima_Curve_Sep_" + '{:.1E}'.format(OSRNum) + ".png")
    plt.savefig(str(args.directory) +
        "/Minima_Curve_Num_" + '{:.1E}'.format(OSRNum) + ".eps")
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
fig = plt.figure()
ax=fig.add_subplot(111)

plt.scatter((OSRNumList),np.log10(1/np.sqrt(LeftMostMinima)),label=r'$L_2$',s=400)
plt.scatter((OSRNumList),np.log10(1/np.sqrt(RightMostMinima)),marker='x',label=r'$L_1$',s=400)

plt.plot((OSRNumList),np.log10(RightMostMinimaTheory),"--",color='#d81b60',label=r'$L_1$ Theory',linewidth=5)


plt.plot((OSRNumList),np.log10(1.15*np.ones(len(OSRNumList))),'--k',label=r'$L^*$',linewidth=5)

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

plt.plot((OSRNumList),np.log10(Periodic*np.ones(len(OSRNumList))),'--',color='cyan',label=r'$\left(1-\omega\right)K^*$',linewidth=5)


handles, labels = plt.gca().get_legend_handles_labels()

#specify order of items in legend
order = [4,3,2,0,1]


#plt.legend(loc='center right',fontsize=20)

# Put a legend to the right of the current axis
ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
    loc='center left', bbox_to_anchor=(1.01, 0.5),prop={'size': 50},
    markerscale=1)#5)

#plt.grid()

#plt.xlabel(r"\#",fontsize=30)
#plt.ylabel(r"$\log_{10}\frac{L_{min}}{\sqrt{DY}}$",fontsize=30)

#plt.title(r"$\delta$: %0.3f"%(OSRSeparation))

ax.tick_params(axis='both', which='major', labelsize=50)


ax.set_xticks([0,10,20,30,40])
ax.set_xticklabels(
    [r'$0$',r'$10$',r'$20$',r'$30$',r'$40$'])

ax.set_yticks([-1,0,1])
ax.set_yticklabels(
    [r'$-1$',r'$0$',r'$1$'])


#plt.ylim(-4,-0.5)
plt.xlim(0,40)


plt.xticks(fontsize=50,fontname="Arial")
plt.yticks(fontsize=50,fontname="Arial")

"""
fig_width, fig_height = plt.gcf().get_size_inches()
"""
bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
width, height = bbox.width, bbox.height
print("WIDTH BEFORE:",width)
print("HEIGHT BEFORE:",height)

fig.set_figwidth(20)

"""
fig_width, fig_height = plt.gcf().get_size_inches()

print("WIDTH AFTER:",fig_width)
print("HEIGHT AFTER:",fig_height)
"""
bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
width, height = bbox.width, bbox.height
print("WIDTH Before:",width)
print("HEIGHT After:",height)


plt.tick_params(left = False,labelleft=False)

#plt.tight_layout()

plt.savefig(str(args.directory) + "/Locations_Of_Minima.png")
plt.savefig(str(args.directory) + "/Locations_Of_Minima.eps",
            bbox_inches="tight")

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
