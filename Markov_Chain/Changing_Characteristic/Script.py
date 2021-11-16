import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import copy
import os
from scipy.signal import savgol_filter
from scipy.signal import argrelextrema
from fractions import Fraction

starttime = time.time()


##############################################################################
##############################################################################
##############################################################################
#Functions

def LinearPesticide(XYPop,PesticideField,Pesticide_XY,Pesticide_Application_Period,d):
    return (XYPop[PesticideField==1] * (1. - Pesticide_XY +
            (Pesticide_XY/Pesticide_Application_Period) *d))


def FlatPesticide(XYPop,PesticideField,Pesticide_XY):
    return (XYPop[PesticideField==1] * (1.-Pesticide_XY))


##############################################################################
##############################################################################
##############################################################################
#Trivial Variable Params
Years = 5
Days = 200

R_Allele_Ratio = 1e-5

#SystemSizes = np.arange(10,600,10)

OSR_Proportion = 0.5

Characteristics = np.linspace(0.5,5,10)

#System Params
FLATPESTICIDE = True

OSR_Attract = 1.
Refuge_Attract = 1.#0.2         #############################################CHANGESCHANGED####################################################

Pesticide_Application_Period = 10
Pesticide_SS = 1#0.86
Pesticide_RR = 0.003
dR = 48.6
Pesticide_SR = (Pesticide_SS + dR*Pesticide_RR)/(1.+dR)



BreedingTime = Days-1 #int(Days/2)



GradientRange = 2   #Smallest value is 2


##############################################################################
##############################################################################
##############################################################################
#SaveFileCreation
SaveFileName = ("Saved_Plots/Years_%d_PAP_%d_gSS_%0.3f_OSRProportion_%0.3f_RefugeAttract_%0.3f_GradientRange_%d_BreedingTime_%d"%
                (Years,Pesticide_Application_Period,Pesticide_SS,OSR_Proportion,Refuge_Attract,GradientRange,BreedingTime))
if FLATPESTICIDE:
    SaveFileName += "_FlatPesticide"
else:
    SaveFileName += "_LinearPesticide"

if not os.path.isdir(SaveFileName):
    os.mkdir(SaveFileName)

##############################################################################
##############################################################################
##############################################################################



#Derived Params

MinPointyList = []
MinPointxList = []

MinCharacteristics = []

SmoothedMinPointyList = []
SmoothedMinPointxList = []

SmoothedCharacteristics = []

for Characteristic in Characteristics:
    InitialSlopeList = []

    FractionRepresentation = Fraction(OSR_Proportion).limit_denominator()
    print(FractionRepresentation)
    Denominator = FractionRepresentation.denominator
    print(Denominator)

    while Denominator < 5:
        Denominator *=2

    SystemSizes = np.arange(Denominator,400,Denominator)


    for SystemSize in SystemSizes:

        print("SystemSize:", SystemSize)

        OSRWidth = int(SystemSize * OSR_Proportion)
        RefugeSize = int(SystemSize * (1-OSR_Proportion))

        SystemSize = OSRWidth + RefugeSize


        ##############################################################################
        ##############################################################################
        ##############################################################################
        #Transition Matrix
        Domain = np.ones(OSRWidth)*OSR_Attract
        Domain = np.pad(
            Domain,(0,RefugeSize),'constant',constant_values=Refuge_Attract)
        """
        Domain = np.pad(
            Domain, int(RefugeSize/2),'constant',constant_values=Refuge_Attract)
        """
        SystemSize = len(Domain)

        TM = np.zeros((len(Domain),len(Domain)))


        PesticideField = np.ones(OSRWidth)
        PesticideField = np.pad(
            PesticideField,(0,RefugeSize),'constant',constant_values=0)

        """
        PesticideField = np.pad(
            PesticideField,int(RefugeSize/2),'constant',constant_values=0)
        """

        #MigrationList
        #List that, if migration occurs, corresponds to migrating to a spot
        MigrateList = np.zeros((SystemSize,SystemSize))
        for x in range(len(TM)):
            for j in range(len(Domain)):
                dist = abs(j-x)
                if dist > SystemSize/2:
                    dist = SystemSize-dist
                
                force = Domain[j] * np.exp(-(dist**2)/(2*(Characteristic**2)))
                MigrateList[x][j] = force
            
        row_sums = MigrateList.sum(axis=1)
        TM = MigrateList / row_sums[:, np.newaxis]


        ##############################################################################
        ##############################################################################
        ##############################################################################
        #Initialise Population
        #Eigenvector calc
        Vals,Vects = np.linalg.eig(TM.T)

        Vals = Vals.real
        Vects=Vects.real

        idx = Vals.argsort()[::-1]
        Vals = Vals[idx]
        Vects = Vects[:,idx]

        Vects = np.transpose(Vects)


        for k in range(len(Vects[0])):
            Vects[0][k] = np.absolute(Vects[0][k])

        EigenVect = (Vects[0]/sum(Vects[0])).real

        SSPop = EigenVect * (1.-R_Allele_Ratio)**2
        SRPop = EigenVect * 2*R_Allele_Ratio*(1.-R_Allele_Ratio)
        RRPop = EigenVect * R_Allele_Ratio**2


        ##############################################################################
        ##############################################################################
        ##############################################################################
        #Main Process

        #Plotting Lists
        SSPop_Prop_Daily= []
        SRPop_Prop_Daily = []
        RRPop_Prop_Daily = []

        SSPop_Daily = []
        SRPop_Daily = []
        RRPop_Daily = []

        R_Ratio_Prop_Daily = []

        SSPop_Yearly = []


        R_Ratio_Yearly = []


        SSPop_Yearly.append(SSPop)
        R_Ratio_Yearly.append(np.sum(2*RRPop + SRPop)/np.sum(2*(
            SSPop + SRPop + RRPop)))


        #InitialValues:
        Tot = SSPop + SRPop + RRPop
        SSPop_Prop_Daily.append(SSPop/Tot)
        SRPop_Prop_Daily.append(SRPop/Tot)
        RRPop_Prop_Daily.append(RRPop/Tot)

        R_Ratio_Prop_Daily.append((2*RRPop + SRPop)/(2*Tot))

        SSPop_Daily.append(SSPop)
        SRPop_Daily.append(SRPop)
        RRPop_Daily.append(RRPop)



        for y in range(Years):
            
            for d in range(Days):
                #print("SS in middle:",SSPop[int(len(SSPop)/2)])
                ######################################################################
                ######################################################################
                ######################################################################
                #Death
                if d < Pesticide_Application_Period:
                    beforesssum = np.sum(2*RRPop + SRPop)/(2*(np.sum(SSPop+SRPop+RRPop)))



                    BeforeSS = np.sum(SSPop)
                    BeforeSR = np.sum(SRPop)
                    BeforeRR = np.sum(RRPop)
                   
                    #print("Before SS:",BeforeSS) 
                    #print("Before SR:",BeforeSR)
                    #print("Before RR:",BeforeRR)


         
                    if not FLATPESTICIDE:
                        SSPop[PesticideField==1] = LinearPesticide(SSPop,PesticideField,Pesticide_SS,Pesticide_Application_Period,d)
                        SRPop[PesticideField==1] = LinearPesticide(SRPop,PesticideField,Pesticide_SR,Pesticide_Application_Period,d)
                        RRPop[PesticideField==1] = LinearPesticide(RRPop,PesticideField,Pesticide_RR,Pesticide_Application_Period,d)
                    
                    else:
                        SSPop[PesticideField==1] = FlatPesticide(SSPop,PesticideField,Pesticide_SS)
                        SRPop[PesticideField==1] = FlatPesticide(SRPop,PesticideField,Pesticide_SR)
                        RRPop[PesticideField==1] = FlatPesticide(RRPop,PesticideField,Pesticide_RR)
                    

                    aftersssum = np.sum(2*RRPop + SRPop)/(2*(np.sum(SSPop+SRPop+RRPop)))

                                
                    AfterSS = np.sum(SSPop)
                    AfterSR = np.sum(SRPop)
                    AfterRR = np.sum(RRPop)

                    #print("After SS:",AfterSS)
                    #print("After SR:",AfterSR)
                    #print("After RR:",AfterRR)

                    #print("SS Ratio:",AfterSS/BeforeSS)
                    #print("SR Ratio:",AfterSR/BeforeSR)
                    #print("RR Ratio:",AfterRR/BeforeRR)


                    #print(aftersssum/beforesssum)
                    
                ######################################################################
                ######################################################################
                ######################################################################
                #Breeding
                if d == BreedingTime:

                    #print("BREED")

                    B_SSPop = SSPop**2 + 2 * 0.5 * SSPop*SRPop + 0.25*SRPop**2
                    
                    B_SRPop = (2 * 0.5 * SSPop*SRPop + 2 * SSPop * RRPop + 
                               0.5 * SRPop**2 + 2 * 0.5 * RRPop*SRPop)

                    B_RRPop = RRPop**2 + 2 * 0.5 * RRPop*SRPop + 0.25*SRPop**2
                
                
                    R_Before = (2*RRPop+SRPop)/(2*(RRPop+SRPop+SSPop))
                    R_After = (2*B_RRPop+B_SRPop)/(2*(B_RRPop+B_SRPop+B_SSPop))

                    #print("R Ratio:",R_After/R_Before)            

                    total = B_SSPop + B_SRPop + B_RRPop
         
                    SSPop = B_SSPop/total
                    SRPop = B_SRPop/total
                    RRPop = B_RRPop/total
                    
                
                ######################################################################
                ######################################################################
                ######################################################################
                #Migration
                Tot = SSPop + SRPop + RRPop
                #print("Before Migration, R Ratio:\t",np.sum(2*RRPop + SRPop)/np.sum(2*Tot))
                SSPop = np.matmul(SSPop,TM)
                SRPop = np.matmul(SRPop,TM)
                RRPop = np.matmul(RRPop,TM)


               


                Tot = SSPop + SRPop + RRPop
                SSPop_Prop_Daily.append(SSPop/Tot)
                SRPop_Prop_Daily.append(SRPop/Tot)
                RRPop_Prop_Daily.append(RRPop/Tot)
                
                SSPop_Daily.append(SSPop)
                SRPop_Daily.append(SRPop)
                RRPop_Daily.append(RRPop)

                R_Ratio_Prop_Daily.append((2*RRPop + SRPop)/(2*Tot))

                #print("After Migration, R Ratio:\t",np.sum(2*RRPop + SRPop)/np.sum(2*Tot))

            SSPop_Yearly.append(SSPop)
            R_Ratio_Yearly.append(np.sum(2*RRPop + SRPop)/np.sum(2*(
                    SSPop + SRPop + RRPop)))
        ##############################################################################
        ##############################################################################
        ##############################################################################
        #Plotting for yearly SSPop
        """
        xlist = np.arange(len(SSPop))
        plt.figure(95)
        for i in range(len(SSPop_Yearly)):
            plt.plot(xlist,SSPop_Yearly[i],label='%d'%(i))
        plt.legend(loc='upper right')
        plt.grid()
        #plt.ylim(0,max(EndState))
        #plt.show()   

        print(SSPop_Yearly[-1])
        print(sum(SSPop_Yearly[-1]))
        """
        ##############################################################################
        ##############################################################################
        ##############################################################################

        ##############################################################################
        ##############################################################################
        ##############################################################################
        #R allele Yearly

        #Initial slope
        x = np.arange(len(R_Ratio_Yearly))
        y = np.log10(R_Ratio_Yearly)

        x_init = x[0:GradientRange]#x[1:4]
        y_init = y[0:GradientRange]#y[1:4]

        print("x_init, y_init: ",x_init,y_init)

        x_final = x[-3:]
        y_final = y[-3:]

        print("x_final, y_final: ",x_final,y_final)

        m_init, b_init = np.polyfit(x_init, y_init, 1)
        m_final, b_final = np.polyfit(x_final, y_final, 1)


        InitialSlopeList.append(m_init)





    #########################################################################
    #########################################################################
    #########################################################################

    #SmoothedInitialSlopeList =  savgol_filter(InitialSlopeList, 5, 3) # window size 51, polynomial order 3
    windowwidth = int(0.25*len(InitialSlopeList))
    if (windowwidth%2 == 0):
            windowwidth +=1
    print("Window width",windowwidth)
    SmoothedInitialSlopeList =  savgol_filter(InitialSlopeList, windowwidth, 3) # window size half system, polynomial order 3


    minpoint = 0
    for i in range(1,len(InitialSlopeList)-1):
        if (InitialSlopeList[i+1] > InitialSlopeList[i]) and (InitialSlopeList[i-1] > InitialSlopeList[i]):
            MinPointyList.append(InitialSlopeList[i])
            MinPointxList.append(SystemSizes[i])
            MinCharacteristics.append(Characteristic)
            break


    minpointindex = argrelextrema(SmoothedInitialSlopeList, np.less)[0]
    if len(minpointindex) > 0:
        SmoothedMinPointyList.append(SmoothedInitialSlopeList[minpointindex][-1])
        SmoothedMinPointxList.append(SystemSizes[minpointindex][-1])

        SmoothedCharacteristics.append(Characteristic)


    #########################################################################
    #########################################################################
    #########################################################################


    plt.figure(1)
    plt.semilogy(
        SystemSizes,
        InitialSlopeList,
        label='Data')
    plt.semilogy(
        SystemSizes,
        SmoothedInitialSlopeList,
        label='Smoothed')


    for i in minpointindex:
        plt.semilogy(SystemSizes[i],SmoothedInitialSlopeList[i],'o',markersize=5,markeredgecolor="black", markerfacecolor="red")

    plt.legend(loc='lower right')
    plt.grid()
    plt.xlabel("System Width")
    plt.ylabel("Initial R Allele Slope")

    plt.title("Characteristic: %0.3f"%(Characteristic))
    plt.savefig(SaveFileName + "/InitialSlopePlot_Characteristic_%0.3f.png"%(Characteristic)) 
    plt.close()


    #########################################################################
    #########################################################################
    #########################################################################



m, b = np.polyfit(SmoothedCharacteristics, SmoothedMinPointxList, 1)

SmoothedCharacteristics = np.asarray(SmoothedCharacteristics)

plt.figure(2)
plt.plot(
    MinCharacteristics,
    MinPointxList,
    label="Data")
plt.plot(
    SmoothedCharacteristics,
    SmoothedMinPointxList,
    label="Smoothed")

plt.plot(
    SmoothedCharacteristics,
    m * SmoothedCharacteristics + b,
    label='Fitted, Slope %0.3f, Int %0.3f'%(m,b))

plt.legend(loc='lower right')

plt.title("OSR Proportion: %0.3f"%(OSR_Proportion))

plt.grid()
plt.xlabel("Characteristic Migration Length")
plt.ylabel("System size of minimum point")

#plt.title("OSR Proportion: %0.3f"%(OSR_Proportion))
plt.savefig(SaveFileName + "/MinInitialSlopePlot_SystemSize.png")
plt.close()




plt.figure(3)
plt.plot(
    MinCharacteristics,
    MinPointyList,
    label="Data")

plt.plot(
    SmoothedCharacteristics,
    SmoothedMinPointyList,
    label="Smoothed")

plt.legend(loc='lower right')

plt.grid()
plt.xlabel("Characteristic Migration Length")
plt.ylabel("Min Point Slope Magnitude")

#plt.title("OSR Proportion: %0.3f"%(OSR_Proportion))
plt.savefig(SaveFileName + "/MinInitialSlopePlot_Magnitude.png")
#plt.show()
plt.close()




endtime = time.time()
print("Time taken:",endtime - starttime)

