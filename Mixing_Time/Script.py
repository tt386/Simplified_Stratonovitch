import numpy as np

import numpy as np
import matplotlib.pyplot as plt

import sys

import time

import os

starttime = time.time()
##############################################################################
##############################################################################
##############################################################################
#Trivial System Params
Days = 200

SystemSizes = np.arange(10,400,10)
OSR_Proportions = np.arange(0.1,1.0,0.05)




#System Params    
OSR_Attract = 1.
Refuge_Attract = 0.2

Characteristic = 5



##############################################################################
##############################################################################
##############################################################################
#SaveFileCreation
SaveFileName = ("Saved_Plots/Days_%d_Characteristic_%d"%
                (Days,Characteristic))

if not os.path.isdir(SaveFileName):
    os.mkdir(SaveFileName)


##############################################################################
##############################################################################
##############################################################################

UpperBoundIntersect = []
LowerBoundIntersect = []


for OSR_Proportion in OSR_Proportions:
    print("OSR Proportion:",OSR_Proportion)

    #DataLists
    UpperBoundList = []
    LowerBoundList = []

    UpperBoundList_diff = []

    LowerBoundList_levin= []
    UpperBoundList_levin=[]

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
            Domain, int(RefugeSize/2),'constant',constant_values=Refuge_Attract)

        SystemSize = len(Domain)

        TM = np.zeros((len(Domain),len(Domain)))


        PesticideField = np.ones(OSRWidth)
        PesticideField = np.pad(
            PesticideField,int(RefugeSize/2),'constant',constant_values=0)


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
        EigenVals,EigenVects = np.linalg.eig(TM.T)

        EigenVals = sorted(EigenVals,reverse=True)

        EigenVects = np.transpose(EigenVects)

        for k in range(len(EigenVects[0])):
            EigenVects[0][k] = np.absolute(EigenVects[0][k])

        EigenVect = (EigenVects[0]/sum(EigenVects[0])).real

     
        epsilon = 0.25
        
        t_min = (EigenVals[1]/(1-EigenVals[1])) * np.log(1/(2*epsilon))
        
        t_max = (1/(1-EigenVals[1])) * np.log(1/(min(EigenVects[0])*epsilon))
        
        
        UpperBoundList.append(t_max)
        LowerBoundList.append(t_min)
        
        UpperBoundList_diff.append((1/(1-EigenVals[1]))*np.log(1/(2*epsilon*np.sqrt(min(EigenVects[0])))))
        
        UpperBoundList_levin.append((1/(1-EigenVals[1])) * np.log(4/min(EigenVects[0])))
        LowerBoundList_levin.append((EigenVals[1]/(1-EigenVals[1])) * np.log(2))
        






    for i in range(len(UpperBoundList)):
        if (UpperBoundList[i] < Days) and (UpperBoundList[i+1] > Days):
            UpperBoundIntersect.append(UpperBoundList[i])

        if (LowerBoundList[i] < Days) and (LowerBoundList[i+1] > Days):
            LowerBoundIntersect.append(LowerBoundList[i])





    plt.figure(1)
    plt.plot(SystemSizes,UpperBoundList_diff,label='Different Def Upper')
    plt.plot(SystemSizes,UpperBoundList,label='Upper Bound')
    plt.plot(SystemSizes,LowerBoundList,label='Lower Bound')
    #plt.plot(OSRWidthList,LowerBoundList_levin,label='Levin Lower')
    #plt.plot(OSRWidthList,UpperBoundList_levin,label='Levin Upper')

    plt.plot(SystemSizes,np.ones(len(SystemSizes))*Days,label='Days Per Year')
    plt.legend(loc='upper left')

    plt.title("OSR Proportion: %0.3f"%(OSR_Proportion))
    plt.xlabel("System Width")
    plt.ylabel("Bounds (time steps)")
    plt.grid()
    plt.savefig(SaveFileName + "/MixingTimes_OSRProportion_%0.3f.png"%(OSR_Proportion))

    plt.close()




plt.figure(2)
plt.plot(OSR_Proportions,UpperBoundIntersect,label='Upper Bound')
plt.plot(OSR_Proportions,LowerBoundIntersect,label='Lower Bound')

plt.legend(loc='lower right')

plt.title("Days: %d"%(Days))

plt.ylabel("Intersection System Size")
plt.xlabel("OSR Proportion")
plt.grid()
plt.savefig(SaveFileName + "/MixingTimes_Intersections.png")
plt.close()





endtime = time.time()
print("Time taken:",endtime - starttime)
