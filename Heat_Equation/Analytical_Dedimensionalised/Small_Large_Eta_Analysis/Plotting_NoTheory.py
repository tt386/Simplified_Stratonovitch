import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
import time

from scipy.stats import linregress

from scipy.optimize import curve_fit

import sys

import scipy.integrate as integrate
import scipy.special as special

from scipy.signal import argrelextrema

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


################
#Create new saved directory 
NoTheoryDir = str(args.directory) + "/NoTheory"
if not os.path.isdir(NoTheoryDir):
    os.mkdir(NoTheoryDir)
    print("Created Directory")

################



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
    if i != "NoTheory":
        print(i)
        if os.path.isdir(os.path.join(args.directory,i)):
            dirlist.append(os.path.join(args.directory,i))



HeatMap_rhoList = []
HeatMap_etaList = []
HeatMap_slopeList = []


etaList = 0
timetaken = 0
PAP = 0
rhoList = []

Phi = 0
dx=0
for i in dirlist:
    with np.load(os.path.join(i, filename)) as data:
        SlopeMatrix.append(data['SlopeList'])
        minmatrix.append(data['minlist'])
        rhoList.append(float(data['rho']))
        etaList = data["etaList"]
        timetaken = data["timetaken"]
        PAP = data["PAP"]
        Phi = data["Phi"] 
        dx = data["dx"]

    if rhoList[-1] != 0:
        for j in range(len(etaList)):
            HeatMap_rhoList.append(rhoList[-1])
            HeatMap_etaList.append(etaList[j])
            HeatMap_slopeList.append(SlopeMatrix[-1][j])


#SlopeMatrix = np.zeros((len(CharacteristicList),len(SystemSizeList)))

#Correctly order SlopeMatrxi and OSRProportion
"""
rhoList,SlopeMatrix,minmatrix = zip(*sorted(zip(rhoList,SlopeMatrix,minmatrix)))
"""
zipped_lists = zip(rhoList, SlopeMatrix,minmatrix)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
rhoList, SlopeMatrix, minmatrix= [ list(tuple) for tuple in  tuples]

rhoList = np.asarray(rhoList)
SlopeMatrix = np.asarray(SlopeMatrix)
minmatrix = np.asarray(minmatrix)
print(rhoList)

#sys.exit()


MinimumList = []

max_minpointlist = []

SecondOrderDerivList = []
SecondOrderDerivList_rho = []

SecondOrderDerivList_min = []
SecondOrderDerivList_rho_min = []

#Minima Slopes:
for c in range(len(SlopeMatrix)):
    Refuge_Proportion = rhoList[c]
    print("Starting minima slopes for rho=",Refuge_Proportion)
    slopelist = SlopeMatrix[c]


    rho = Refuge_Proportion
    b = np.sqrt(4*etaList*(1-PAP))
    Phiprime = 1-Phi

    #########################################################################
    #########################################################################
    #########################################################################
    def SmallEtaIntegrand(b,rho,xlist,Phi):
        from_rho = 0.5*(1+special.erf((rho-xlist)/b))
        from_1 = 0.5*(1-special.erf((1-xlist)/b))
        from_0 = 0.5*(special.erf((xlist)/b) - 1)
        from_1_Plus_rho = 0.5*(-1+special.erf((rho+1-xlist)/b))
        c = from_rho + from_1 + from_0 + from_1_Plus_rho

        return Phi/(Phi + c)
    #########################################################################
    #########################################################################
    #########################################################################
    def LargeEtaIntegrand(b,rho,xlist,Phi):
        if rho > 0:
            ACon = []
            for x in xlist:
                z = (x-rho)/b
                first = np.imag(np.exp(1j*np.pi*x/rho) *special.erf(z+1j*b*np.pi/(2*rho)) )
            
                z = x/b
                second = np.imag(np.exp(1j*np.pi*x/rho) *special.erf(z+1j*b*np.pi/(2*rho)) )
            
                ACon.append(  2/np.pi * np.exp(-eta*(np.pi/rho)**2) * (second-first))

            #We shall implememt a poor man's periodicity, by having the curve on the right too
            ACon_right = []
            for x in xlist:
                X = x-1 #Translate curve to the right
                
                z = (X-rho)/b
                first = np.imag(np.exp(1j*np.pi*X/rho) *special.erf(z+1j*b*np.pi/(2*rho))) 
            
                z = X/b
                second = np.imag(np.exp(1j*np.pi*X/rho) *special.erf(z+1j*b*np.pi/(2*rho)) )
            
                ACon_right.append(  2/np.pi * np.exp(-eta*(np.pi/rho)**2) * (second-first))

            ACon = np.asarray(ACon)
            ACon_right = np.asarray(ACon_right)
            return Phi/(Phi+ACon+ACon_right)

        else:
            return Phi*np.ones(len(xlist))



    #Change from 'Slope' to endstate number
    slopelist = np.exp(slopelist + np.log(Phi))

    plt.figure()
    plt.loglog(etaList,slopelist,label="Data")

    plt.legend(loc='lower left')

    plt.xlabel("eta")
    plt.ylabel("Log of initial yearly R Allele Slope")

    plt.title("Refuge Proportion " + '{:.1E}'.format(Refuge_Proportion))
    plt.grid()
    plt.savefig(NoTheoryDir +
        "/Minima_Curve_rho_" + '{:.1E}'.format(Refuge_Proportion) + ".png")
    plt.close()
    

    #Plot of slopes of the curves
    plt.figure()
    plt.loglog(etaList,abs(np.gradient(slopelist)),label="Data")

    plt.legend(loc='lower left')

    plt.xlabel("eta")
    plt.ylabel("Slope")

    plt.title("Refuge Proportion " + '{:.1E}'.format(Refuge_Proportion))

    plt.grid()

    plt.savefig(NoTheoryDir +
        "/Gradient_Curve_rho_" + '{:.1E}'.format(Refuge_Proportion) + ".png")
    plt.close()


    #Second order derivative
    plt.figure()
    plt.loglog(etaList,abs(np.gradient(np.gradient(slopelist))),label="Data")

    plt.legend(loc='lower left')

    plt.title("Refuge Proportion " + '{:.1E}'.format(Refuge_Proportion))

    plt.xlabel("eta")
    plt.ylabel("Second Derivative")

    plt.grid()

    plt.savefig(NoTheoryDir +
        "/SecondOrderDeriv_Curve_rho_" + '{:.1E}'.format(Refuge_Proportion) + ".png")
    plt.close()


    #finding max points
    if rho > 0:
        SecondOrderDerivMaxList = [argrelextrema(abs(np.gradient(np.gradient(slopelist))), np.greater)]
        #print("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
        #print(SecondOrderDerivMaxList)
        #print(SecondOrderDerivMaxList[0])
        #print(SecondOrderDerivMaxList[0][0])
        #print("BBBBBBBBBBBBBBBBBBBBBBBBBBB")
        try:
            #print(SecondOrderDerivMaxList)
            #print(SecondOrderDerivMaxList[0][0][1])
            SecondOrderDerivList.append(etaList[SecondOrderDerivMaxList[0][0][1]])
            SecondOrderDerivList_rho.append(rho)
        except:
            print("no maxima found")

        #print("SecondOrderDerivList:",SecondOrderDerivList)


    #finding min points
    if rho > 0:
        SecondOrderDerivMinList = [argrelextrema(abs(np.gradient(np.gradient(slopelist))), np.less)]
        #print("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
        #print("Second order min:",SecondOrderDerivMinList)
        #print("Second order min:",SecondOrderDerivMinList[0])
        #print("Second order min:",SecondOrderDerivMinList[0][0])
        #print("BBBBBBBBBBBBBBBBBBBBBBBBBBB")
        try:
            #print(SecondOrderDerivMinList)
            #print(SecondOrderDerivMinList[0][0][0])
            SecondOrderDerivList_min.append(etaList[SecondOrderDerivMinList[0][0][0]])
            SecondOrderDerivList_rho_min.append(rho)
        except:
            print("no minima found")

        #print("SecondOrderDerivList:",SecondOrderDerivList)


    #Find the minumum of the curve
    minRAllele = 10000000
    minSystemSize = -1
    for i in range(len(slopelist)):
        if slopelist[i] < minRAllele:
            minRAllele = slopelist[i]
            minSystemSize = etaList[i]

    MinimumList.append(minSystemSize)
    print("minRAllele:",minRAllele)
    print("minSystemSize",minSystemSize)


############################################################################
############################################################################
############################################################################
print("LOWER KINK Rho List",SecondOrderDerivList_rho_min)
print("LOWER KINK Eta List",SecondOrderDerivList_min)


print("UPPER KINK Rho List",SecondOrderDerivList_rho)
print("UPPER KINK Eta List",SecondOrderDerivList)

############################################################################
############################################################################
############################################################################
x = 1.-np.asarray(SecondOrderDerivList_rho)
y = np.asarray(SecondOrderDerivList)

x_slope = np.log(x[x>0.3])
y_slope = np.log(y[x>0.3])

slope, intercept, r, p, se = linregress(x_slope, y_slope)


theorykinks = (x)**2/(16*(1-PAP) * special.erfinv(1-Phi)**2)

plt.figure()
plt.loglog(x,y,label="Data")
plt.loglog(x,theorykinks,"--",label='Theory')
plt.legend(loc='upper left')
plt.grid()
plt.xlabel("1-rho")
plt.ylabel("eta at kink")
plt.title("How initial kink changes with rho")

plt.xlim(0.3,1)
plt.savefig(NoTheoryDir + "/LowerKink_ChangingRho.png")
plt.close()

############################################################################
############################################################################
############################################################################


x = np.asarray(SecondOrderDerivList_rho_min)
y = np.asarray(SecondOrderDerivList_min)

x_slope = np.log(x)
y_slope = np.log(y)

slope, intercept, r, p, se = linregress(x_slope, y_slope)

theorykinks = (x**2/(np.pi**2 * PAP)) * np.log(8*x/(Phi*np.pi**2))

plt.figure()
plt.loglog(SecondOrderDerivList_rho_min,SecondOrderDerivList_min,label="Data")
plt.loglog(x,theorykinks,"--",label='Theory')
plt.legend(loc='upper left')
plt.grid()
plt.xlabel("rho")
plt.ylabel("eta at Second order min")
plt.title("How Second order min changes with rho")

plt.savefig(NoTheoryDir + "/UpperKink_ChangingRho.png")
plt.close()


############################################################################
############################################################################
############################################################################



fig,ax=plt.subplots()
ax.plot(rhoList*(1-rhoList),MinimumList)
ax.set_xlabel("rho(1-rho)")
ax.set_ylabel("eta at minimum")
ax.grid()

for i in range(len(rhoList)):
    ax.annotate("rho: %0.2f"%(rhoList[i]),(np.asarray(rhoList*(1-rhoList))[i],MinimumList[i]))

plt.savefig(NoTheoryDir + "/MinimumChangingRho.png")

plt.close()

fig,ax=plt.subplots()
ax.semilogy(rhoList*(1-rhoList),MinimumList)
ax.set_xlabel("rho(1-rho)")
ax.set_ylabel("eta at minimum")
ax.grid()

for i in range(len(rhoList)):
    ax.annotate("rho: %0.2f"%(rhoList[i]),(np.asarray(rhoList*(1-rhoList))[i],MinimumList[i]))

plt.savefig(NoTheoryDir + "/MinimumChangingRhoSemilogy.png")

plt.close()

fig,ax=plt.subplots()
ax.loglog(rhoList*(1-rhoList),MinimumList)
ax.set_xlabel("rho(1-rho)")
ax.set_ylabel("eta at minimum")
ax.grid()

for i in range(len(rhoList)):
    ax.annotate("rho: %0.2f"%(rhoList[i]),(np.asarray(rhoList*(1-rhoList))[i],MinimumList[i]))

plt.savefig(NoTheoryDir + "/MinimumChangingRhoLoglog.png")

plt.close()


#############################################################################
#############################################################################
#############################################################################
t = 1-PAP

testtheorylist = ((rhoList**2 * (1-rhoList)**2) / (16*np.pi**2 * PAP*t)*
                    special.erfinv(1-Phi)**-2 *
                    np.log(8*rhoList/(Phi*np.pi**2))                    
                )**0.5

"""
MAG = 1e5
def fitfunc(x,a,b,c,d):
    return a * (MAG-x)**b * (x*np.log(x*c))**d

popt, pcov = curve_fit(fitfunc, MAG*np.asarray(rhoList), MAG*np.asarray(MinimumList),maxfev=10000)
"""




plt.figure()
plt.loglog(rhoList,MinimumList,label='data')#,label="Slope = %0.3f"%(result.slope))
plt.loglog(rhoList,testtheorylist,label='Theory')
#plt.loglog(rhoList,fitfunc(MAG*rhoList,*popt)/MAG,label='Theory Fit')
plt.legend(loc='lower left')

plt.xlabel("Refuge Proportion")
plt.ylabel("Eta at Minimum")

plt.title("PAP: %0.5f"%(PAP))
plt.grid()
plt.savefig(NoTheoryDir + "/ComparingMinima.png")
plt.close()

endtime = time.time()

print("Time taken,",endtime-starttime)

print("MINIMALIST RHO",rhoList.tolist())
print("MINIMALIST ETA",MinimumList)






##################################
#Scatter Plot

BoundRhoList = np.linspace(min(HeatMap_rhoList),max(HeatMap_rhoList),100)
x = BoundRhoList
LowerKink = (1-x)**2/(16*(1-PAP) * special.erfinv(1-Phi)**2)
UpperKink = (x**2/(np.pi**2 * PAP)) * np.log(8*x/(Phi*np.pi**2))



HeatMap_rhoList = np.log(np.asarray(HeatMap_rhoList))
HeatMap_etaList = np.log(np.asarray(HeatMap_etaList))
HeatMap_slopeList = np.log(np.asarray(HeatMap_slopeList))

cp = ax.tricontourf(HeatMap_rhoList.ravel(), HeatMap_etaList.ravel(), HeatMap_slopeList.ravel(),10,cmap='coolwarm')




fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.scatter(HeatMap_rhoList,HeatMap_etaList,c=HeatMap_slopeList,s=150,cmap='coolwarm')

plt.plot(np.log(BoundRhoList),np.log(UpperKink),label='Upper Kink')
plt.plot(np.log(BoundRhoList),np.log(LowerKink),label='Lower Kink')

#ax.set_yscale("log")
#ax.set_xscale("log")

plt.legend(loc='lower left',fontsize=20)

ax.set_ylabel('LOG(ETA)',fontsize=15)
ax.set_xlabel('LOG(RHO)',fontsize=15)


cbar = plt.colorbar(cp)
cbar.set_label(label='LOG(SLOPE)',size=15)

cbar.ax.tick_params(labelsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)

plt.xticks(rotation=45)
plt.yticks(rotation=45)


plt.tight_layout()
plt.grid(True)

plt.savefig(NoTheoryDir + '/Scatter.png')
plt.close()

