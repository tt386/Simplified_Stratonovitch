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

    theoryendofyear = 2*((1-rho)/2 + (rho/2)*(Phi/(Phi+1)) + b*np.sqrt(np.pi)*(Phi*np.log((1+Phi)/Phi) - 1/2 - (Phi/(Phi+1))/2))

    theoryList = 2*((1-rho)/2 + (rho/2)*(Phi/(Phi+Phiprime)) + b*np.sqrt(np.pi)*(Phi*np.log((Phiprime+Phi)/Phi) - 1/2 - (Phi/(Phi+Phiprime))/2))#np.log(2*((1-rho)/2 + (rho/2)*(Phi/(Phi+Phiprime)) + b*np.sqrt(np.pi)*(Phi*np.log((Phiprime+Phi)/Phi) - 1/2 - (Phi/(Phi+Phiprime))/2))) - np.log(Phi)

    

    AnalyticalSlopeList = []
    for eta in etaList:
        b = np.sqrt(4*eta*(1-PAP))
        AnalyticalSlopeList.append(np.log(2*integrate.quad(lambda x: Phi/(Phi + 0.5+0.5*special.erf((rho-x)/b) + 0.5 * (1-special.erf((1-x)/b))), rho/2, (1+rho)/2))-np.log(Phi))

    print("Length off AnalyticalSlopeList",len(AnalyticalSlopeList))



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


    LargeEtaSlopeList = []
    SmallEtaSlopeList =  []
    for eta in etaList:
        b = np.sqrt(4*eta*(1-PAP))

        Theoryxlist = np.linspace(rho/2,(1+rho)/2,1000)

        LargeEtaSlopeList.append(2*integrate.simps(LargeEtaIntegrand(b,rho,Theoryxlist,Phi),Theoryxlist))#(np.log(2*integrate.simps(LargeEtaIntegrand(b,rho,Theoryxlist,Phi),Theoryxlist)) - np.log(Phi))
        SmallEtaSlopeList.append(2*integrate.simps(SmallEtaIntegrand(b,rho,Theoryxlist,Phi),Theoryxlist))#(np.log(2*integrate.simps(SmallEtaIntegrand(b,rho,Theoryxlist,Phi),Theoryxlist)) - np.log(Phi))
    #########################################################################
    #########################################################################
    #########################################################################

    #Change from 'Slope' to endstate number
    slopelist = np.exp(slopelist + np.log(Phi))

    plt.figure()
    plt.loglog(etaList,slopelist,label="Data")
    #plt.semilogy(etaList,minmatrix[c],label="Min Point")
    plt.loglog(etaList,SmallEtaSlopeList,label='Small Eta')
    plt.loglog(etaList,theoryList,label="Linearised Small Eta")
    plt.loglog(etaList,LargeEtaSlopeList,label='Large Eta')

    """
    #Plot minima and maxima
    val, idx = max((val, idx) for (idx, val) in enumerate(minmatrix[c]))
    plt.semilogy(etaList[idx],val,"ro")

    max_minpointlist.append(etaList[idx])

    val, idx = min((val, idx) for (idx, val) in enumerate(slopelist))
    plt.semilogy(etaList[idx],val,"bo")
    """

    plt.legend(loc='lower left')

    plt.xlabel("eta")
    plt.ylabel("Log of initial yearly R Allele Slope")

    plt.title("Refuge Proportion " + '{:.1E}'.format(Refuge_Proportion))
    plt.grid()
    plt.savefig(str(args.directory) +
        "/Minima_Curve_rho_" + '{:.1E}'.format(Refuge_Proportion) + ".png")
    plt.close()
    

    #Plot of slopes of the curves
    plt.figure()
    plt.loglog(etaList,abs(np.gradient(slopelist)),label="Data")
    plt.loglog(etaList,abs(np.gradient(SmallEtaSlopeList)),label='Small Eta')
    plt.loglog(etaList,abs(np.gradient(theoryList)),label="Linearised Small Eta")
    plt.loglog(etaList,abs(np.gradient(LargeEtaSlopeList)),label='Large Eta')

    plt.legend(loc='lower left')

    plt.xlabel("eta")
    plt.ylabel("Slope")

    plt.title("Refuge Proportion " + '{:.1E}'.format(Refuge_Proportion))

    plt.grid()

    plt.savefig(str(args.directory) +
        "/Gradient_Curve_rho_" + '{:.1E}'.format(Refuge_Proportion) + ".png")
    plt.close()


    #Second order derivative
    plt.figure()
    plt.loglog(etaList,abs(np.gradient(np.gradient(slopelist))),label="Data")
    plt.loglog(etaList,abs(np.gradient(np.gradient(SmallEtaSlopeList))),label='Small Eta')
    plt.loglog(etaList,abs(np.gradient(np.gradient(theoryList))),label="Linearised Small Eta")
    plt.loglog(etaList,abs(np.gradient(np.gradient(LargeEtaSlopeList))),label='Large Eta')

    plt.legend(loc='lower left')

    plt.title("Refuge Proportion " + '{:.1E}'.format(Refuge_Proportion))

    plt.xlabel("eta")
    plt.ylabel("Second Derivative")

    plt.grid()

    plt.savefig(str(args.directory) +
        "/SecondOrderDeriv_Curve_rho_" + '{:.1E}'.format(Refuge_Proportion) + ".png")
    plt.close()


    #finding max points
    if rho > 0:
        SecondOrderDerivMaxList = [argrelextrema(abs(np.gradient(np.gradient(slopelist))), np.greater)]
        print("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
        print(SecondOrderDerivMaxList)
        print(SecondOrderDerivMaxList[0])
        print(SecondOrderDerivMaxList[0][0])
        print("BBBBBBBBBBBBBBBBBBBBBBBBBBB")
        try:
            print(SecondOrderDerivMaxList)
            print(SecondOrderDerivMaxList[0][0][1])
            SecondOrderDerivList.append(etaList[SecondOrderDerivMaxList[0][0][1]])
            SecondOrderDerivList_rho.append(rho)
        except:
            print("no maxima found")

        print("SecondOrderDerivList:",SecondOrderDerivList)


    #finding min points
    if rho > 0:
        SecondOrderDerivMinList = [argrelextrema(abs(np.gradient(np.gradient(slopelist))), np.less)]
        print("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
        print("Second order min:",SecondOrderDerivMinList)
        print("Second order min:",SecondOrderDerivMinList[0])
        print("Second order min:",SecondOrderDerivMinList[0][0])
        print("BBBBBBBBBBBBBBBBBBBBBBBBBBB")
        try:
            print(SecondOrderDerivMinList)
            print(SecondOrderDerivMinList[0][0][0])
            SecondOrderDerivList_min.append(etaList[SecondOrderDerivMinList[0][0][0]])
            SecondOrderDerivList_rho_min.append(rho)
        except:
            print("no minima found")

        print("SecondOrderDerivList:",SecondOrderDerivList)


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




x = 1.-np.asarray(SecondOrderDerivList_rho)
y = np.asarray(SecondOrderDerivList)

x_slope = np.log(x[x>0.3])
y_slope = np.log(y[x>0.3])

slope, intercept, r, p, se = linregress(x_slope, y_slope)

plt.figure()
plt.loglog(x,y,label="Slope: %0.5f"%(slope))
plt.legend(loc='upper left')
plt.grid()
plt.xlabel("1-rho")
plt.ylabel("eta at kink")
plt.title("How initial kink changes with rho")

plt.xlim(0.3,1)
plt.savefig(str(args.directory) + "/KinkLocations.png")
plt.close()




x = np.asarray(SecondOrderDerivList_rho_min)
y = np.asarray(SecondOrderDerivList_min)

x_slope = np.log(x)
y_slope = np.log(y)

slope, intercept, r, p, se = linregress(x_slope, y_slope)

plt.figure()
plt.loglog(SecondOrderDerivList_rho_min,SecondOrderDerivList_min,label="Slope: %0.5f"%(slope))
plt.legend(loc='upper left')
plt.grid()
plt.xlabel("rho")
plt.ylabel("eta at Second order min")
plt.title("How Second order min changes with rho")

plt.savefig(str(args.directory) + "/SecondOrderMinimumKink.png")
plt.close()






fig,ax=plt.subplots()
ax.plot(rhoList*(1-rhoList),MinimumList)
ax.set_xlabel("rho(1-rho)")
ax.set_ylabel("eta at minimum")
ax.grid()

for i in range(len(rhoList)):
    ax.annotate("rho: %0.2f"%(rhoList[i]),(np.asarray(rhoList*(1-rhoList))[i],MinimumList[i]))

plt.savefig(str(args.directory) + "/MinimumChangingRho.png")

plt.close()

fig,ax=plt.subplots()
ax.semilogy(rhoList*(1-rhoList),MinimumList)
ax.set_xlabel("rho(1-rho)")
ax.set_ylabel("eta at minimum")
ax.grid()

for i in range(len(rhoList)):
    ax.annotate("rho: %0.2f"%(rhoList[i]),(np.asarray(rhoList*(1-rhoList))[i],MinimumList[i]))

plt.savefig(str(args.directory) + "/MinimumChangingRhoSemilogy.png")

plt.close()

fig,ax=plt.subplots()
ax.loglog(rhoList*(1-rhoList),MinimumList)
ax.set_xlabel("rho(1-rho)")
ax.set_ylabel("eta at minimum")
ax.grid()

for i in range(len(rhoList)):
    ax.annotate("rho: %0.2f"%(rhoList[i]),(np.asarray(rhoList*(1-rhoList))[i],MinimumList[i]))

plt.savefig(str(args.directory) + "/MinimumChangingRhoLoglog.png")

plt.close()


"""
plt.figure()
plt.plot(MinimumList,max_minpointlist,"ro")
plt.xlabel("eta at min slope")
plt.ylabel("eta at max minpoint")
plt.grid()
plt.savefig(str(args.directory) + "/Correlation.png")

plt.close()
"""
"""
#Fitting
MAGNTIDUEINCREASE = 10**10
def fit_func(x,b,c,d):
    return  b *( (1*MAGNTIDUEINCREASE-x)**c * x**d)

def fit_func2(x,a,b,c,d,e):
    return a*(MAGNTIDUEINCREASE-x)**b + c*x**d + e


def fit_func3(x,a,b,c):
    return a * ()


firstplateau = 0
lastplateau = 0

for i in range(len(rhoList)):
    if rhoList[i] > 0.01:
        firstplateau = i
        break

for i in range(len(rhoList)):
    if rhoList[i] > 0.98:
        lastplateau = i
        break
    

fitMinimumList = np.asarray(MinimumList[firstplateau:lastplateau])*MAGNTIDUEINCREASE
fitrhoList = rhoList[firstplateau:lastplateau]*MAGNTIDUEINCREASE

sigmaweights = np.ones(len(fitrhoList))
#sigmaweights[:20] = 0.01
#sigmaweights[20:] = 0.01
sigmaweights*=0.01
popt,pcov = curve_fit(fit_func,fitrhoList,fitMinimumList,maxfev=100000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
b,c,d = popt


fitrhoList = fitrhoList/MAGNTIDUEINCREASE
fittedmin1 =  b*( (1 - fitrhoList)**c * fitrhoList**d) * MAGNTIDUEINCREASE**(c+d-1)

print("DifferenceList: ",fittedmin1 - fitMinimumList/MAGNTIDUEINCREASE)

#Second fit
fitrhoList*=MAGNTIDUEINCREASE
popt,pcov = curve_fit(fit_func2,fitrhoList,fitMinimumList,maxfev=100000000,sigma=sigmaweights)

print(popt)
print("Errors:",np.sqrt(np.diag(pcov)))
a,b,c,d,e = popt

fitrhoList = fitrhoList/MAGNTIDUEINCREASE
fittedmin2 = a * MAGNTIDUEINCREASE**(b-1) * (1-fitrhoList)**b + c*MAGNTIDUEINCREASE**(d-1) * fitrhoList**d + e*MAGNTIDUEINCREASE**-1
"""

t = 1-PAP

testtheorylist = (1-np.sqrt(1+(np.pi/rhoList)**2 *PAP* (1-rhoList)**2 / t))/(-4*(np.pi/rhoList)**2 * PAP)

plt.figure()
plt.loglog(rhoList,MinimumList,label='data')#,label="Slope = %0.3f"%(result.slope))
plt.loglog(rhoList,testtheorylist,label='Theory')
#plt.loglog(fitrhoList,fittedmin1,label='fit: product')
#plt.loglog(fitOSRRatioList,fittedmin2,label='fit: sum')
plt.legend(loc='lower left')

plt.xlabel("Refuge Proportion")
plt.ylabel("Eta at Minimum")

plt.title("PAP: %0.5f"%(PAP))
plt.grid()
plt.savefig(str(args.directory) + "/ComparingMinima.png")
plt.close()

endtime = time.time()

print("Time taken,",endtime-starttime)
