# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 11:36:56 2022

@author: tt386
"""


import numpy as np
import copy
import matplotlib.pyplot as plt

from scipy import integrate
import sys
##############################################################################
##############################################################################
##############################################################################
#For the very first step
def A(n,SRatio):
    """Returns the Fourier Coefficients for the Dicichlet Heat Eqn"""
    if n%2 == 0:
        return 0
    else:
        #print("First  Step A",n,":",(1-SRatio)*4/(np.pi * n))
        return (1-SRatio)*4/(np.pi * n)

def u1(x,t,l,k,Phi):
    """Numerical result for Dirichlet Heat Eqn"""
    if  (x < l):
        tot = 0
        for n in range(1,100):
            tot += A(n,Phi) * np.sin(n*np.pi*x/l) * np.exp(-k*t*(n*np.pi/l)**2)
            #print(x,"First Step Integral",n,":",A(n,Phi) * np.sin(n*np.pi*x/l) * np.exp(-k*t*(n*np.pi/l)**2))
        return tot
    else:
        return 0
##############################################################################
##############################################################################
##############################################################################
#For the second step
def A0(L,init,xlist):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    print("A0",(1/(2*L))*integrate.simps(init, xlist))
    return (1/(2*L))*integrate.simps(init, xlist)

def An(n,init,xlist,L):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    tempinit = np.asarray(copy.copy(init))

    for x in range(len(xlist)):
        tempinit[x] *= np.cos(n*np.pi*xlist[x]/L)

    return (1/L)*integrate.simps(tempinit, xlist)

def Bn(n,init,xlist,L):
    """Returns the Fourier Coefficients for the Periodic Heat Eqn"""
    tempinit = np.asarray(copy.copy(init))

    for x in range(len(xlist)):
        tempinit[x] *= np.sin(n*np.pi*xlist[x]/L)

    return (1/L)*integrate.simps(tempinit, xlist)

def u2(t,L,k,init,xlist):
    """Numerical result for Periodic Heat Eqn"""
    A0Val = A0(L,init,xlist)
    AnList = []
    BnList = []

    #print("A0",A0Val)

    for n in range(1,200):
        AnList.append(An(n,init,xlist,L))
        BnList.append(Bn(n,init,xlist,L))

    

    resultslist = []
    
    #ratiolist = []
    for x in xlist:
        tot = A0Val

        for n in range(1,len(AnList)):
            tot += (np.exp(-k*t*(np.pi*n/L)**2) *
                (AnList[n-1] * np.cos(n*np.pi*x/L) +
                BnList[n-1] * np.sin(n*np.pi*x/L)))
        #print("x:",x,"A0",A0Val,"Tot:",tot)
        #ratiolist.append(A0Val/tot)
        resultslist.append(tot)
    #print(ratiolist)
    
    #plt.figure()
    #plt.loglog(xlist,ratiolist)
    #plt.grid()
    #plt.show()
    
    return resultslist









Year = 0.1#200#0.1#200
PAP = 0.5
T = Year*(1-PAP)
Phi = 1e-5 #Initial Conditions
k = 1e-5

OSR_Proportion = 0.8

dx =1e-4#1 #1e-4#1
dt =1e-4#1 #1e-4#1

SystemSize = 2e-3

Refuge_Proportion = 1- OSR_Proportion


FigName = ("YearLength_%0.3f,PAP_%0.3f,k_%0.3f,OSRProp_%0.3f_"%
            (Year,PAP,k,OSR_Proportion) + "dt_" + '{:.1E}'.format(dt)
            + "_dx_" + '{:.1E}'.format(dx) + "_SysSize_" +
            '{:.1E}'.format(SystemSize) + ".png")

#########################
print("(Year,System Size) =",Year,SystemSize)

######################################################
######################################################
######################################################
###Analytical
#Create system
OSRWidth = (1-Refuge_Proportion)*SystemSize
RefugeWidth = Refuge_Proportion*SystemSize
xlist = np.arange(int(SystemSize/dx))*dx
ylist = []


#Step 1
#Solve the analytical result for the system of Refuge with Dirichlet
#Boundary Conditions
for x in range(len(xlist)):
    ylist.append(u1(xlist[x],PAP*Year,RefugeWidth,k,Phi))
    #if x >  5: sys.exit()

#Find area, similar to the number of susceptible pests
step1area = integrate.simps(ylist,xlist)


#Steo 2
#Solve the analytical result for the system fo Refuge, OSR and
#Periodic Boundary conditions
L = SystemSize/2
yprimelist = np.asarray(u2(T,L,k,ylist,xlist))

#Ensure no negative values
yprimelist[yprimelist<0] = 0
######################################################
######################################################
######################################################
###Sim
days = int(Year/dt)
MigrationsPerDay = 1

F = k * (dt)/(dx)**2

systemsites = int(SystemSize/dx)

OSRWidth = int(systemsites * OSR_Proportion)
RefugeSize = int(systemsites * (1-OSR_Proportion))

systemsites = OSRWidth + RefugeSize

#Making lists
S_Domain = np.ones(systemsites)*(1-Phi)
R_Domain = np.ones(systemsites)*(Phi)


PesticideDomain = np.ones(OSRWidth)
PesticideDomain = np.pad(PesticideDomain,(0,RefugeSize),mode='constant')

PesticideDomain = np.flip(PesticideDomain)

EndOfPAPDomain = 0

####Process#######################
for t in range(days):
    #print(t)
    if t < PAP*days:
        S_Domain[PesticideDomain == 1] = 0
        EndOfPAPDomain = copy.copy(S_Domain)
    for tm in range(MigrationsPerDay):

        S_Domain = (S_Domain +
                    F*(np.roll(S_Domain,-1) - 2 * S_Domain +
                       np.roll(S_Domain,1)))
        

######################################################
######################################################
######################################################
#Plotting

plt.figure()
plt.plot(xlist,ylist,label='Analytical')
plt.plot(np.arange(len(S_Domain))*dx,EndOfPAPDomain,label='Numerical')
plt.legend(loc='upper center')
plt.grid()
plt.savefig("EndOfPAP_"+FigName)
plt.close()


plt.figure()
plt.plot(xlist,ylist,'--b',label='After PAP: Analytical')
plt.plot(np.arange(len(S_Domain))*dx,EndOfPAPDomain,'--r',label='After PAP: Numerical')

plt.plot(xlist,yprimelist,'b',label='End of Year: Analytical')
plt.plot(np.arange(len(S_Domain))*dx,S_Domain,'r',label='End of Year: Numerical')
plt.legend(loc='upper center')
plt.grid()
plt.savefig("EndOfYear_"+FigName)
