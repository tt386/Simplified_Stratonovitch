import numpy as np
import matplotlib.pyplot as plt

from scipy import integrate
from scipy import optimize
from scipy import special

from scipy.optimize import curve_fit

import copy

import sys

import functools
#############################################
#############################################
#############################################


def Core(ParamDict):
    #############################
    ###Unpack ParamDict##########
    #############################
    L = ParamDict["L"]
    eta = ParamDict["eta"]
    PAP = ParamDict["PAP"]
    Phi = ParamDict["Phi"]
    xlist = ParamDict["xlist"]
    OSRNum = ParamDict["OSRNum"]
    OSRSeparationList = ParamDict["OSRSeparationList"]
    N = ParamDict["N"]

    #############################
    ###Functions#################
    #############################
    #Function for PAP Dist
    def A(n,RRatio):
        """Returns the Fourier Coefficients for the Dicichlet Heat Eqn"""
        if n%2 == 0:
            return 0
        else:
            return (1-RRatio)*4/(np.pi * n)

    @functools.lru_cache()
    def u1(x,eta,d,PAP,Phi,N):
        """Numerical result for Dirichlet Heat Eqn"""
        tot = 0
        """
        for n in range(1,N):
            tot += (A(n,Phi) * np.sin(n*np.pi*x/d) *
                    np.exp(-eta*PAP*(n*np.pi/d)**2))
        """
        for n in range(1,N,2):
            tot += ((1-Phi)*4/(np.pi * n) * np.sin(n*np.pi*x/d) *
                    np.exp(-eta*PAP*(n*np.pi/d)**2))


        return tot








    def PAP_Dist(x,L,eta,t,d,OSRNum,N):
        L = 1

        W = L/OSRNum    #Width of OSR

        Upperx = L + (OSRNum-1)*d#(OSRNum-1)*(1+d) + 1
        
        modx = x%(W + d)
        
        if x < 0:
            return (1-Phi)*special.erf((-x)/np.sqrt(4*eta*t))
        elif x > Upperx:
            return (1-Phi)*special.erf((x-Upperx)/np.sqrt(4*eta*t))
        elif modx > W:
            return u1((modx-W),eta,d,t,Phi,N)
        else:
            return 0

    def Sense(r,eta,t):
        return (1./np.sqrt(4*np.pi*eta*t)*np.exp(-r**2/(4*eta*t)))



    ###################################################
    ###Main Process####################################
    ###################################################
    PAPMatrix = []
    YearMatrix = []
    EndYearRatioList = []
    for OSRSeparation in OSRSeparationList:
        print("Sep:",OSRSeparation,", OSRNum",OSRNum,", Eta:",eta)
        PAPylist = []
        for x in xlist:
            PAPylist.append(PAP_Dist(x,L,eta,PAP,OSRSeparation,OSRNum,N)) 
            
        
        T = 1-PAP
        Yearylist = []
        for x in xlist:
            combinedfun = (lambda r: 
                PAP_Dist(r,L,eta,PAP,OSRSeparation,OSRNum,N) * Sense(x-r,eta,T))
            
            #Split integral into upper and lower to eradiate error
            upper = integrate.quad(combinedfun,x,np.inf)[0]
            lower = integrate.quad(combinedfun,-np.inf,x)[0]
            Yearylist.append(upper+lower)
            
            
        Yearylist = np.asarray(Yearylist)
           
        rlist = np.ones(len(xlist)) * Phi 
        
        PAPMatrix.append(PAPylist)
        YearMatrix.append(Yearylist)
        
        EndYearRatioList.append(
            integrate.simps(rlist/(rlist+Yearylist)-rlist,xlist))
        
    ###################################################
    ###Package and Return Results######################
    ###################################################
    ResultsDict = {
        "PAPMatrix": PAPMatrix,
        "YearMatrix": YearMatrix,
        "RDifference": EndYearRatioList
        }

    return ResultsDict


