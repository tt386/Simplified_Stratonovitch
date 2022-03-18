import numpy as np
import matplotlib.pyplot as plt

from scipy import integrate
from scipy import optimize
from scipy import special

from scipy.optimize import curve_fit

import copy

import sys
#############################################
#############################################
#############################################


def Core(ParamDict):
    #############################
    ###Unpack ParamDict##########
    #############################
    etalist = ParamDict["etalist"]
    PAP = ParamDict["PAP"]
    Phi = ParamDict["Phi"]
    xlist = ParamDict["xlist"]
    OSRNum = ParamDict["OSRNum"]
    OSRSeparation = ParamDict["OSRSeparation"]

    #############################
    ###Functions#################
    #############################
    #Function for PAP Dist
    def A(n,SRatio):
        """Returns the Fourier Coefficients for the Dicichlet Heat Eqn"""
        if n%2 == 0:
            return 0
        else:
            return (1-SRatio)*4/(np.pi * n)

    def u1(x,eta,rho,PAP,Phi,N):
        """Numerical result for Dirichlet Heat Eqn"""
        if  (x < rho):
            tot = 0
            for n in range(1,N):
                tot += (A(n,Phi) * np.sin(n*np.pi*x/rho) *
                        np.exp(-eta*PAP*(n*np.pi/rho)**2))
            return tot
        else:
            return 0
        
    def PAP_Dist(x,k,t,d,OSRNum):
        Upperx = (OSRNum-1)*(1+d) + 1
        
        modx = x%(1+d)
        
        if x < 0:
            return special.erf((-x)/np.sqrt(4*k*t))
        elif x > Upperx:
            return special.erf((x-Upperx)/np.sqrt(4*k*t))
        elif modx > 1:
            rho = d/(1+d)
            return u1((modx-1)/(d),k,1,t,Phi,50)
        else:
            return 0

    def Sense(r,k,t):
        return (1./np.sqrt(4*np.pi*k*t)*np.exp(-r**2/(4*k*t)))



    ###################################################
    ###Main Process####################################
    ###################################################
    PAPMatrix = []
    YearMatrix = []
    EndYearRatioList = []
    for eta in etalist:
        print("Sep:",OSRSeparation,", Eta:",eta)
        PAPylist = []
        for x in xlist:
            PAPylist.append(PAP_Dist(x,eta,PAP,OSRSeparation,OSRNum)) 
            
        
        T = 1-PAP
        Yearylist = []
        midx = (OSRNum + (OSRNum-1)*OSRSeparation)/2
        for x in xlist:
            combinedfun = (lambda r: 
                PAP_Dist(r,eta,PAP,OSRSeparation,OSRNum) * Sense(x-r,eta,T))
            
            #Split integral into upper and lower to eradiate error
            upper = integrate.quad(combinedfun,x,np.inf)[0]
            lower = integrate.quad(combinedfun,-np.inf,x)[0]
            Yearylist.append(upper+lower)
            
            
        Yearylist = np.asarray(Yearylist)
           
        rlist = np.ones(len(xlist)) * Phi 
        
        PAPMatrix.append(PAPylist)
        YearMatrix.append(Yearylist)
        
        EndYearRatioList.append(
            integrate.simps(rlist/(rlist+Yearylist) - rlist,xlist))
        
    ###################################################
    ###Package and Return Results######################
    ###################################################
    ResultsDict = {
        "PAPMatrix": PAPMatrix,
        "YearMatrix": YearMatrix,
        "RDifference": EndYearRatioList
        }

    return ResultsDict


