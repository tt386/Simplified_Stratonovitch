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
    etalist = ParamDict["etalist"]
    PAP = ParamDict["PAP"]
    Phi = ParamDict["Phi"]
    xlist = ParamDict["xlist"]
    OSRNum = ParamDict["OSRNum"]
    OSRSeparation = ParamDict["OSRSeparation"]
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
        #if  (x < rho):
        tot = 0

        for n in range(1,N,2):
            tot += ((1-Phi)*4/(np.pi * n) * np.sin(n*np.pi*x/d) *
                    np.exp(-eta*PAP*(n*np.pi/d)**2))
            """
            tot += ((1-Phi)*4/(np.pi * n) * np.sin(n*np.pi*x/rho) *
                    np.exp(-eta*PAP*(n*np.pi/d)**2))
            """
        return tot


    def PAP_Dist(x,eta,t,d,OSRNum,N):
        Upperx = (OSRNum-1)*(1+d) + 1
        
        modx = x%(1+d)
        
        if x < 0:
            return (1-Phi)*special.erf((-x)/np.sqrt(4*eta*t))
        elif x > Upperx:
            return (1-Phi)*special.erf((x-Upperx)/np.sqrt(4*eta*t))
        elif modx > 1:
            return u1((modx-1),eta,d,t,Phi,N)
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
    for eta in etalist:
        print("Phi",Phi,"PAP",PAP,"Sep:",OSRSeparation,", OSRNum",OSRNum,", Eta:",eta)
        PAPylist = []

        #It will be much better to have both for loops put together.
        #for x in xlist:
        #    PAPylist.append(PAP_Dist(x,eta,PAP,OSRSeparation,OSRNum,N)) 
            
        
        T = 1.-PAP
        Yearylist = []

        for x in xlist:#[:len(xlist)//2]:
            PAPylist.append(PAP_Dist(x,eta,PAP,OSRSeparation,OSRNum,N))

            #print(eta,PAP,OSRSeparation,OSRNum,T)
            combinedfun = (lambda r: 
                PAP_Dist(r,eta,PAP,OSRSeparation,OSRNum,N) * Sense(x-r,eta,T))
            
            #Split integral into upper and lower to eradiate error
            upper = integrate.quad(combinedfun,x,np.inf)[0]#integrate.quad(combinedfun,x,np.inf)[0]
            lower = integrate.quad(combinedfun,-np.inf,x)[0]#integrate.quad(combinedfun,-np.inf,x)[0]
            Yearylist.append(upper+lower)
            
            
        Yearylist = np.asarray(Yearylist)
           
        rlist = np.ones(len(xlist)) * Phi 
        
        PAPMatrix.append(PAPylist)
        YearMatrix.append(Yearylist)
       
        Integrand = rlist/(rlist+Yearylist)-rlist

        #Search for the lower limit which is the first time we
        # fall below 10^-19 (beyond here we get much error)
        UpperLim = len(xlist)
        for x in range(len(xlist)):
            if xlist[x] > (OSRNum + OSRSeparation*(OSRNum-1)):
                if Integrand[x] < 1e-19:
                    UpperLim = x
                    break
    

        EndYearRatioList.append(
            2*integrate.simps(Integrand[:UpperLim],xlist[:UpperLim]))
        
    ###################################################
    ###Package and Return Results######################
    ###################################################
    ResultsDict = {
        "PAPMatrix": PAPMatrix,
        "YearMatrix": YearMatrix,
        "RDifference": EndYearRatioList
        }

    return ResultsDict


