from numpy import *
from typing import Callable

def funGenSpiral(dSamp:float|Callable, kRhoTht:float|Callable, phase:float=0, rhoMax:float=0.5) -> ndarray:
# dSamp:float: sampling interval
# dSamp:Callable: function of sampling interval with respect to rho and theta, e.g. lambda rho, tht: rho + 1
# kRhoTht:float: ratio of rho/theta
# kRhoTht:Callable: function of ratio of rho/theta with respect to rho and theta, e.g. lambda rho, tht: rho + 1
# phase: phase controlling rotation of spiral
# rhoMax: maximum value of rho, 0.5 covers the whole kspace
# return: kspace trajectory: [[kx1, ky1], [kx2, ky2], ..., [kxn, kyn]]
    if isinstance(dSamp, float):
        rho1 = dSamp
    elif isinstance(dSamp, Callable):
        rho1 = dSamp(0, 0)
    else:
        raise TypeError("dSamp:float|Callable")
        
    lstTht = array([0, 2*pi])
    lstRho = array([0, rho1])
    while True:
        if isinstance(dSamp, float):
            dSamp_This = dSamp
        elif isinstance(dSamp, Callable):
            dSamp_This = dSamp(lstRho[-1], lstTht[-1])
        else:
            raise TypeError("dSamp:float|Callable")
        
        if isinstance(kRhoTht, float):
            kRhoTht_This = kRhoTht
        elif isinstance(kRhoTht, Callable):
            kRhoTht_This = kRhoTht(lstRho[-1], lstTht[-1])
        else:
            raise TypeError("kRhoTht:float|Callable")
        
        #        |\
        #        | \
        #    dRho|  \dSamp
        # =k*dTht|   \
        #        ------
        #    lstRho[-1]*dTht
        
        dTht = dSamp_This/sqrt(kRhoTht_This**2 + lstRho[-1]**2)
        dRho = dTht*kRhoTht_This
        
        thtNew = lstTht[-1]+dTht
        rhoNew = lstRho[-1]+dRho
        # append new point
        if(rhoNew < rhoMax):
            lstTht = append(lstTht, thtNew)
            lstRho = append(lstRho, rhoNew)
        else:
            break

    lstKx = lstRho*cos(lstTht + phase)
    lstKy = lstRho*sin(lstTht + phase)

    return array([lstKx, lstKy]).T.copy()#
