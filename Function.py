import numpy as np
from typing import Callable

def genSpiral(dSamp:float|Callable, kRhoTht:float|Callable, phase:float=0, rhoMax:float=0.5) -> np.ndarray:
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
        
    lstTht = np.array([0, 2*np.pi])
    lstRho = np.array([0, rho1])
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
        
        dTht = dSamp_This/np.sqrt(kRhoTht_This**2 + lstRho[-1]**2)
        dRho = dTht*kRhoTht_This
        
        thtNew = lstTht[-1]+dTht
        rhoNew = lstRho[-1]+dRho
        # append new point
        if(rhoNew < rhoMax):
            lstTht = np.append(lstTht, thtNew)
            lstRho = np.append(lstRho, rhoNew)
        else:
            break

    lstKx = lstRho*np.cos(lstTht + phase)
    lstKy = lstRho*np.sin(lstTht + phase)

    return np.array([lstKx, lstKy]).T.copy()#

def genRadial(lstTht:np.ndarray, lstRho:np.ndarray) -> np.ndarray:
    # shape check
    assert(np.size(lstTht.shape) == 1)
    assert(np.size(lstRho.shape) == 1)
    
    # generate kspace trajectory
    lstKx = np.zeros([lstTht.size*lstRho.size], dtype=np.float64)
    lstKy = np.zeros([lstTht.size*lstRho.size], dtype=np.float64)
    idxPt = 0
    for idxTht in range(lstTht.size):
        lstKx[idxPt:idxPt+lstRho.size] = lstRho*np.cos(lstTht[idxTht])
        lstKy[idxPt:idxPt+lstRho.size] = lstRho*np.sin(lstTht[idxTht])
        idxPt += lstRho.size

    return np.array([lstKx, lstKy]).T.copy()