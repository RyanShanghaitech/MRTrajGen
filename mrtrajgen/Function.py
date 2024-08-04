from numpy import *
from typing import Callable

def genSpiral(getDeltaK:Callable, getDrhoDtht:Callable, phase:float|int=0, rhoMax:float|int=0.5) -> ndarray:
    """
    # description:
    generate spiral sampling trajectory

    # parameter:
    `getDeltaK`:Callable: function of sampling interval with respect to rho and theta, e.g. `lambda rho, tht: dNyq` for spiral with constant Nyquist sampling interval
    `getDrhoDtht`:Callable: function of dRho/dTheta with respect to rho and theta, e.g. `lambda rho, tht: b` for Archimedean spiral, `lambda rho, tht: rho + 1` for variable density spiral
    `phase`: phase controlling rotation of spiral
    `rhoMax`: maximum value of rho, 0.5 covers the whole kspace

    # return:
    kspace trajectory: [[kx1, ky1], [kx2, ky2], ..., [kxn, kyn]]
    """
    lstTht = array([0, 2*pi])
    lstRho = array([0, getDeltaK(0, 0)])
    while True:
        dK = getDeltaK(lstRho[-1], lstTht[-1])
        dRhoTht = getDrhoDtht(lstRho[-1], lstTht[-1])
        
        dTht = dK/sqrt(dRhoTht**2 + lstRho[-1]**2)
        dRho = dTht*dRhoTht
        
        # append new point
        thtNew = lstTht[-1]+dTht
        rhoNew = lstRho[-1]+dRho
        if(rhoNew < rhoMax):
            lstTht = append(lstTht, thtNew)
            lstRho = append(lstRho, rhoNew)
        else:
            break

    lstKx = lstRho*cos(lstTht + phase)
    lstKy = lstRho*sin(lstTht + phase)

    return array([lstKx, lstKy]).T.copy()

def genRadial(lstTht:ndarray, lstRho:ndarray) -> ndarray:
    """
    # description:
    generate radial sampling trajectory

    # parameter:
    `lstTht`: list of theta of spokes
    `lstRho`: list of rho of spokes

    # return:
    kspace trajectory: [[kx1, ky1], [kx2, ky2], ..., [kxn, kyn]]
    """
    # shape check
    assert(size(lstTht.shape) == 1)
    assert(size(lstRho.shape) == 1)
    
    # generate kspace trajectory
    lstKx = zeros([lstTht.size, lstRho.size], dtype=float64)
    lstKy = zeros([lstTht.size, lstRho.size], dtype=float64)
    for idxTht in range(lstTht.size):
        lstKx[idxTht,:] = lstRho*cos(lstTht[idxTht])
        lstKy[idxTht,:] = lstRho*sin(lstTht[idxTht])

    return array([lstKx, lstKy]).transpose([1,2,0])

def genCartesian(numPt:int|float, max:int|float=0.5) -> ndarray:
    """
    # description:
    generate cartesian sampling trajectory

    # parameter:
    `numPt`: number of point in one dimension
    `max`: maximum coordinate value, 0.5 for kspace

    # return:
    trajectory: [[x1, y1], [x2, y2], ..., [xn, yn]]
    """
    # shape check
    lstKx, lstKy = meshgrid(
        linspace(-max, max, numPt, endpoint=False),
        linspace(-max, max, numPt, endpoint=False))
    return array([lstKx.flatten(), lstKy.flatten()]).T
    