from numpy import *
from typing import Callable
from .Utility import getSlewRate_Pix

def genSpiral_DeltaK(getDeltaK:Callable, getDrhoDtht:Callable, phase:int|float=0, rhoMax:float|int=0.5) -> ndarray:
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
        quoDrhoDtht = getDrhoDtht(lstRho[-1], lstTht[-1])
        
        dTht = dK/sqrt(quoDrhoDtht**2 + lstRho[-1]**2)
        dRho = dTht*quoDrhoDtht
        
        # append new point
        thtNew = lstTht[-1]+dTht
        rhoNew = lstRho[-1]+dRho
        print(f"rhoNew={rhoNew}")
        if(rhoNew < rhoMax):
            lstTht = append(lstTht, thtNew)
            lstRho = append(lstRho, rhoNew)
        else:
            break

    lstKx = lstRho*cos(lstTht + phase)
    lstKy = lstRho*sin(lstTht + phase)
    lstKxKy = array([lstKx, lstKy]).T

    return lstKxKy

def genSpiral_Slewrate(getSlewRate:Callable, getDrhoDtht:Callable, dt:int|float, phase:int|float=0, rhoMax:float|int=0.5, gamma:float|int=42.58e6) -> tuple[ndarray, ndarray, ndarray]:
    lstTht = array([0, 2*pi])
    lstRho = array([0, gamma*(getSlewRate(0, 0)*dt)*dt/2])
    lstKx = lstRho*cos(lstTht + phase)
    lstKy = lstRho*sin(lstTht + phase)
    lstDk = lstRho.copy()
    grad = getSlewRate(0, 0)*dt*array([cos(phase), sin(phase)])
    while True:
        dK = 3*lstDk[-1]
        limSlewRate = getSlewRate(lstRho[-1], lstTht[-1])
        quoDrhoDtht = getDrhoDtht(lstRho[-1], lstTht[-1])
        while True:
            dTht = dK/sqrt(quoDrhoDtht**2 + lstRho[-1]**2)
            dRho = dTht*quoDrhoDtht
            
            thtNew = lstTht[-1]+dTht
            rhoNew = lstRho[-1]+dRho
            kxNew = rhoNew*cos(thtNew + phase)
            kyNew = rhoNew*sin(thtNew + phase)

            sr = getSlewRate_Pix(array([kxNew, kyNew]), array([lstKx[-1], lstKy[-1]]), grad, dt)
            sr = sqrt(sr[0]**2 + sr[1]**2)

            errSlewRate = sr - limSlewRate
            if abs(errSlewRate) < 1e-1*limSlewRate:
                break
            else:
                print(f"sr={sr}, lim={limSlewRate}")
                biasDk = -1e-3*sign(errSlewRate)*sqrt(abs(errSlewRate))
                dK += biasDk

        # append new point
        print("")
        print(f"rhoNew={rhoNew}")
        if(rhoNew < rhoMax):
            lstTht = append(lstTht, thtNew)
            lstRho = append(lstRho, rhoNew)
            lstKx = append(lstKx, kxNew)
            lstKy = append(lstKy, kyNew)
            lstDk = append(lstDk, dK)

            gradMean = array([lstKx[-1]-lstKx[-2], lstKy[-1]-lstKy[-2]])/gamma/dt
            grad = grad + 2*(gradMean - grad)
        else:
            break

    lstKxKy = array([lstKx, lstKy]).T

    return lstKxKy, lstRho, lstDk

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
    lstKx = [lstRho*cos(lstTht[idxTht]) for idxTht in range(lstTht.size)]
    lstKy = [lstRho*sin(lstTht[idxTht]) for idxTht in range(lstTht.size)]

    return array([lstKx, lstKy]).transpose([1,2,0])

def genCart(numPt:int|float, max:int|float=0.5) -> ndarray:
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
    