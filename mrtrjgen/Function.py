from numpy import *
from typing import Callable
from .Utility import tranGrad2Traj_MinSR
from . import ext
import math

def genSpiral(getD0RhoTht:Callable, getD1RhoTht:Callable, getD2RhoTht:Callable, srlim:int|float, glim:int|float, dt:int|float, kmax:int|float, oversamp:int=2, flagDebugInfo:bool=False) -> tuple[ndarray, ndarray]:
    '''
    # description
    generate spiral trajectory, subject to slew rate

    # parameter
    `getD0RhoTht`: function of 0th order derivation of rho with respect to theta, in `/pix`
    `getD1RhoTht`: function of 1st order derivation of rho with respect to theta, in `/pix/rad`
    `getD2RhoTht`: function of 2nd order derivation of rho with respect to theta, in `/pix/rad^2`
    `srlim`: slew rate limit, in `Hz/pix/s`
    `grad`: gradient limit, in `Hz/pix`
    `dt`: time between 2 adjacent points, in `s`
    `kmax`: maximum value of k, in `/pix`, typically `0.5`
    `oversamp`: oversampling ratio when solving numerial equation
    `flagDebugInfo`: whether print debug info

    # return
    kspace trajectory: [[kx1, ky1], [kx2, ky2], ..., [kxn, kyn]], in `/pix`
    gradient list: [[gx1, gy1], [gx2, gy2], ..., [gxn, gyn]], in `Hz/pix`
    '''
    sovQDE = lambda a, b, c: (-b+sqrt(max(b**2-4*a*c, 0)))/(2*a)
    lstTht = empty([0], dtype=float64)
    lstRho = empty([0], dtype=float64)

    d0ThtTime = 0
    d1ThtTime = 0
    d2ThtTime = 0
    d0RhoTht = getD0RhoTht(d0ThtTime) # ; assert(d0RhoTht == 0)
    d1RhoTht = getD1RhoTht(d0ThtTime)
    d2RhoTht = getD2RhoTht(d0ThtTime)

    lstTht = append(lstTht, d0ThtTime)
    lstRho = append(lstRho, d0RhoTht)
    idxPt = 1

    while d0RhoTht < kmax:
        sr = srlim*(1 - exp(-idxPt/oversamp))
        
        a = d0RhoTht**2 + d1RhoTht**2
        b = 2*d0RhoTht*d1RhoTht*d1ThtTime**2 + 2*d1RhoTht*d2RhoTht*d1ThtTime**2
        c = d0RhoTht**2*d1ThtTime**4 - 2*d0RhoTht*d2RhoTht*d1ThtTime**4 + 4*d1RhoTht**2*d1ThtTime**4 + d2RhoTht**2*d1ThtTime**4 - sr**2

        d2ThtTime = sovQDE(a, b, c)
        d1ThtTime += d2ThtTime*(dt/oversamp)
        if d1ThtTime*dt*d0RhoTht > glim*dt: d1ThtTime = min(d1ThtTime, glim/d0RhoTht) # dk >= d1ThtTime*dt*d0RhoTht
        d0ThtTime += d1ThtTime*(dt/oversamp)
        d0RhoTht = getD0RhoTht(d0ThtTime)
        d1RhoTht = getD1RhoTht(d0ThtTime)
        d2RhoTht = getD2RhoTht(d0ThtTime)

        lstTht = append(lstTht, d0ThtTime)
        lstRho = append(lstRho, d0RhoTht)
        idxPt += 1

        if flagDebugInfo and lstRho.size%1000 == 0: print(f"rho = {d0RhoTht:.2f}/{kmax:.2f}")

    lstRho = lstRho[::oversamp]
    lstTht = lstTht[::oversamp]
    lstTraj = array([
        lstRho*cos(lstTht),
        lstRho*sin(lstTht)]).T
    lstGrad = array([
        (lstTraj[1:,0] - lstTraj[:-1,0])/dt,
        (lstTraj[1:,1] - lstTraj[:-1,1])/dt]).T
    lstGrad = concatenate([array([[0, 0]]), lstGrad], axis=0)

    return lstTraj, lstGrad

def genRadial(lstTht:ndarray, lstRho:ndarray) -> ndarray:
    """
    # description:
    generate radial sampling trajectory

    # parameter:
    `lstTht`: list of theta of spokes, in `rad`
    `lstRho`: list of rho of spokes, in `/pix`

    # return:
    kspace trajectory: [[kx1,ky1], [kx2,ky2], ..., [kxn,kyn]], in `/pix`
    """
    # shape check
    assert(size(lstTht.shape) == 1)
    assert(size(lstRho.shape) == 1)
    
    # generate kspace trajectory
    lstKx = [lstRho*cos(lstTht[idxTht]) for idxTht in range(lstTht.size)]
    lstKy = [lstRho*sin(lstTht[idxTht]) for idxTht in range(lstTht.size)]

    return array([lstKx, lstKy]).transpose([1,2,0])

def genCart(numPt:int|float, max:int|float=0.5, numDim:int=2) -> ndarray:
    """
    # description:
    generate Cartesian sampling trajectory

    # parameter:
    `numPt`: number of point in one dimension
    `max`: maximum coordinate value

    # return:
    trajectory: [[kx1,ky1], [kx2,ky2], ..., [kxn,kyn]]
    """
    tupLstK = meshgrid(
        *[linspace(-max, max, numPt, endpoint=False) for _ in range(numDim)],
        indexing="ij")
    return array([lstK.flatten() for lstK in tupLstK]).T

def genSpiral3DTypeA(sr:int|float, numPix:int|float, tht0:int|float, phi0:int|float, uTht:int|float, uPhi:int|float, dt:int|float=10e-6, kmax:int|float=0.5) -> ndarray:
    """
    # description:
    generate Spiral3D-TypeA trajectory

    # parmaeter:
    `sr`: desired slewrate
    `numPix`: matrix size of acquired image
    `tht0`: initial phase of theata
    `phi0`: initial phase of phi
    `uTht`: undersamp ratio of theta
    `uPhi`: undersamp ratio of phi
    `dt`: temporal resolution of trajectory
    `kmax`: maximum of k, typically 0.5

    # return:
    trajectory: [[kx1,ky1,kz1], [kx2,ky2,kz2], ..., [kxn,kyn,kzn]]
    """
    return ext.GenSpiral3D(sr, numPix, tht0, phi0, uTht, uPhi, dt, kmax)