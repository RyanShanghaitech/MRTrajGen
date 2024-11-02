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
    kspace trajectory: [[kx1, ky1], [kx2, ky2], ..., [kxn, kyn]], in `/pix`
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
    generate cartesian sampling trajectory

    # parameter:
    `numPt`: number of point in one dimension
    `max`: maximum coordinate value

    # return:
    trajectory: [[x1, y1], [x2, y2], ..., [xn, yn]]
    """
    tupLstK = meshgrid(
        *[linspace(-max, max, numPt, endpoint=False) for _ in range(numDim)],
        indexing="ij")
    return array([lstK.flatten() for lstK in tupLstK]).T

def genSpiral3DTypeA(uTht:float|int, uPhi:float|int, tht0:float|int, phi0:float|int, sr:float|int, numPix:float|int, dt:float|int):
    sovQDE = lambda a, b, c: (-b+sqrt(max(b**2-4*a*c, 0)))/(2*a)

    u_tht = uTht
    u_phi = uPhi
    Np = numPix
    s = sr

    rho = (0.01*s)*dt*dt/2 # since solver of d1phi doesn't allow phi=0, we have to calculate the initial point
    phi = (2*pi*Np/u_phi)*rho
    tht = u_phi/(2*u_tht)*(phi**2)
    
    tht_d1 = 0
        
    lstKx = [0, rho*sin(tht + tht0)*cos(phi + phi0)]
    lstKy = [0, rho*sin(tht + tht0)*sin(phi + phi0)]
    lstKz = [0, rho*cos(tht + tht0)]

    while rho <= 0.5:
        a = (1/32)*u_tht*(8*tht**3*u_phi*(2*tht*u_phi + u_tht*math.cos(2*tht0 + 2*tht) + 3*u_tht) + 24*tht**2*u_phi**2 + 4*tht**2*u_tht**2*math.sin(tht0 + tht)**2 + 6*tht*u_phi*u_tht*math.sin(tht0 + tht)**2 + u_phi**2)/(Np**2*pi**2*tht**3*u_phi)
        b = 0
        c = -s**2
    
        tht_d1 = sovQDE(a,b,c)
        tht += tht_d1*dt
        
        phi = sqrt(2*u_tht/u_phi)*sqrt(tht)
        rho = u_phi/(2*pi*Np)*phi

        kx = rho*sin(tht + tht0)*cos(phi + phi0)
        ky = rho*sin(tht + tht0)*sin(phi + phi0)
        kz = rho*cos(tht + tht0)

        lstKx.append(kx)
        lstKy.append(ky)
        lstKz.append(kz)
        
    arrKx = array(lstKx)
    arrKy = array(lstKy)
    arrKz = array(lstKz)

    return array([arrKx, arrKy, arrKz]).T

def genSpiral3DTypeA_Cpp(sr:int|float, numPix:int|float, tht0:int|float, phi0:int|float, uTht:int|float, uPhi:int|float, dt:int|float, kmax:int|float) -> ndarray:
    return ext.GenSpiral3D(sr, numPix, tht0, phi0, uTht, uPhi, dt, kmax)