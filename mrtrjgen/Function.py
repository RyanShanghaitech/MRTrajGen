from numpy import *
from typing import Callable
from .Utility import tranGrad2Traj_MinSR
from . import ext
import math

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

def genSpiral2D(numPix:int|float, u:int|float, tht0:int|float, kmax:int|float, sr:int|float, dt:int|float=10e-6, ov:int|float=1e2) -> ndarray:
    """
    # description:
    generate Spiral3D-TypeA trajectory

    # parmaeter:
    `numPix`: matrix size of acquired image
    `u`: undersamp ratio
    `tht0`: initial phase of theata
    `kmax`: maximum of k, typically 0.5
    `sr`: desired slewrate
    `dt`: temporal resolution of trajectory
    `ov`: oversampling factor when recur

    # return:
    trajectory: [[kx0,ky0], [kx1,ky1], ..., [kxn,kyn]]
    gradient: [[gx0,gy0], [gx1,gy1], ..., [gxn,gyn]]
    """
    return ext.GenSpiral2D(numPix, u, tht0, kmax, sr, dt, ov)

def genSpiral3DTypeA(numPix:int|float, uTht:int|float, uPhi:int|float, tht0:int|float, phi0:int|float, kmax:int|float, sr:int|float, dt:int|float=10e-6, ov:int|float=1) -> ndarray:
    """
    # description:
    generate Spiral3D-TypeA trajectory

    # parmaeter:
    `numPix`: matrix size of acquired image
    `uTht`: undersamp ratio of theta
    `uPhi`: undersamp ratio of phi
    `tht0`: initial phase of theata
    `phi0`: initial phase of phi
    `kmax`: maximum of k, typically 0.5
    `sr`: desired slewrate
    `dt`: temporal resolution of trajectory
    `ov`: oversampling factor when recur

    # return:
    trajectory: [[kx0,ky0,kz0], [kx1,ky1,kz1], ..., [kxn,kyn,kzn]]
    gradient: [[gx0,gy0,gz0], [gx1,gy1,gz1], ..., [gxn,gyn,gzn]]
    """
    return ext.GenSpiral3D_A(numPix, uTht, uPhi, tht0, phi0, kmax, sr, dt)

def genSpiral3DTypeB(numPix:int|float, uTht:int|float, uPhi:int|float, tht0:int|float, phi0:int|float, kmax:int|float, sr:int|float, dt:int|float=10e-6, ov:int|float=1) -> ndarray:
    """
    # description:
    generate Spiral3D-TypeB trajectory

    # parmaeter:
    `numPix`: matrix size of acquired image
    `uTht`: undersamp ratio of theta
    `uPhi`: undersamp ratio of phi
    `tht0`: initial phase of theata
    `phi0`: initial phase of phi
    `kmax`: maximum of k, typically 0.5
    `sr`: desired slewrate
    `dt`: temporal resolution of trajectory
    `ov`: oversampling factor when recur

    # return:
    trajectory: [[kx0,ky0,kz0], [kx1,ky1,kz1], ..., [kxn,kyn,kzn]]
    gradient: [[gx0,gy0,gz0], [gx1,gy1,gz1], ..., [gxn,gyn,gzn]]
    """
    return ext.GenSpiral3D_B(numPix, uTht, uPhi, tht0, phi0, kmax, sr, dt)
