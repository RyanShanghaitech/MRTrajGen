from numpy import *
from . import ext

def genRadial(arrTht:ndarray, arrRho:ndarray) -> ndarray:
    """
    # description:
    generate Radial sampling trajectory

    # parameter:
    `arrTht`: list of theta of spokes, in `rad`
    `arrRho`: list of rho of spokes, in `/pix`

    # return:
    kspace trajectory: [[kx1,ky1], [kx2,ky2], ..., [kxn,kyn]], in `/pix`
    """
    # shape check
    assert(size(arrTht.shape) == 1)
    assert(size(arrRho.shape) == 1)
    
    # generate kspace trajectory
    lstKx = [arrRho*cos(arrTht[idxTht]) for idxTht in range(arrTht.size)]
    lstKy = [arrRho*sin(arrTht[idxTht]) for idxTht in range(arrTht.size)]

    return array([lstKx, lstKy]).transpose([1,2,0])

def genCart(nPix:int|float, max:int|float=0.5, numDim:int=2) -> ndarray:
    """
    # description:
    generate Cartesian sampling trajectory

    # parameter:
    `nPix`: number of point in one dimension
    `max`: maximum coordinate value

    # return:
    trajectory: [[kx1,ky1], [kx2,ky2], ..., [kxn,kyn]]
    """
    tupLstK = meshgrid(
        *[linspace(-max, max, nPix, endpoint=False) for _ in range(numDim)],
        indexing="ij")
    return array([lstK.flatten() for lstK in tupLstK]).T

def genSpiral2D(nPix:int|float, nSp:int|float, ovVds:int|float, tht0:int|float, kmax:int|float, sr:int|float, dt:int|float=10e-6, ov:int|float=1e2) -> tuple[ndarray, ndarray]:
    """
    # description:
    generate Spiral2D trajectory

    # parmaeter:
    `nPix`: number of pixel
    `nSp`: number of spiral
    `ovVds`: oversampling factor of Variable Density Spiral, 1 for Constant Density Spiral
    `tht0`: initial phase of theata
    `kmax`: maximum of k, typically 0.5
    `sr`: desired slewrate
    `dt`: temporal resolution of trajectory
    `ov`: oversampling factor when recur

    # return:
    trajectory: [[kx0,ky0], [kx1,ky1], ..., [kxn,kyn]]
    gradient: [[gx0,gy0], [gx1,gy1], ..., [gxn,gyn]]
    """
    thtMax = (2*pi)*(nPix/2)/nSp
    dKtht1 = (kmax/thtMax)/ovVds
    dKtht2 = ((nSp/(2*pi*nPix))**2 - dKtht1**2)/2
    return ext.GenSpiral2D(dKtht1, dKtht2, tht0, kmax, sr, dt, ov)

def genSpiral3DTypeA(nPix:int|float, uTht:int|float, uPhi:int|float, tht0:int|float, phi0:int|float, kmax:int|float, sr:int|float, dt:int|float=10e-6, ov:int|float=1) -> tuple[ndarray, ndarray]:
    """
    # description:
    generate Spiral3D-TypeA trajectory

    # parmaeter:
    `nPix`: matrix size of acquired image
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
    return ext.GenSpiral3D_A(nPix, uTht, uPhi, tht0, phi0, kmax, sr, dt)

def genSpiral3DTypeB(nPix:int|float, uTht:int|float, uPhi:int|float, tht0:int|float, phi0:int|float, kmax:int|float, sr:int|float, dt:int|float=10e-6, ov:int|float=1) -> tuple[ndarray, ndarray]:
    """
    # description:
    generate Spiral3D-TypeB trajectory

    # parmaeter:
    `nPix`: matrix size of acquired image
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
    return ext.GenSpiral3D_B(nPix, uTht, uPhi, tht0, phi0, kmax, sr, dt)
