from numpy import *
from matplotlib.pyplot import *

def tranTraj2Grad_Ideal(lstTraj:ndarray, dt:float|int, gamma:float|int=42.58e6) -> ndarray:
    """
    # description:
    transform trajectory in k-space to gradient,
    slew rate will be considered infinite

    # parameter:
    `lstTraj`: list of trajectory in k-space, `shape[0]` is the number of points, `shape[1]` is the dimension of k-space, in `/pix`, maximum should be `0.5`
    `dt`: time between two points in trajectory, in `s`
    `gamma`: gyromagnetic ratio, in `Hz/T`

    # return:
    list of gradient in k-space, in `T/pix`

    # note:
    multiply the returned value by `(1 pix)/(length per pix in m)` to get the actual gradient in `m`.
    """

    lstDeltaK = (lstTraj[1:,:] - lstTraj[:-1,:]) # /pix
    lstGrad = lstDeltaK/(gamma*dt) # T/pix

    return lstGrad

def tranTraj2Grad_MinSR(lstTraj:ndarray, dt:float|int, gamma:float|int=42.58e6) -> ndarray:
    """
    # description:
    transform trajectory in k-space to gradient,
    ramp time will be considerd equal to `dt`, and slew rate will be assumed equal to `(Gnxt - Gpre)/(dt)`, which is the minimum slew rate to reach desired gradient amplitude in `dt`. 

    # parameter:
    `lstTraj`: list of trajectory in k-space, `shape[0]` is the number of points, `shape[1]` is the dimension of k-space, in `/pix`, maximum should be `0.5`
    `dt`: time between two points in trajectory, in `s`
    `gamma`: gyromagnetic ratio, in `Hz/T`

    # return:
    list of gradient in k-space, in `T/pix`

    # note:
    multiply the returned value by `(1 pix)/(length per pix in m)` to get the actual gradient in `m`.
    """
    
    lstGrad_Ideal = tranTraj2Grad_Ideal(lstTraj, dt, gamma)
    lstGrad = zeros(lstGrad_Ideal.shape, dtype=float64)
    for idxGrad in range(lstGrad_Ideal.shape[0]):
        if idxGrad == 0:
            lstGrad[idxGrad,:] = 2*lstGrad_Ideal[idxGrad,:]
        else:
            lstGrad[idxGrad,:] = 2*lstGrad_Ideal[idxGrad,:] - lstGrad_Ideal[idxGrad-1,:]

    return lstGrad

def tranTraj2Grad_MaxSR(lstTraj:ndarray, dt:float|int, sr:float|int, gamma:float|int=42.58e6) -> ndarray:
    """
    # description:
    transform trajectory in k-space to gradient,
    slew rate will be considered constant

    # parameter:
    `lstTraj`: list of trajectory in k-space, `shape[0]` is the number of points, `shape[1]` is the dimension of k-space, in `/pix`, maximum should be `0.5`
    `dt`: time between two points in trajectory, in `s`
    `sr`: slew rate, in `T/pix/s`
    `gamma`: gyromagnetic ratio, in `Hz/T`

    # return:
    list of gradient in k-space, in `T/pix`

    # note:
    multiply the returned value by `(1 pix)/(length per pix in m)` to get the actual gradient in `m`.
    """
    
    lstGrad_Ideal = tranTraj2Grad_Ideal(lstTraj, dt, gamma)
    lstGrad = zeros(lstGrad_Ideal.shape, dtype=float64)
    for idxGrad in range(lstGrad_Ideal.shape[0]):
        prodGradTime = lstGrad_Ideal[idxGrad,:]*dt
        if idxGrad == 0:
            gradPrev = zeros(lstGrad_Ideal.shape[1], dtype=float64)
        else:
            gradPrev = lstGrad[idxGrad-1,:]
        s = sign(lstGrad_Ideal[idxGrad,:] - gradPrev)
        s[s == 0] = 1
        lstGrad[idxGrad,:] = gradPrev - s*sr*dt + s*sqrt((s*sr*dt - gradPrev)**2 - gradPrev**2 + prodGradTime*2*sr*s)

    return lstGrad

def getSlewRateCircle(dk:float|int, dt:float|int, rho:float|int, gamma:float|int=42.58e6) -> ndarray:
    """
    # description:
    get the maximum slew rate for a given `dk`, `dt`, `rho` in spiral trajectory

    # parameter:
    `dk`: distance between two points in k-space, in `/pix`
    `dt`: time between two points in trajectory, in `s`
    `rho`: distance between `k` and origin, range from `0` to `0.5`, in `/pix`
    `gamma`: gyromagnetic ratio, in `Hz/T`

    # return:
    maximum slew rate, in `T/pix/s`

    # note:
    multiply the returned value by `(1 pix)/(length per pix in m)` to get the actual slew rate in `T/m/s`.
    """

    return (dk**2)/(gamma*rho*dt**2) # derivation of `d2G/dt2`

def getSlewRate(lstGrad:ndarray, dt:int|float) -> ndarray:
    """
    # description:
    get the slew rate of a gradient list

    # parameters:
    `lstGrad`: list represents the gradient, dim0 for points, dim1 for axis
    `dt`: time interval between two adjacent points

    # return:
    list of magnitude of slew rate
    """
    lstSlewRate = (lstGrad[:,:] - concatenate([[[0, 0]], lstGrad[:-1,:]]))/dt # unit: T/pix/s
    
    return sqrt(lstSlewRate[:,0]**2 + lstSlewRate[:,1]**2)