from numpy import *
from matplotlib.pyplot import *
from scipy.signal import fftconvolve, convolve

def tranTraj2Grad_Ideal(lstTraj:ndarray, dt:float|int) -> ndarray:
    """
    # description:
    transform trajectory in k-space to gradient,
    slew rate will be considered infinite

    # parameter:
    `lstTraj`: list of trajectory in k-space, `shape[0]` is the number of points, `shape[1]` is the dimension of k-space, in `/pix`, maximum should be `0.5`
    `dt`: time between two points in trajectory, in `s`

    # return:
    list of gradient in k-space, in `Hz/pix`

    # note:
    multiply the returned value by `(1 pix)/(length per pix in m)` to get the actual gradient in `m`.
    """

    lstTraj = concatenate([zeros([1, lstTraj.shape[1]]), lstTraj])
    lstDeltaK = (lstTraj[1:,:] - lstTraj[:-1,:]) # /pix
    lstGrad = lstDeltaK/dt # Hz/pix

    return lstGrad

def tranTraj2Grad_MinSR(lstTraj:ndarray, dt:float|int) -> ndarray:
    """
    # description:
    transform trajectory in k-space to gradient,
    ramp time will be considerd equal to `dt`, and slew rate will be assumed equal to `(Gnxt - Gpre)/(dt)`, which is the minimum slew rate to reach desired gradient amplitude in `dt`. 

    # parameter:
    `lstTraj`: list of trajectory in k-space, `shape[0]` is the number of points, `shape[1]` is the dimension of k-space, in `/pix`, maximum should be `0.5`
    `dt`: time between two points in trajectory, in `s`

    # return:
    list of gradient in k-space, in `Hz/pix`

    # note:
    multiply the returned value by `(1 pix)/(length per pix in m)` to get the actual gradient in `m`.
    """
    
    lstGrad_Ideal = tranTraj2Grad_Ideal(lstTraj, dt)
    lstGrad = zeros(lstGrad_Ideal.shape, dtype=float64)
    for idxGrad in range(lstGrad_Ideal.shape[0]):
        if idxGrad == 0:
            lstGrad[idxGrad,:] = 2*lstGrad_Ideal[idxGrad,:]
        else:
            lstGrad[idxGrad,:] = 2*lstGrad_Ideal[idxGrad,:] - lstGrad[idxGrad-1,:]

    return lstGrad

def tranTraj2Grad_MaxSR(lstTraj:ndarray, dt:float|int, sr:float|int) -> ndarray:
    """
    # description:
    transform trajectory in k-space to gradient,
    slew rate will be considered constant

    # parameter:
    `lstTraj`: list of trajectory in k-space, `shape[0]` is the number of points, `shape[1]` is the dimension of k-space, in `/pix`, maximum should be `0.5`
    `dt`: time between two points in trajectory, in `s`
    `sr`: slew rate, in `Hz/pix/s`

    # return:
    list of gradient in k-space, in `Hz/pix`

    # note:
    multiply the returned value by `(1 pix)/(length per pix in m)` to get the actual gradient in `m`.
    """
    
    lstGrad_Ideal = tranTraj2Grad_Ideal(lstTraj, dt)
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

def tranGrad2Slewrate(lstGrad:ndarray, dt:int|float) -> ndarray:
    """
    # description:
    get the slew rate of a gradient list

    # parameters:
    `lstGrad`: list represents the gradient, dim0 for points, dim1 for axis, in `Hz/pix`
    `dt`: time interval between two adjacent points, in `s`

    # return:
    list of magnitude of slew rate (Hz/pix/s)
    """
    lstSlewRate = (lstGrad[1:,:] - lstGrad[:-1,:])/dt # unit: Hz/pix/s
    
    return sqrt(lstSlewRate[:,0]**2 + lstSlewRate[:,1]**2)

def tranGrad2Traj_MinSR(lstGrad:ndarray, dt:int|float) -> ndarray:
    """
    # description:
    get the trajectory of a gradient list

    # parameters:
    `lstGrad`: list represents the gradient, dim0 for points, dim1 for axis, in `Hz/pix`
    `dt`: time interval between two adjacent points, in `s`

    # return:
    list of trajectory in k-space (/pix)
    """
    lstTraj = array([[0, 0]], dtype=float64)
    for idxGrad in range(1, lstGrad.shape[0]):
        dkx = (lstGrad[idxGrad, 0] + lstGrad[idxGrad-1, 0])*dt/2
        dky = (lstGrad[idxGrad, 1] + lstGrad[idxGrad-1, 1])*dt/2
        lstTraj = append(lstTraj, array([[lstTraj[-1,0]+dkx, lstTraj[-1,1]+dky]]), axis=0)
    return lstTraj

def getSlewRate_Circle(dk:float|int, dt:float|int, rho:float|int) -> ndarray:
    """
    # description:
    get the maximum slew rate for a given `dk`, `dt`, `rho` in spiral trajectory

    # parameter:
    `dk`: distance between two points in k-space, in `/pix`
    `dt`: time between two points in trajectory, in `s`
    `rho`: distance between `k` and origin, range from `0` to `0.5`, in `/pix`

    # return:
    maximum slew rate, in `Hz/pix/s`

    # note:
    multiply the returned value by `gamma*(1 pix)/(length per pix in m)` to get the actual slew rate, in `T/m/s`.
    """

    return (dk**2)/(rho*dt**2) # derivation of `d2G/dt2`

def getSlewRate_Pix(k1:ndarray, k0:ndarray, grad0:ndarray, dt:int|float) -> ndarray:
    '''
    # description:
    get the slew rate of two adjacent points in k-space

    # parameter:
    `k1`: k-space point 1, in `/pix`
    `k0`: k-space point 0, in `/pix`
    `grad0`: gradient at point 0, in `Hz/pix`
    `dt`: time interval between two adjacent points, in `s`

    # return:
    slew rate, in `Hz/pix/s`
    '''
    assert(k1.shape == (2,) and k0.shape == (2,) and grad0.shape == (2,))
    gradMean = (k1 - k0)/dt
    return (gradMean - grad0)/(dt/2)

def copyTraj(traj:ndarray, numCopy:int, intv:int|float=None):
    """
    # description:
    copy trajectory

    # parameter:
    `traj`: trajectory to copy
    `numCopy`: number of copies
    `intv`: interval between two copies, in `rad`

    # return:
    copied trajectory
    """
    assert(traj.ndim == 2) # there should be only one trajectory
    assert(traj.shape[1] == 2) # only support 2D trajectory
    if intv is None:
        intv = (2*pi)/numCopy
    trjNew = zeros([numCopy, traj.shape[0], traj.shape[1]], dtype=float64)
    for idxCopy in range(numCopy):
        trjNew[idxCopy,:,0] = cos(intv*idxCopy)*traj[:,0] - sin(intv*idxCopy)*traj[:,1]
        trjNew[idxCopy,:,1] = sin(intv*idxCopy)*traj[:,0] + cos(intv*idxCopy)*traj[:,1]
    return trjNew

def intpTraj(arrG:ndarray, dtGrad:int|float, dtADC:int|float) -> tuple[ndarray,ndarray]:
    """
    # description:
    interpolate gradient waveform and derive trajectory

    # parameter
    `arrG`: array of gradient waveform
    `dtGrad`, `dtADC`: temporal resolution of gradient system and ADC

    # return:
    interpolated trajectory and gradient
    """
    nGrad, nDim = arrG.shape
    nADC = int(dtGrad/dtADC)*(nGrad - 1)
    arrG_Resamp = zeros([nADC,nDim], dtype=float64)
    for iDim in range(nDim):
        arrG_Resamp[:,iDim] = interp(dtADC*arange(nADC) + dtADC/2, dtGrad*arange(nGrad), arrG[:,iDim])
    arrDk = zeros_like(arrG_Resamp)
    arrDk[0,:] = (0 + arrG_Resamp[0,:])*dtADC/2
    arrDk[1:] = (arrG_Resamp[:-1] + arrG_Resamp[1:])*dtADC/2
    arrK = cumsum(arrDk,axis=0)
    return arrK, arrG_Resamp

def delayGrad(arrG:ndarray, tau:int|float) -> ndarray:
    """
    # description:
    delay the input gradient waveform by time constant tau

    # parameter
    `arrG`: array of single gradient waveform
    `tau`: time constant in RL circuit transfer function

    # return:
    delayed gradient waveform
    """
    assert arrG.ndim == 2, "only single gradient waveform is supported."
    if tau == 0: return arrG.copy() # avoid divided-by-0 later
    nPt, nAx = arrG.shape

    # perform oversample to get better impluse response profile
    ov = clip(10/tau, 1, 1e3).astype(int64) # the smaller the ov, the bigger oversampling is needed
    arrG_ov = zeros([nPt*ov,nAx], dtype=arrG.dtype)
    for iAx in range(nAx):
        arrG_ov[:,iAx] = interp(linspace(0,nPt,nPt*ov,0), linspace(0,nPt,nPt,0), arrG[:,iAx]) # oversample
    nPt *= ov
    tau *= ov

    # derive impluse response of RL circuit
    arrG_ov = concatenate([arrG_ov, zeros_like(arrG_ov)], axis=0)
    arrT = linspace(0,2*nPt,2*nPt,0) + 0.5
    arrImpResRL = (1/tau)*exp(-arrT/tau)
    if abs(arrImpResRL.sum() - 1) > 1e-2: raise ValueError(f"arrImpResRL.sum() = {arrImpResRL.sum():.2f} (supposed to be 1) (tau too small or too large)")
    
    # perform convolution between input waveform and impulse response
    for iAx in range(nAx):
        # arrG_ov[:,iAx] = fftconvolve(arrG_ov[:,iAx], arrImpResRL, mode="same")
        arrG_ov[:,iAx] = fft.ifft(fft.fft(arrG_ov[:,iAx])*fft.fft(arrImpResRL)).real

    # compensate the aliased signal produced by exp(-t/tau)
    arrG_ov = arrG_ov[:nPt,:] - arrG_ov[nPt:,:]*(arrG_ov[0,:]-arrG[0,:])/(arrG_ov[nPt,:])

    # de-oversample
    arrG = arrG_ov[::ov]

    return arrG