from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import *

import mrtrjgen

sr = 100 # desired slew rate
fov = 0.25
nPix = 256
u = 24
dtGrad = 10e-6 # temporal resolution of gradient coil
dtADC = 2.5e-6 # temporal resolution of ADC
gamma = 42.5756e6 # UIH

# calculate trajectory
lstArrK = []
lstArrG = []
lstArrSR = []
scale = 1/gamma*nPix/fov
for idxTht in range(u):
    _, arrG = mrtrjgen.genSpiral2D(0.5/(2*pi)/8, 1e-3, 2*pi*idxTht/u, 0.5, sr*gamma*fov/nPix, dtGrad, 1e2)
    arrK, _ = mrtrjgen.intpTraj(arrG, dtGrad, dtADC)
    arrSR = (arrG[1:,:] - arrG[:-1,:])/dtGrad
    lstArrK.append(arrK)
    lstArrG.append(arrG*scale)
    lstArrSR.append(asarray(sqrt((arrSR**2).sum(axis=-1)))*scale)
    
# plot
fig = figure()
ax = fig.add_subplot(111)
for idxTR in range(u):
    ax.plot(lstArrK[idxTR][:,0], lstArrK[idxTR][:,1], ".-")
ax.axis("equal")
draw()

markersize = 2
linewidth = 1
fig = figure(figsize=(12,6), dpi=150)
gs = GridSpec(3,6,fig)
axK = fig.add_subplot(gs[0:3,0:3])
axGx = fig.add_subplot(gs[0:1,3:6])
axGy = fig.add_subplot(gs[1:2,3:6])
axSR = fig.add_subplot(gs[2:3,3:6])

subplots_adjust(left=0.0, right=0.9, top=0.9, bottom=0.1, hspace=0.4, wspace=0.4)

while 1:
    for idxTR in range(u):
        axK.clear()
        axGx.clear()
        axGy.clear()
        axSR.clear()
        axK.plot(lstArrK[idxTR][:,0], lstArrK[idxTR][:,1], ".-", markersize=markersize, linewidth=linewidth)
        axGx.plot(lstArrG[idxTR][:,0], ".-", markersize=markersize, linewidth=linewidth)
        axGy.plot(lstArrG[idxTR][:,1], ".-", markersize=markersize, linewidth=linewidth)
        axSR.plot(lstArrSR[idxTR][:], ".-", markersize=markersize, linewidth=linewidth)
        axK.set_xlim([-0.5,0.5])
        axK.set_ylim([-0.5,0.5])
        axGx.set_ylim([-50e-3,50e-3])
        axGy.set_ylim([-50e-3,50e-3])
        axSR.set_ylim([0,2*sr])
        axK.set_title("Kxy (/pix)")
        axGx.set_title("Gx (T/m)")
        axGy.set_title("Gy (T/m)")
        axSR.set_title("Slew Rate (T/m/s)")
        show(block=0)
        pause(1e-6)
