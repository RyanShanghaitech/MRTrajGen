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
lstArrKxy = []
lstArrGxy = []
lstArrSR = []
scale = 1/gamma*nPix/fov
for idxTht in range(u):
    _, arrGxy = mrtrjgen.genSpiral2D(0.5/(2*pi)/8, 1e-3, 2*pi*idxTht/u, 0.5, sr*gamma*fov/nPix, dtGrad, 1e2)
    arrKxy, _ = mrtrjgen.intpTraj(arrGxy, dtGrad, dtADC)
    arrSR = (arrGxy[1:,:] - arrGxy[:-1,:])/dtGrad
    lstArrKxy.append(arrKxy)
    lstArrGxy.append(arrGxy*scale)
    lstArrSR.append(asarray(sqrt((arrSR**2).sum(axis=-1)))*scale)
    
# plot
fig = figure()
ax = fig.add_subplot(111)
for idxTR in range(u):
    ax.plot(lstArrKxy[idxTR][:,0], lstArrKxy[idxTR][:,1], ".-")
ax.axis("equal")
draw()

markersize = 2
linewidth = 1
fig = figure(figsize=(12,6), dpi=150)
gs = GridSpec(3,6,fig)
axKxy = fig.add_subplot(gs[0:3,0:3])
axGx = fig.add_subplot(gs[0:1,3:6])
axGy = fig.add_subplot(gs[1:2,3:6])
axS = fig.add_subplot(gs[2:3,3:6])

subplots_adjust(left=0.0, right=0.9, top=0.9, bottom=0.1, hspace=0.4, wspace=0.4)

while 1:
    for idxTR in range(u):
        axKxy.clear()
        axGx.clear()
        axGy.clear()
        axS.clear()
        axKxy.plot(lstArrKxy[idxTR][:,0], lstArrKxy[idxTR][:,1], ".-", markersize=markersize, linewidth=linewidth)
        axGx.plot(lstArrGxy[idxTR][:,0], ".-", markersize=markersize, linewidth=linewidth)
        axGy.plot(lstArrGxy[idxTR][:,1], ".-", markersize=markersize, linewidth=linewidth)
        axS.plot(lstArrSR[idxTR][:], ".-", markersize=markersize, linewidth=linewidth)
        axKxy.set_xlim([-0.5,0.5])
        axKxy.set_ylim([-0.5,0.5])
        axGx.set_ylim([-50e-3,50e-3])
        axGy.set_ylim([-50e-3,50e-3])
        axS.set_ylim([0,2*sr])
        axKxy.set_title("Kxy (/pix)")
        axGx.set_title("Gx (T/m)")
        axGy.set_title("Gy (T/m)")
        axS.set_title("Slew Rate (T/m/s)")
        show(block=0)
        pause(1e-6)
