from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import *

import mrtrjgen

sr = 100 # desired slew rate
fov = 0.25
numPix = 128
u = 8
dt = 10e-6
gamma = 42.5756e6 # UIH

# calculate trajectory
lstArrKxy = []
lstArrGxy = []
lstArrSR = []
scale = 1/gamma*numPix/fov
for idxTht in range(u):
    arrKxy, arrGxy = mrtrjgen.genSpiral2D(numPix, u, 2*pi*idxTht/u, 0.5, sr*gamma*fov/numPix, dt, 1e2)
    arrSR = (arrGxy[1:,:] - arrGxy[:-1,:])/dt
    lstArrKxy.append(arrKxy)
    lstArrGxy.append(arrGxy*scale)
    lstArrSR.append(asarray(sqrt((arrSR**2).sum(axis=-1)))*scale)
    
# plot
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
        axGx.set_ylim([-30e-3,30e-3])
        axGy.set_ylim([-30e-3,30e-3])
        axS.set_ylim([0,2*sr])
        axKxy.set_title("Kxy (/pix)")
        axGx.set_title("Gx (T/m)")
        axGy.set_title("Gy (T/m)")
        axS.set_title("Slew Rate (T/m/s)")
        show(block=0)
        pause(1e-6)
