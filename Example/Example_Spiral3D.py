from numpy import *
from matplotlib.pyplot import *
from matplotlib import gridspec
from mpl_toolkits.mplot3d import *

import mrtrjgen

sr = 100 # desired slew rate
fov = 0.25
numPix = 128
uTht = 32
uPhi = 32
dt = 10e-6
gamma = 42.5756e6 # UIH

lstArrKxyz = []
lstArrGxyz = []
lstArrS = []
scale = 1/gamma*numPix/fov
for idxTht in range(uTht):
    for idxPhi in range(uPhi):
        arrKxyz, arrGxyz = mrtrjgen.genSpiral3DTypeA(sr*gamma*fov/numPix, numPix, 2*pi*idxTht/uTht, 2*pi*idxPhi/uPhi, uTht, uPhi, dt, 0.5)
        arrS = (arrGxyz[1:,:] - arrGxyz[:-1,:])/dt
        lstArrKxyz.append(arrKxyz)
        lstArrGxyz.append(arrGxyz*scale)
        lstArrS.append(asarray(sqrt((arrS**2).sum(axis=-1)))*scale)

# plot
markersize = 2
linewidth = 1
fig = figure(figsize=(14,8), dpi=150)
gs = GridSpec(4,7,fig)
axKxyz = fig.add_subplot(gs[0:4,0:4], projection="3d")
axGx = fig.add_subplot(gs[0:1,4:7])
axGy = fig.add_subplot(gs[1:2,4:7])
axGz = fig.add_subplot(gs[2:3,4:7])
axS = fig.add_subplot(gs[3:4,4:7])

subplots_adjust(left=0.0, right=0.9, top=0.9, bottom=0.1, hspace=0.4, wspace=0.4)

for idxTR in range(uTht*uPhi):
    axKxyz.clear()
    axGx.clear()
    axGy.clear()
    axGz.clear()
    axS.clear()
    axKxyz.plot(lstArrKxyz[idxTR][:,0], lstArrKxyz[idxTR][:,1], lstArrKxyz[idxTR][:,2], ".-", markersize=markersize, linewidth=linewidth)
    axGx.plot(lstArrGxyz[idxTR][:,0], ".-", markersize=markersize, linewidth=linewidth)
    axGy.plot(lstArrGxyz[idxTR][:,1], ".-", markersize=markersize, linewidth=linewidth)
    axGz.plot(lstArrGxyz[idxTR][:,2], ".-", markersize=markersize, linewidth=linewidth)
    axS.plot(lstArrS[idxTR][:], ".-", markersize=markersize, linewidth=linewidth)
    axKxyz.set_xlim([-0.5,0.5])
    axKxyz.set_ylim([-0.5,0.5])
    axKxyz.set_zlim([-0.5,0.5])
    axGx.set_ylim([-30e-3,30e-3])
    axGy.set_ylim([-30e-3,30e-3])
    axGz.set_ylim([-30e-3,30e-3])
    axS.set_ylim([0,2*sr])
    axKxyz.set_title("Kxyz (/pix)")
    axGx.set_title("Gx (T/m)")
    axGy.set_title("Gy (T/m)")
    axGz.set_title("Gz (T/m)")
    axS.set_title("Slew Rate (T/m/s)")
    draw()
    pause(1e-7)
