from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import *

import mrtrjgen

# parameters
genSpiral3D = mrtrjgen.genSpiral3DTypeA # proposed Spiral-3D
# genSpiral3D = mrtrjgen.genSpiral3DTypeB # elementary approximation of Seiffert-Spiral

sr = 100 # desired slew rate
fov = 0.5
nPix = 256
nTht = 64
nPhi = 64
dtGrad = 10e-6
dtADC = 2e-6
gamma = 42.5756e6 # UIH

# calculate trajectory
lstArrK = []
lstArrG = []
lstArrSR = []
scale = 1/gamma*nPix/fov
for iPhi in range(nPhi):
    for iTht in range(nTht):
        tht0 = (2*pi)*(iTht/nTht) #+ (2*pi/nTht)*(iPhi/nPhi)
        phi0 = (2*pi)*(iPhi/nPhi)
        arrK, arrG = genSpiral3D(nPix, nTht, nPhi, tht0, phi0, 0.5, sr*gamma*fov/nPix, dtGrad)
        # arrK, _ = mrtrjgen.intpTraj(arrG, dtGrad, dtADC)
        arrSR = (arrG[1:,:] - arrG[:-1,:])/dtGrad
        lstArrK.append(arrK)
        lstArrG.append(arrG*scale)
        lstArrSR.append(asarray(sqrt((arrSR**2).sum(axis=-1)))*scale)

# plot
markersize = 2
linewidth = 1
fig = figure(figsize=(14,8), dpi=150)
gs = GridSpec(4,7,fig)
axK = fig.add_subplot(gs[0:4,0:4], projection="3d")
axGx = fig.add_subplot(gs[0:1,4:7])
axGy = fig.add_subplot(gs[1:2,4:7])
axGz = fig.add_subplot(gs[2:3,4:7])
axSR = fig.add_subplot(gs[3:4,4:7])

subplots_adjust(left=0.0, right=0.9, top=0.9, bottom=0.1, hspace=0.4, wspace=0.4)

for iTR in range(len(lstArrK)):
    axK.clear()
    axGx.clear()
    axGy.clear()
    axGz.clear()
    axSR.clear()

    axK.plot(lstArrK[iTR][:,0], lstArrK[iTR][:,1], lstArrK[iTR][:,2], ".-", markersize=markersize, linewidth=linewidth)
    axGx.plot(lstArrG[iTR][:,0], ".-", markersize=markersize, linewidth=linewidth)
    axGy.plot(lstArrG[iTR][:,1], ".-", markersize=markersize, linewidth=linewidth)
    axGz.plot(lstArrG[iTR][:,2], ".-", markersize=markersize, linewidth=linewidth)
    axSR.plot(lstArrSR[iTR][:], ".-", markersize=markersize, linewidth=linewidth)

    axK.set_xlim([-0.5,0.5])
    axK.set_ylim([-0.5,0.5])
    axK.set_zlim([-0.5,0.5])
    axGx.set_ylim([-40e-3,40e-3])
    axGy.set_ylim([-40e-3,40e-3])
    axGz.set_ylim([-40e-3,40e-3])
    axSR.set_ylim([0,2*sr])

    axK.set_title("K (/pix)")
    axGx.set_title("Gx (T/m)")
    axGy.set_title("Gy (T/m)")
    axGz.set_title("Gz (T/m)")
    axSR.set_title("Slew Rate (T/m/s)")

    show(block=0)
    pause(1e-3)

# figure()
# for iTR in range(len(lstArrK)):
#     arrRho = asarray(lstArrK[iTR])
#     arrRho = sqrt(sum(arrRho**2, axis=-1))
#     plot(arrRho, ".-")