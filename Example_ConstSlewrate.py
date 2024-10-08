from numpy import *
from matplotlib.pyplot import *
from mrtrjgen import *

# parameters
numPix = 128
widFov = 0.2 # m
gamma = 42.58e6
dt = 10e-6
sr = 100*gamma*(widFov/numPix) # Hz/pix/s
flagVariableDensity = 0
if flagVariableDensity:
    turbo = 1 # variable density spiral
    evoRhoTht = 1.1
else:
    turbo = 64 # homogeneous spiral
    evoRhoTht = 1
kmax = 0.5 # /pix
quoRhoTht = kmax/(2*pi)/(numPix/(2*turbo))

# generate trajectory
lstTraj, lstGrad = genSpiral_Slewrate(
    lambda tht: quoRhoTht*(tht + (evoRhoTht-1)/2*tht**2),
    lambda tht: quoRhoTht*(1 + (evoRhoTht-1)*tht),
    lambda tht: quoRhoTht*(evoRhoTht-1),
    sr, inf, dt, kmax, 10, True)

numCopy = turbo
lstTraj = copyTraj(lstTraj, numCopy)

# derive slew rate
lstSR = tranGrad2Slewrate(lstGrad, dt)

# plot
figure()

subplot(131)
for idxTrj in range(numCopy):
    plot(lstTraj[idxTrj,:,0], lstTraj[idxTrj,:,1], marker=".")
axis("equal"); title("kx-ky (/m)")

subplot(132)
plot(lstSR/gamma*(numPix/widFov), marker=".")
ylim([0, 2*sr/gamma*(numPix/widFov)]); title("slew rate (T/m/s)")

subplot(133)
plot(lstSR/gamma*(numPix/widFov), marker=".")
title("slew rate (T/m/s)")

figure()
plot(lstGrad[:,0]/gamma*(numPix/widFov), marker=".")
plot(lstGrad[:,1]/gamma*(numPix/widFov), marker=".")
title("gradient (T/m)")

show()