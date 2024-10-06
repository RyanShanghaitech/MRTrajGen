from numpy import *
from matplotlib.pyplot import *
from mrtrjgen import *

# parameters
numPix = 128
widFov = 0.4 # m
gamma = 42.58e6
dt = 10e-6
sr = 100*(widFov/numPix) # T/pix/s
flagVariableDensity = 0
if flagVariableDensity:
    turbo = 1 # variable density spiral
    evoRhoTht = 0.1
else:
    turbo = 48 # homogeneous spiral
    evoRhoTht = 0
kmax = 0.5 # /pix
quoRhoTht = kmax/(2*pi)/(numPix/(2*turbo))

# generate trajectory
lstTraj, lstGrad = genSpiral_Slewrate(
    lambda tht: quoRhoTht*(tht + evoRhoTht/2*tht**2),
    lambda tht: quoRhoTht*(1 + evoRhoTht*tht),
    lambda tht: quoRhoTht*evoRhoTht,
    sr, inf, dt, kmax, 10, True)

print(lstTraj[:10,:])
print(lstGrad[:10,:])

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
plot(lstSR*(numPix/widFov), marker=".")
ylim([0, 2*sr*(numPix/widFov)]); title("slew rate (T/m/s)")

subplot(133)
plot(lstSR*(numPix/widFov), marker=".")
title("slew rate (T/m/s)")

figure()
plot(lstGrad[:,0]*(numPix/widFov), marker=".")
plot(lstGrad[:,1]*(numPix/widFov), marker=".")
title("gradient (T/m)")

show()