from numpy import *
from matplotlib.pyplot import *
from mrtrajgen import *

sizPix = 1e-3 # m
sizIm = 16

# generate trajectory of spiral
getDeltaK = lambda rho, tht: 1/sizIm
getDrhoDtht = lambda rho, tht: 0.5/((sizIm/2)*(2*pi))
trjSpiral = genSpiral(getDeltaK, getDrhoDtht) # derive trajectory

# calculate gradient list and slew rate list
lstGrad = tranTraj2Grad_MinRamp(trjSpiral, 10e-6)
lstSlewRate = getSlewRate(lstGrad, 10e-6)

# show trajectories
figure()

subplot(221)
plot(trjSpiral[:,0], trjSpiral[:,1], marker='.')
axis("equal"); title("Spiral")

subplot(223)
plot(lstGrad*(1/sizPix), marker=".")
title("Gradient")

subplot(224)
plot(lstSlewRate*(1/sizPix), marker=".")
title("Slew Rate")

subplots_adjust(hspace=0.4, wspace=0.2)

show()