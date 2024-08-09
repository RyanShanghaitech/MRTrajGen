from numpy import *
from matplotlib.pyplot import *
from mrtrajgen import *

sizPix = 1e-3 # m
sizImg = 16
turbo = 1

# generate trajectory of spiral
trjSpiral = genSpiral_DeltaK(
    lambda rho, tht: 1/sizImg,
    lambda rho, tht: turbo*0.5/((sizImg/2)*(2*pi))) # derive trajectory

# calculate gradient list and slew rate list
lstGrad = tranTraj2Grad_MinSR(trjSpiral, 10e-6)
lstSlewRate = getSlewRate(lstGrad, 10e-6)

# show trajectories
trjSpiral *= 1/sizPix
lstGrad *= 1/sizPix
lstSlewRate *= 1/sizPix
figure()
plot(trjSpiral[:,0], trjSpiral[:,1], marker='.')
axis("equal"); title("Spiral (/m)")
figure()
plot(lstGrad, marker=".")
title("Gradient (T/m)")
figure()
plot(lstSlewRate, marker=".")
title("Slew Rate (T/m/s)")
subplots_adjust(hspace=0.4, wspace=0.2)

show()