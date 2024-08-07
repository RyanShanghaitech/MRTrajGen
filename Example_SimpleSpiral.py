from numpy import *
from matplotlib.pyplot import *
from mrtrajgen import *

sizPix = 1e-3 # m
sizImg = 16
turbo = 8

# generate trajectory of spiral
trjSpiral,_,_ = genSpiral(
    lambda rho, tht: 1/sizImg,
    lambda rho, tht: turbo*0.5/((sizImg/2)*(2*pi)),
    linspace(0, 2*pi, turbo)) # derive trajectory

# calculate gradient list and slew rate list
lstGrad = tranTraj2Grad_MinSR(trjSpiral[0,:,:], 10e-6)
lstSlewRate = getSlewRate(lstGrad, 10e-6)

# show trajectories
figure()

subplot(221)
plot(trjSpiral.reshape([-1,2])[:,0], trjSpiral.reshape([-1,2])[:,1], marker='.')
axis("equal"); title("Spiral")
subplot(223)
plot(lstGrad*(1/sizPix), marker=".")
title("Gradient")
subplot(224)
plot(lstSlewRate*(1/sizPix), marker=".")
title("Slew Rate")
subplots_adjust(hspace=0.4, wspace=0.2)

show()