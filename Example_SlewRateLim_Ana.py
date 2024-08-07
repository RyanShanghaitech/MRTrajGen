from numpy import *
from matplotlib.pyplot import *
from mrtrajgen import *

# parameter definition
gamma = 42.58e6
sizImg = 128
sizPix = 1e-3
dt = 10e-6
sr = 640 # slew rate in T/m/s
turbo = 1
rhoMax = 0.5
srLim = getSlewRateCircle(1/sizImg, dt, rhoMax) # use slew rate at the boundary of kspace as slew rate limit
print(f"srLim={srLim}")

# generate trajectory of spiral
rho1 = gamma*(srLim*dt)*dt/2 # rho of first point
getDeltaK = lambda rho, tht: rho1 if rho == 0 else min(sqrt(srLim*gamma*rho*(dt)**2), 1/sizImg) # function of sampling interval, with respect to rho and theta, used to make slew rate constant
getDrhoDtht = lambda rho, tht: turbo*0.5/(sizImg*pi) # function of drho/dtheta, with respect to rho and theta
trjSpiral,_,_ = genSpiral(getDeltaK, getDrhoDtht, [0], rhoMax) # derive trajectory
trjSpiral = trjSpiral.reshape([-1, 2])

print(f"pts of spiral: {turbo*trjSpiral.shape[0]}")
print(f"pts of cartes: {sizImg*sizImg}")

# calculate gradient of spiral
lstGrad_Ideal = tranTraj2Grad_Ideal(trjSpiral, dt)*(1/sizPix)
lstGrad_MaxSR = tranTraj2Grad_MaxSR(trjSpiral, dt, sr*(sizPix/1))*(1/sizPix)
lstGrad_MinSR = tranTraj2Grad_MinSR(trjSpiral, dt)*(1/sizPix)

# show gradient
figure()

subplot(311)
plot(lstGrad_Ideal[:,0], label="Gx", marker='.')
plot(lstGrad_Ideal[:,1], label="Gy", marker='.')
legend()
title("lstGrad_Ideal")

subplot(312)
plot(lstGrad_MaxSR[:,0], label="Gx", marker='.')
plot(lstGrad_MaxSR[:,1], label="Gy", marker='.')
legend()
title("lstGrad_MaxSR")

subplot(313)
plot(lstGrad_MinSR[:,0], label="Gx", marker='.')
plot(lstGrad_MinSR[:,1], label="Gy", marker='.')
legend()
title("lstGrad_MinSR")

# show slew rate
lstSlewRateAbs_MinSR = getSlewRate(lstGrad_MinSR*(sizPix/1), dt)*(1/sizPix)

figure()
subplot(221)
plot(trjSpiral[:,0], trjSpiral[:,1], marker='.')
axis("equal"); title("Spiral")
subplot(222)
plot(trjSpiral[:,0], trjSpiral[:,1], marker='.')
axis("equal"); xlim([-0.03, 0.03]); ylim([-0.03, 0.03])
subplot(223)
plot(lstGrad_MinSR, marker=".")
title("Gradient")
subplot(224)
plot(lstSlewRateAbs_MinSR, marker=".")
title("Slew Rate")
subplots_adjust(hspace=0.4, wspace=0.2)

show()