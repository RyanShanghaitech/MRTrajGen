from numpy import *
from matplotlib.pyplot import *
from mrtrajgen import *

# parameter definition
gamma = 42.58e6
sizIm = 128
sizPix = 1e-3
dt = 10e-6
sr = 640 # slew rate in T/m/s
turbo = 1
rhoMax = 0.5
srLimit = getSlewRateCircle(1/sizIm, dt, rhoMax) # use slew rate at the boundary of kspace as slew rate limit

print(f"srLimit={srLimit*(1/1e-3)}")
rho1 = gamma*(srLimit*dt)*dt/2 # rho of first point
print(f"rho1={rho1}")

# generate trajectory of spiral
getDeltaK = lambda rho, tht: rho1 if rho == 0 else min(sqrt(srLimit*gamma*rho*(dt)**2), 1/sizIm) # function of sampling interval, with respect to rho and theta, used to make slew rate constant
getDrhoDtht = lambda rho, tht: turbo*0.5/(sizIm*pi) # function of rho/theta, with respect to rho and theta
trjSpiral = genSpiral(getDeltaK, getDrhoDtht, 0, rhoMax) # derive trajectory

print(f"pts of spiral: {turbo*trjSpiral.shape[0]}")
print(f"pts of cartes: {sizIm*sizIm}")

# show trajectories
figure()
plot(trjSpiral[:,0], trjSpiral[:,1], marker='.')
axis("equal"); title("Spiral")

# calculate gradient of spiral
lstGrad_Ideal = tranTraj2Grad_Ideal(trjSpiral, dt)*(1/sizPix)
lstGrad_MinRamp = tranTraj2Grad_MinRamp(trjSpiral, dt)*(1/sizPix)
lstGrad_MaxRamp = tranTraj2Grad_MaxRamp(trjSpiral, dt, sr*(sizPix/1))*(1/sizPix)

# show gradient
figure()

subplot(221)
plot(lstGrad_Ideal[:,0], label="Gx", marker='.')
plot(lstGrad_Ideal[:,1], label="Gy", marker='.')
legend()
title("lstGrad_Ideal")

subplot(222)
plot(lstGrad_MinRamp[:,0], label="Gx", marker='.')
plot(lstGrad_MinRamp[:,1], label="Gy", marker='.')
legend()
title("lstGrad_MinRamp")

subplot(223)
plot(lstGrad_MaxRamp[:,0], label="Gx", marker='.')
plot(lstGrad_MaxRamp[:,1], label="Gy", marker='.')
legend()
title("lstGrad_MaxRamp")

subplot(224)
plot(lstGrad_MaxRamp[:,0] - lstGrad_Ideal[:,0], label="Gx", marker='.')
plot(lstGrad_MaxRamp[:,1] - lstGrad_Ideal[:,1], label="Gy", marker='.')
legend()
title("diff(Ideal, MaxRamp)")

# show slew rate
lstSlewRateAbs_MinRamp = getSlewRate(lstGrad_MinRamp*(sizPix/1), dt)*(1/sizPix)

figure()
subplot(221)
plot(trjSpiral[:,0], trjSpiral[:,1], marker='.')
axis("equal"); title("Spiral")
subplot(222)
plot(trjSpiral[:,0], trjSpiral[:,1], marker='.')
axis("equal"); xlim([-0.03, 0.03]); ylim([-0.03, 0.03])
subplot(223)
plot(lstGrad_MinRamp, marker=".")
title("Gradient")
subplot(224)
plot(lstSlewRateAbs_MinRamp, marker=".")
title("Slew Rate")
subplots_adjust(hspace=0.4, wspace=0.2)

figure()
plot(lstSlewRateAbs_MinRamp, marker=".")
title("Slew Rate")

show()