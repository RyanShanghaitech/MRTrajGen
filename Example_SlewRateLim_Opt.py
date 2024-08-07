from mrtrajgen import *
from numpy import *
from matplotlib.pyplot import *
from scipy.spatial import voronoi_plot_2d
from scipy.interpolate import interp1d

sizImg = 128
turbo = 8
gamma = 42.58e6
sizPix = (300/256)*1e-3
dt = 8e-6
srLim = 0.1 # T/pix/s

print(f"srLim = {srLim:.4f} T/pix/s")

lstRho = linspace(0, 0.5, 10000)
lstDk = sqrt(srLim*gamma*lstRho*dt**2)
getDkRho = interp1d(append(lstRho, 0.5), append(lstDk, lstDk[-1]))
rho1 = gamma*(srLim*dt)*dt/2
def getDk(rho, tht):
    if rho == 0:
        dk = rho1
    else:
        dk = getDkRho(rho)
    return clip(dk, rho1, 1/sizImg)
def getDrDt(rho, tht):
    return turbo*0.5/(2*pi)/(sizImg/2)
figure()
for idxIt in range(64):
    trjSpiral, lstRho, lstDk = genSpiral(getDk, getDrDt)
    trjSpiral = trjSpiral.reshape([-1,2])
    print(f"numPt={trjSpiral.shape[0]}")

    lstGrad = tranTraj2Grad_MinSR(trjSpiral, dt)
    lstSR = getSlewRate(lstGrad, dt)

    lstDk_Fix = sqrt(srLim/lstSR)
    lstDk_Fix = 1+(lstDk_Fix-1)*1e-1
    lstDk = clip(lstDk*lstDk_Fix, rho1, 1/sizImg)
    getDkRho = interp1d(concatenate(([0], lstRho[:-1], [0.5])), concatenate(([0], lstDk[1:], [min(gamma*srLim*dt, 1/sizImg)])))

    print(f"idxIt={idxIt}")
    if idxIt%16 == 15:
        subplot(121)
        plot(lstRho, lstSR, marker=".")
        title(f"lstSR_{idxIt}")
        subplot(122)
        plot(lstRho, lstDk, marker=".")
        title(f"lstDk_{idxIt}")
        show(block=False)
        pause(1e-1)

figure()
subplot(121)
plot(lstRho, lstSR, marker=".")
title(f"lstSR_End")
subplot(122)
plot(lstRho, lstDk, marker=".")
plot(lstRho, (1/sizImg)*ones_like(lstDk), marker=".")
title(f"lstDk_End")

figure()
plot(trjSpiral[:,0], trjSpiral[:,1], marker=".")
axis("equal"); title("trjSpiral")

show()