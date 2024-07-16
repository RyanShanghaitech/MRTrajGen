from numpy import *
from matplotlib.pyplot import *
from mrtrajgen import *

sizIm = 128

getD = lambda rho, tht: 1/sizIm # function of sampling interval with respect to rho and theta
getK = lambda rho, tht: (8*rho+1)*0.5/(2*pi)/(sizIm/2) # function of ratio of rho/theta
trjSpiral = genSpiral(getD, getK) # derive trajectory

trjRadial = genRadial(linspace(0, pi, int64(sizIm*pi/8)), linspace(-0.5, 0.5, sizIm)) # derive trajectory

figure()
subplot(121)
scatter(trjSpiral[:,0], trjSpiral[:,1], marker='.')
axis("equal"); title("Spiral")
subplot(122)
scatter(trjRadial[:,0], trjRadial[:,1], marker='.')
axis("equal"); title("Radial")
show()