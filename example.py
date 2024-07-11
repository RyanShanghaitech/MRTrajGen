from numpy import *
from matplotlib.pyplot import *
from packTrajGen import modTrajGen

sizIm = 128

funD = lambda rho, tht: 1/sizIm # function of sampling interval with respect to rho and theta
funK = lambda rho, tht: (8*rho+1)*0.5/(2*pi)/(sizIm/2) # function of ratio of rho/theta
traj = modTrajGen.funGenSpiral(funD, funK) # derive trajectory

figure()
plot(traj[:,0], traj[:,1], marker='.')
axis("equal")
show()