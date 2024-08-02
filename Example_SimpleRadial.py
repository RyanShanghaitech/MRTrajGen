from numpy import *
from matplotlib.pyplot import *
from mrtrajgen import *

sizIm = 16

# generate trajectory of radial
trjRadial = genRadial(linspace(0, 2*pi, int64(sizIm*pi/1)), linspace(0, 0.5, sizIm//2)) # derive trajectory

# show trajectories
figure()
for idxTraj in range(trjRadial.shape[0]):
    plot(trjRadial[idxTraj,:,0], trjRadial[idxTraj,:,1], marker='.')
axis("equal"); title("Radial")

show()