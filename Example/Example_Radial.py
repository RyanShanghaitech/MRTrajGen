from numpy import *
from matplotlib.pyplot import *
from mrtrjgen import *

nPix = 128

# generate trajectory of radial
arrK = genRadial(linspace(0, 2*pi, int64(nPix*pi/1)), linspace(0, 0.5, nPix//2)) # derive trajectory

# show trajectories
figure()
for idxTraj in range(arrK.shape[0]):
    plot(arrK[idxTraj,:,0], arrK[idxTraj,:,1], marker='.')
axis("equal"); title("Radial")

show()