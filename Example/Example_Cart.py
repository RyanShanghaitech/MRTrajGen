from numpy import *
from matplotlib.pyplot import *
from mrtrjgen import *

nPix = 128

# generate trajectory of radial
arrK = genCart(nPix) # derive trajectory
nPE, nRO, _ = arrK.shape

# show trajectories
figure()
title("Cartesian")
for iPE in range(nPE):
    plot(arrK[iPE,:,0], arrK[iPE,:,1], ".-")
    axis("equal")
    ylim(nPix//2, -nPix//2)
    xlim(-nPix//2, nPix//2)
    draw()
    pause(0.1)

show()