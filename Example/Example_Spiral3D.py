from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import *
import os
import mrtrjgen

gamma = 42.58e6
fov = 0.25
numPix = 128

for idxTht in range(32):
    for idxPhi in range(32):
        arrKxyz = mrtrjgen.genSpiral3DTypeA(100*gamma*fov/numPix, numPix, 2*pi*idxTht/32, 2*pi*idxPhi/32, 32, 32, 10e-6, 0.5+2/128)
        save(f"./Resource/arrKxyz_tht{idxTht:02d}_phi{idxPhi:02d}.npy", arrKxyz)

figure()
subplot(111,projection="3d")
plot(arrKxyz[:,0], arrKxyz[:,1], arrKxyz[:,2], ".-")

show()