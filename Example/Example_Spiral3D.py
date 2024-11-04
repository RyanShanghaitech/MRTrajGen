from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import *
import os
import mrtrjgen

gamma = 42.5756e6 # UIH
pi = arccos(-1)
fov = 0.25
numPix = 128

for idxTht in range(32):
    for idxPhi in range(32):
        arrKxyz, arrGxyz = mrtrjgen.genSpiral3DTypeA(100*gamma*fov/numPix, numPix, 2*pi*idxTht/32, 2*pi*idxPhi/32, 32, 32, 2.5e-6, 0.5)
        save(f"../Resource/arrKxyz_tht{idxTht:02d}_phi{idxPhi:02d}.npy", arrKxyz)

figure()
subplot(111, projection="3d")
plot(arrKxyz[:,0], arrKxyz[:,1], arrKxyz[:,2], ".-")

figure()
subplot(311)
plot(arrGxyz[:,2], ".-")
subplot(312)
plot(arrGxyz[:,1], ".-")
subplot(313)
plot(arrGxyz[:,0], ".-")

show()