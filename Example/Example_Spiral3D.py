from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import *

import mrtrjgen

gamma = 42.5756e6 # UIH
fov = 0.25
numPix = 128
uTht = 32
uPhi = 32
dt = 10e-6

for idxTht in range(uTht):
    for idxPhi in range(uPhi):
        arrKxyz, arrGxyz = mrtrjgen.genSpiral3DTypeB(100*gamma*fov/numPix, numPix, 2*pi*idxTht/uTht, 2*pi*idxPhi/uPhi, uTht, uPhi, dt, 0.5)
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