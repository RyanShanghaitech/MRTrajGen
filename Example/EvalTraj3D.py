from numpy import *
from matplotlib.pyplot import *

numPix = 128

arrCart = zeros([numPix,numPix,numPix], int64)

numTht = 32; numPhi = 32
for idxTht in range(numTht):
    print(f"idxTht={idxTht}")
    for idxPhi in range(numPhi):
        arrKxyz = load(f"./Resource/arrKxyz_tht{idxTht:02d}_phi{idxPhi:02d}.npy")
        for k in arrKxyz:
            k = int64((k+0.5)*numPix + 0.5).clip(0, numPix-1)
            arrCart[k[2],k[1],k[0]] = 1

figure()
ax = imshow(arrCart[0,:,:], vmin=0, vmax=1)
for z in range(numPix):
    print(f"z={z}")
    ax.set_data(arrCart[z,:,:])
    draw()
    pause(1/10)