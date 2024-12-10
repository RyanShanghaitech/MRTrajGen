from numpy import *
from matplotlib.pyplot import *
import mrtrjgen

dtGrad = 10e-6
dtADC = 2.5e-6
sr = 100
gamma = 42.58e6
fov = 0.5
nPix = 256

_, arrG = mrtrjgen.genSpiral2D(0.5/(2*pi)/8, 1e-3, 0, 0.5, sr*gamma*fov/nPix, dtGrad, 1e2)
nPt, _ = arrG.shape
arrG_Delay = mrtrjgen.delayGrad(arrG, 1/4)

arrK_Delay, _ = mrtrjgen.intpTraj(arrG_Delay, dtGrad, dtADC)
arrK, _ = mrtrjgen.intpTraj(arrG, dtGrad, dtADC)

figure()
plot(arrK[:,0], arrK[:,1], ".-")
plot(arrK_Delay[:,0], arrK_Delay[:,1], ".-")
axis("equal")

figure()
subplot(211)
plot(arrG[:,0], ".-")
plot(arrG_Delay[:,0], ".-")
subplot(212)
plot(arrG[:,1], ".-")
plot(arrG_Delay[:,1], ".-")

show()