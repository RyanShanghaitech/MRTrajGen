from numpy import *
from scipy.ndimage import convolve1d

a = array([1, 2, 3, 4, 5], dtype=float64)
a = convolve1d(a, concatenate([zeros([7]), (1/2)**arange(8)]), mode="constant", cval=0)
print(a)