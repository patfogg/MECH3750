# MECH3750: ICT Week 1
#
# Author: P Fogg 27-07-2021
#
# Methods 
import numpy as np
import matplotlib.pyplot as plt

def f1(x):
    return x*x*x
# create an array for x and f(x)
h = 1
x = np.arange(-5,6,h,dtype=float)
f = f1(x)

# calculate forwards difference.

forwardD = np.zeros_like(f)
forwardD[:-1] = (f[1:]-f[:-1])/h


fig,ax = plt.subplots()

ax.plot(x, f,
        x[:-1], forwardD[:-1])
ax.grid()
plt.show()
