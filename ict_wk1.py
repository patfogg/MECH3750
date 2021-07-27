# MECH3750: ICT Week 1
#
# Author: P Fogg 27-07-2021
#
# Calculating the first Derivative using a couple methods
import numpy as np
import matplotlib.pyplot as plt


# Functions to find first div's of:

def f1(x):
    return x*x*x
def f2(x):
    return 3*x**2-2*x
def f3(x):
    "Input in degrees"
    return np.sin(x*np.pi/180)

# create an array for x and f(x)
h = 1
x = np.arange(-5,6,h,dtype=float)
print("Function of interest? (a,b or c)")
f = input()
if f == 'a':
    f = f1(x)
elif f == 'b':
    f = f2(x)
elif f == 'c':
    f = f3(x)

# calculate forwards difference.

forwardD = np.zeros_like(f)
forwardD[:-1] = (f[1:]-f[:-1])/h

# bit lost from here down, seems like there all the same really.
'--------------------------------'
backwardD = np.zeros_like(f)
backwardD[1:] = (f[:-1]-f[1:])/h

centralD = np.zeros_like(f)
centralD = (f[1:]-f[:-1])/(2*h)
'---------------------------------'
fig,ax = plt.subplots()

ax.plot(x, f,
        x[:-1], forwardD[:-1],
        x[1:], backwardD[1:],
        x, centralD)
ax.legend('f(x)','forwardD','backwardD','centralD')
ax.grid()
#plt.show()
