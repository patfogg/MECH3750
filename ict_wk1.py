"""
MECH3750: ICT Week 1

Author: P Fogg 27-07-2021

Calculating the first Derivative using a couple methods
"""
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
# Known divs

def f1prime(x):
    return 3*x**2
def f2prime(x):
    return 6*x-2
def f3prime(x):
    return np.cos(x*np.pi/180)

# create an array for x and f(x)
h = 1
x = np.arange(-5,6,h,dtype=float)
print("Function of interest? (a,b or c)")
f = input()
if f == 'a':
    f = f1(x)
    fd = f1prime(x)
elif f == 'b':
    f = f2(x)
    fd = f2prime(x)
elif f == 'c':
    f = f3(x)
    fd = f3prime(x)

# calculate forwards difference.

forwardD = np.zeros_like(f)
forwardD[:-1] = (f[1:]-f[:-1])/h

# bit lost from here down, seems like there all the same really. 28/7 - This is true, it's the index of the LHS that is of importance
'--------------------------------'
backwardD = np.zeros_like(f)
backwardD[1:] = (f[1:]-f[:-1])/h

centralD = np.zeros_like(f)
centralD[1:-1] = (f[2:]-f[:-2])/(2*h)
'---------------------------------'

# Calculate the error



errors = np.zeros((3,x.size))
errors[0,:-1] = fd[:-1]-forwardD[:-1] 
errors[1,1:] = fd[1:]-backwardD[1:]
errors[2,1:-1] = fd[1:-1]-centralD[1:-1]
fig,ax = plt.subplots()

ax.plot(x, f,
        x[:-1], forwardD[:-1],'-.',
        x[1:], backwardD[1:],'.',
        x[1:-1], centralD[1:-1], '--')

ax.legend(['f(x)','forwardD','backwardD','centralD'])
ax.grid()

fig2,ax2 = plt.subplots()

ax2.plot(x[:-1], errors[0,:-1],'-.',
        x[1:], errors[1,1:],'.',
        x[1:-1], errors[2,1:-1], '--')
ax2.legend(['forwardD','backwardD','centralD'])
ax2.grid()
plt.show()
