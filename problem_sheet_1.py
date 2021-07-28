"""
Problem Sheet 1

 Author: P Fogg

28-07-2021
"""
import numpy as np
import matplotlib.pyplot as plt



# Input parameters
L = 0.001 #m
d = 5e-6 #m
k = 200 #W/mK
hbar = 1000 #W/m2K

Ta = 293.15 #K
T0 = 353.15 #K
T4 = 343.15 #K

# Calculated parameters
P = np.pi*d
A = np.pi*d**2 /4

print('Number of discretisation points?    ')
n = int(input())
#n=5
dx = L/n

beta = np.sqrt((hbar*P)/(k*A))
sigma = -2 -beta**2 * dx**2

# Analytical soln

def AnalyitcalTheta(x):
    return ( (T4-Ta)*np.sinh(beta*x) + (T0-Ta)*np.sinh(beta*(L-x)) )/(np.sinh(beta*L))
    
# Set up arrays
diagonalL = np.array([0]+[1]*(n-2))
diagonal = np.array([1]+[sigma]*(n-2)+[1])
diagonalU = np.array([1]*(n-2)+[0])


M1 = np.diag(diagonalL, k=1)
M2 = np.diag(diagonal, k=0)
M3 = np.diag(diagonalU, k=-1)

M = M1+M2+M3

# Tutor method
"""
n = 5

M = np.diag(np.array([1]+[sigma]*(n-2)+[1])) +\
    np.diag(np.array([1]*(n-2)+[0]),k=-1) +\
    np.diag(np.array([0]+[1]*(n-2)),k=1);
"""    
# np.lin alg
b = np.array([T0-Ta]+[0]*(n-2)+[T4-Ta])
invM = np.linalg.inv(M)
theta = np.linalg.solve(M,b) # Not working
"""
invM = np.linalg.inv(M)
theta = np.matmul(invM,b)
"""

# Error/grids

x = np.linspace(0,L,n)
thetaAnalytical = AnalyitcalTheta(x)

fig,ax = plt.subplots()

ax.plot(x, thetaAnalytical, '.',
         x, theta, '-.')
ax.grid()
ax.set_title('PBL1: Temperature Distibution in a Hot Wire')
ax.set_ylabel('Temperature [K]')
ax.set_xlabel('Wire length [m]')
ax.legend(['Analytical Solution','Numerical Solution'])
ax.annotate(['dx=',str(L/n),'m'],xy=(0.0002,57.5))
plt.show()
