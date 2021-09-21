"""
MECH3750 Week 9 ict

Author: P Fogg
"""

import numpy as np
import matplotlib.pyplot as plt

TOL = np.finfo(float).resolution #1/(x+TOL) : x->0

Lx = 0.1 #m
Ly = 0.075 #m
T0 = 323.15 #K
rho = 7.86 #kg/m3
cp = 490 #J/kg/K
k = 54 #W/m/K


#a.
D = k/(rho*cp)
# From Alex 
dx = dy = 0.025
dt = 100

sigma_x = D*dt/dx**2
sigma_y = D*dt/dy**2

nx = int((Lx+TOL) // dx) +1
ny = int((Ly+TOL) // dy) +1

T =(T0) * np.ones(nx*ny)

#BC Config
bot = np.array([110,100,90,80,70]) + 273.15
top = np.array([0,10,20,30,40]) + 273.15
left = np.array([110,65,25,0]) + 273.15
right = np.array([70,65,50,40]) + 273.15

T[0:nx] = bot
T[-nx:] = top
T[0:nx*(ny+1):nx] = left
T[nx-1:nx*ny:nx] = right


imp_matrix = np.zeros((nx*ny,nx*ny),dtype=float)

for i in range(ny):
    for j in range(nx):
        s = i*nx + j
        #skip boundary nodes
        if i==0 or i==ny-1:
            imp_matrix[s,s] = 1
            continue
        if j==0 or j==nx-1:
            imp_matrix[s,s] = 1
            continue
            
        # Finite difference equ
        
        imp_matrix[s,   s]     = (1+2*sigma_x + 2*sigma_y)
        imp_matrix[s,   s-1]   = -sigma_x
        imp_matrix[s,   s+1]   = -sigma_x
        imp_matrix[s,   s-nx]  = -sigma_y
        imp_matrix[s,   s+nx]  = -sigma_y

recording = []

recording.append(T)

for _ in range(dt):
    T_new = np.linalg.solve(imp_matrix, T)
    T = T_new.copy()
    recording.append(T)

plt.imshow(T.reshape(ny, nx) - 273.15, origin='lower')
#plt.show()  


#b.
#change to use sparese matrix
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve
dx = dx/2
dy = dy/2

sigma_x = D*dt/dx**2
sigma_y = D*dt/dy**2

nx = int((Lx+TOL) // dx) +1
ny = int((Ly+TOL) // dy) +1
print(f'The grid is {nx} wide, and {ny} tall')
T = (T0) * np.ones(nx*ny) #redefine T
# Assume BC's are linear between points
bot = np.array([110,105,100,95,90,85,80,75,70]) + 273.15
top = np.array([0,5,10,15,20,25,30,35,40]) + 273.15
left = np.array([110,87.5,65,45,25,12.5,0]) + 273.15
right = np.array([70,65,60,55,50,45,40]) + 273.15

T[0:nx] = bot
T[-nx:] = top
T[0:nx*(ny+1):nx] = left
T[nx-1:nx*ny:nx] = right

imp_matrix = (
    (1+2*sigma_x + 2*sigma_y) * sps.eye(nx*ny, k=0)
    -sigma_x * sps.eye(nx*ny, k=-1)
    -sigma_x * sps.eye(nx*ny, k=+1)
    -sigma_y * sps.eye(nx*ny, k=-nx)
    -sigma_y * sps.eye(nx*ny, k=+nx)
)

recording = []

recording.append(T)

for _ in range(dt):
    T_new = spsolve(imp_matrix,T)
    T = T_new.copy()
    recording.append(T)

plt.imshow(T.reshape(ny, nx) - 273.15, origin='lower')
plt.show()  

