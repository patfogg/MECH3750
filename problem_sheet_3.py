"""
MECH3750: PBL Week 7

Author: P Fogg

Heat transfer in slender rods
"""

import numpy as np
import matplotlib.pyplot as plt
# find intial values
n = 100
L = 1
def T0(x, n):
    L = 1
    dx = L/(n-1)
    return 30*(np.cos((x*np.pi*dx)/L)+2)

Tint = []
for i in range(100):
    Tint.append(T0(i,100))

# Function for general case

def T(dt,n,nt):
    L = 1
    rho = 7860
    Cp = 490
    k = 54
    dx = L/(n-1)

    
    # Big solution matrix
    sigma = dt/dx**2 * k/(rho*Cp)
    


    M1 = np.diag([sigma]*(n-1), k=1)
    M2 = np.diag([1-2*sigma]*n, k=0)
    M3 = np.diag([sigma]*(n-1), k=-1)

    M = M1+M2+M3

    M[0,:] = 0
    M[0,0]=1 
    M[-1,:] = 0
    M[-1,-1] = 1
  
    Told = Tint
    
    Tnew = M @ Told
    x = np.linspace(0,L,n)

    for j in range(nt):
        T = M @ Told
        Told = T

    
    return x, T, sigma

print('dt = 1000')
#print(T(1000,100))
#print('dt = 2000')
#print(T(2000,100))
#print('dt = 3000')
#print(T(3000,100))
#print('dt = 10000')
#print(T(10000,100))

sol = T(100,100,10)
sol1 = T(1,100,10)
sol2 = T(0.1,100,100)



fig, ax = plt.subplots(2,2)

ax[0,0].plot(sol[0],sol[1])
ax[0,0].grid()
ax[0,0].set_title('dt = 100, sigma = '+ str(round(sol[2],2)))
ax[0,0].set_ylabel('Temperature [K]')
ax[0,0].set_xlabel('Length [m]')

ax[0,1].plot(sol1[0],sol1[1])
ax[0,1].grid()
ax[0,1].set_title('dt = 1, sigma = '+ str(round(sol1[2],2)))
ax[0,1].set_ylabel('Temperature [K]')
ax[0,1].set_xlabel('Length [m]')

ax[1,1].plot(sol2[0],sol2[1])
ax[1,1].grid()
ax[1,1].set_title('dt = 0.1, sigma = '+ str(round(sol2[2],2)))
ax[1,1].set_ylabel('Temperature [K]')
ax[1,1].set_xlabel('Length [m]')

plt.show()


