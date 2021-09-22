"""
MECH3750: PBL Week 9

Author: P Fogg
"""

import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    # Constraints
    L = 48/100 #m
    c = np.sqrt(4/10**4) #c**2 = [m2/s2]
    n = 8 
    
    dx = L/n
    dt = 1
    
    sigma = c * dt/dx
    # Boundary Conditions
        # Dirichlet BC -> u{M+1} = u{m}
    def U0func(x):
        L = 0.48
        return (x*(x-L))/L**2
    def U1func(n):
        return 0
    U0 = []
    x = np.linspace(0,L,n)
    for i in range(n):
        U0.append(U0func(x[i]))
    U1 = np.zeros(n)
    # Build Matrix
    M = np.diag(np.full(n, 2-2*sigma**2)) +\
        np.diag(np.full(n-1,sigma**2),k=1) +\
        np.diag(np.full(n-1,sigma**2),k=-1)
    M[0,:] = 0; M[0,0] = 1; M[-1,:]=0; M[-1,-1] = 1
        
    # from lect. notes
        # d2y/dt2 = c**2 d2y/dx2
        # gen. CD form for first time-step
        # u{1}=1/2 * M * U{0} + dt * U{1}
        # gen . CD form for subsequent steps
        # u{m+1} = M * u{m} - u{m-1}
    nt = 48
    u = np.dot(1/2*M,U0) + dt*U1
    for i in range(1,nt+1):
        u_new = M @ u - u_old
        u = u_new
        u_old = u
