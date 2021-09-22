"""
MECH3750: PBL Week 9

Author: P Fogg
"""

import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    # Constraints
    L = 48 #cm
    c = np.sqrt(4) #c**2 = [cm2/s2]
    n = 50

    dx = L/n
    
    def Dt(CFL, dx, c):
        dt = CFL*dx/c
        return dt
    sigma = 0.9
    dt = Dt(sigma, dx, c)
    print(f'Time step = {dt}')
    # Boundary Conditions
        # Dirichlet BC -> u{M+1} = u{m}
    def U0func(x,flag):
        U_initial = np.zeros_like(x)
        if flag == 'a':
            U_initial = (x*(x-L))/L**2
        elif flag == 'b':
            U_initial[np.where(x<=L/4)] = x[np.where(x<=L/4)]/6            
            #while x <= L/4:
            #    U_initital = 8/L*x
            U_initial[np.where((x>L/4) & (x<5*L/8))] = 4 - x[np.where((x>L/4) & (x<5*L/8))]/6
            #while L/4 < x < 5*L/4:
            #    U_inititial = 32 - 8/L*x
            U_initial[np.where(x>=5*L/8)] = x[np.where(x>=5*L/8)]/18 - 8/3
            #while x >= 5*L/8:
            #    U_inititial = 1/3/L -29/24
        elif flag == 'c':
            U_initial = np.zeros_like(x)
        elif flag == 'd':
            U_initial[np.where(x<=L/2)] = x[np.where(x<=L/2)]/24
            U_initial[np.where(x>L/2)] = 1-x[np.where(x>L/2)]/48
        return U_initial
    

    def U1func(x,flag):
        U1_initial = np.zeros_like(x)
        if flag == 'a':
            U1_inital = U1_initial
        elif flag =='b':
            U1_initial = U1_initial
        elif flag == 'c':
            U1_initial = (x*(x-L))/L**2
        elif flag == 'd':
            U1_initial[np.where(x<=L/2)] = -x[np.where(x<=L/2)]/24
            U1_initial[np.where(x>L/2)] = -1*(1-x[np.where(x>L/2)]/48)
            
        return U1_initial
    
    x = np.linspace(0,L,n)
 
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
    nt = 4
    def stringSolver(x,flag):
        u = np.dot(1/2*M,U0func(x,flag)) + dt*U1func(x,flag)
        u_old = U0func(x,flag)
        record = np.zeros((nt,n))
        record[0,:] = U0func(x,flag)
        
        for i in range(1,nt):
            u_new = M @ u - u_old
            u_old = u
            u = u_new
            record[i,:] = u
        return record 


    recordA = stringSolver(x,'a')
    recordB = stringSolver(x,'b')
    recordC = stringSolver(x,'c')
    recordD = stringSolver(x,'d')
    
    fig, ax = plt.subplots(2,2)
    for j in range(nt):
        ax[0,0].plot(x,recordA[j,:],'--',
                     x,recordA[-1,:],'tab:blue')
        ax[0,0].grid()
        ax[0,1].plot(x,recordB[j,:],'--',
                     x,recordB[-1,:])
        ax[0,1].grid()
        ax[1,0].plot(x,recordC[j,:],'--',
                     x,recordC[-1,:])
        ax[1,0].grid()
        ax[1,1].plot(x,recordD[j,:],'--',
                     x,recordD[-1,:])
        ax[1,1].grid()
        
    plt.show()




    
