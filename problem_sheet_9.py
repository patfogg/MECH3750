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
    n = 8 

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
        if flag == 'a':
            U_initial = (x*(x-L))/L**2
        elif flag == 'b':
            U_initial = np.where(x<=L/4,x,8/L*x)            
            #while x <= L/4:
            #    U_initital = 8/L*x
            U_initial = np.where((x>L/4) & (x<5*L/8),x,32 - 8/L*x)
            #while L/4 < x < 5*L/4:
            #    U_inititial = 32 - 8/L*x
            U_initial = np.where(x>=5*L/8,x,1/3/L -29/24)
            #while x >= 5*L/8:
            #    U_inititial = 1/3/L -29/24
        elif flag == 'c':
            U_initial = np.zeros_like(x)
        elif flag == 'd':
            U_initial = np.zeros_like(x)
            U_initial[np.where(x<=L/2)] = x[np.where(x<=L/2)]/24
            U_initial[np.where(x>L/2)] = 1-x[np.where(x>L/2)]/24
        return U_initial
    

    def U1func(n):
        return 0
    U0 = np.zeros(n)
    x = np.linspace(0,L,n)
    #print(U0func(x,'a'))
    #print(U0func(x,'b'))
    #print(U0func(x,'c'))
    print(U0func(x,'d'))
    plt.plot(x,U0func(x,'d'))
    plt.show()
    #for i in range(n):
    #    U0[i]=U0func(x[i])
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
    nt = 3
    u = np.dot(1/2*M,U0) + dt*U1
    u_old = U0
    record = np.zeros((nt,n))
    record[0,:] = U0
    
    for i in range(1,nt):
        u_new = M @ u - u_old
        u_old = u
        u = u_new
        record[i,:] = u
        
    fig, ax = plt.subplots()
    for j in range(len(record)):
        ax.plot(x,record[j,:],'--')
        ax.grid()


    plt.show()




    
