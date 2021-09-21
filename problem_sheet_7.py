"""
PBL Sheet 7

Author: P Fogg
"""
import numpy as np
import matplotlib.pyplot as plt

def Q1():
    "Explicit, Forward Difference transient heat solver for slender rod"
    
    T0 = 90 + 273.15 #K
    TL = 30 + 273.15 #K

    def T0(x):
        return 30 * (np.cos((np.pi*x)/1)+2)

    n = 5  # Number of nodes
    dt = 1000 #s
    def explicit_FD(dx,dt,nt):
        "dx = grid spaceing [m]"
        "dt = time-step [s]"
        "nt = number of time steps"
        # Constraints
        L = 1.0 #m
        rho = 7860 #kg/m3
        cp = 490 #J/kg/K
        k = 54 #W/m/K
        n = int(L/dx)
        D = k/(rho*cp)
    
        sigma =D * dt/dx**2

        x = np.linspace(0,L,n)
    
        # Build Temperature Matrix
        T = np.ones(5) * 273.15
        # Configure Boundary Conditions
        for i in range(n):
            T[i] = T0(float(x[i]))
        print('Initial Temperature array')
        Tint = T
        print(T)
        # Explicit Matrix
        M = np.diag(np.full(n,1-2*sigma)) + \
            np.diag(np.full(n-1,sigma),k=-1) + \
            + np.diag(np.full(n-1,sigma),k=1)
        M[0,:] = 0; M[0,0] = 1; M[-1,:] = 0; M[-1,-1] = 1

        recording = []
        # Loop over time
        for _ in range(nt+1):
            Tnew = M @ T
            T = Tnew
            recording.append(T)
        print(f'Temperature array after {nt} time steps')
        print(f'  with a time step of {dt}')
        print(T)
        return T, Tint, x
    print('a)+b)')
    explicit_FD(0.2,1000,100)
    print('c)')
    explicit_FD(0.2,2000,100)
    explicit_FD(0.2,3000,100)

    T,Tint,x = explicit_FD(0.2,1000,100)
    
    F_0 = 'T(x,0)=30(cos((n pi x)/L)+2)'
    fig,ax = plt.subplots()
    fig.suptitle('Heat Conduction in a Steel Rod: '+F_0,fontsize=8)
    ax.plot(x,Tint,'b-',x,T,'g-',linewidth=0.5)
    ax.legend(['T0','EX'],fontsize=6)
    ax.set_xlabel('X-Coordinate (m)',fontsize=6)
    ax.set_ylabel('Temperature (K)',fontsize=6)
    ax.label_outer()
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    ax.grid()
    plt.show()
    
#-------------------------#
def Q2():
    "Implicit, Backwards Difference transient heat solver for slender rod"
        
    T0 = 90 + 273.15 #K
    TL = 30 + 273.15 #K

    def T0(x):
        return 30 * (np.cos((np.pi*x)/1)+2)

    n = 5  # Number of nodes
    dt = 1000 #s
    def implicit_FD(dx,dt,nt):
        "dx = grid spaceing [m]"
        "dt = time-step [s]"
        "nt = number of time steps"
        # Constraints
        L = 1.0 #m
        rho = 7860 #kg/m3
        cp = 490 #J/kg/K
        k = 54 #W/m/K
        n = int(L/dx)
        D = k/(rho*cp)
    
        sigma =D * dt/dx**2

        x = np.linspace(0,L,n)
    
        # Build Temperature Matrix
        T = np.ones(5) * 273.15
        # Configure Boundary Conditions
        for i in range(n):
            T[i] = T0(float(x[i]))
        print('Initial Temperature array')
        Tint = T
        print(T)
        # Implicit Matrix
        M = np.diag(np.full(n,1+2*sigma)) + \
            np.diag(np.full(n-1,-sigma),k=-1) + \
            + np.diag(np.full(n-1,-sigma),k=1)
        M[0,:] = 0; M[0,0] = 1; M[-1,:] = 0; M[-1,-1] = 1

        recording = []
        # Loop over time
        for _ in range(nt+1):
            Tnew = np.linalg.solve(M, T)
            T = Tnew
            recording.append(T)
        print(f'Temperature array after {nt} time steps')
        print(f'  with a time step of {dt}')
        print(T)
        return T, Tint, x
    
    T,Tint,x = implicit_FD(0.2,1000,100)
    F_0 = 'T(x,0)=30(cos((n pi x)/L)+2)'
    fig,ax = plt.subplots()
    fig.suptitle('Heat Conduction in a Steel Rod: '+F_0,fontsize=8)
    ax.plot(x,Tint,'b-',x,T,'g-',linewidth=0.5)
    ax.legend(['T0','IM'],fontsize=6)
    ax.set_xlabel('X-Coordinate (m)',fontsize=6)
    ax.set_ylabel('Temperature (K)',fontsize=6)
    ax.label_outer()
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    ax.grid()
    plt.show()
    

"""
    fig,axs = plt.subplots(2,1)
    fig.suptitle('Heat Conduction in a Steel Rod: '+F_0,fontsize=8)
    axs[1].plot(xa,T_0,'b-',xa,T_EX[new,:],'g-',linewidth=0.5)
    axs[1].legend(['T0','EX'],fontsize=6)
    axs[0].plot(xa,T_0,'b-',xa,T_IM[new,:],'y-',linewidth=0.5)
    axs[0].legend(['T0','IM'],fontsize=6)
    for ax in axs.flat:
        ax.set_xlabel('X-Coordinate (m)',fontsize=6)
        ax.set_ylabel('Temperature (K)',fontsize=6)
        ax.label_outer()
        ax.tick_params(axis='x', labelsize=6)
        ax.tick_params(axis='y', labelsize=6)
        ax.grid()
    fig.savefig('pbl07RodTemp1D.jpg',dpi=300)
    plt.show()
"""
if __name__ == '__main__':
    questions = (Q1,Q2)
    for question in questions:
        input(f'Press `Enter` to run {question.__name__} ')

        plt.close('all')        # <- Close all existing figures
        question()
        plt.show(block=False)   # <- Allow code execution to continue

    input('Press `Enter` to quit the program.')
