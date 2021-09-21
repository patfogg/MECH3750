"""
MECH3750: Engineering Analysis
Part II: Engineering Applications
PBL07: Heat Conduction in Slender Rod

Forward and backward difference solutions of the diffusion equation as a heat
conduction problem in one dimension. It is assumed that the slender rod is 
insulated along its length and that the temperature does not vary across the
rod's cross-section. The PDE that models this situation is,

dT/dt = D(d^2 T/dx^2)

where D is the thermal diffusivity, D = k/(rho*Cp) [m^2/s], and boundaries
exist at x = 0 and x = L. The initial temperature distribution is defined as,

T(x,0) = 30(cos((n pi x)/L) + 2)
"""

__author__ = 'Christopher Leonardi'
__version__ = '0.0.1'
__date__ = '13/09/2021'

import numpy as np
import matplotlib.pyplot as plt

def rod_temp_1d(L, k, rho, Cp, nx, nt, dt, bcT, bcV):
    """
    Calculation of the transient temperature distribution in a 1D rod
    L [m]       : The length of the steel rod in the x-direction
    k [W/mK]    : The thermal conductivity of the steel rod
    rho [kg/m^3]: The density of the steel rod
    Cp [m]      : The specific heat of the steel rod
    nx [-]      : Number of calculation points along rod, endpoints inclusive
    nt [-]      : The number of simulated timesteps
    dt [s]      : The timestep used in the finite difference scheme
    bcT         : List of LHS and RHS boundary types - D = Dirichlet
                                                       N1 = 1st Order Neumann
                                                       N2 = 2nd Order Neumann
                                                       P = Periodic
    bcV         : List of LHS and RHS boundary value - D = Temperature [K]
                                                       N = Heat flux [W]
    """
    
    # Create the data structures used in the explicit and implicit solutions
    xa = np.linspace(0., L, nx)             # 1D array of points along body
    T_EX = np.array([np.zeros_like(xa),np.zeros_like(xa)])
    T_IM = np.zeros_like(T_EX)
    
    # Calculate parameters used in the explicit and implicit solutions
    D = k/(rho*Cp)                          # Thermal diffusivity of the body
    dx = L/(nx-1)                           # Spatial discretisation
    sigma = (D*dt)/(dx**2)                  # Courant number of the system
    I = np.eye(nx)                          # Identity matrix
    D2M = np.diag(np.full(nx,-2)) + \
        np.diag(np.full(nx-1,1),k=-1) + \
        + np.diag(np.full(nx-1,1),k=1)      # 2nd order differential operator
    
    # Calculate the initial temperature distribution in the body
    T_0 = 30*(np.cos(np.pi*xa/L) + 2)
    F_0 = 'T(x,0)=30(cos((n pi x)/L)+2)'
    
    # EXPLICIT IMPLEMENTATION
    
    # Initialise the array used for the explicit update scheme
    T_EX[0,:] = T_0
    if bcT[0] == 'D':
        T_EX[0,0] = bcV[0]
    if bcT[1] == 'D':
        T_EX[0,-1] = bcV[1]
    
    # Loop over the required number of timesteps and perfrom explicit update
    for i in range(1, nt+1):
        # Calculate indexes for the 'new' and 'old' row data
        new = i%2
        old = np.abs(1-new)
        # Explicit temperature update at each internal node
        T_EX[new,1:-1] = sigma*T_EX[old,:-2] + \
            (1-2*sigma)*T_EX[old,1:-1] + sigma*T_EX[old,2:]
        # Assign the boundary conditions at the LHS
        if bcT[0] == 'D':
            T_EX[new,0] = bcV[0]
        elif bcT[0] == 'N1':
            T_EX[new,0] = T_EX[new,1]
        elif bcT[0] == 'N2':
            T_EX[new,0] = (4 * T_EX[new,1] - T_EX[new,2])/3
        elif bcT[0] == 'P':
            T_EX[new,0] = sigma*T_EX[old,-1] + \
                (1-2*sigma)*T_EX[old,0] + sigma*T_EX[old,1]
        else:
            print('Invalid boundary condition specification on LHS')
        # Assign the boundary conditions at the RHS
        if bcT[1] == 'D':
            T_EX[new,-1] = bcV[1]
        elif bcT[1] == 'N1':
            T_EX[new,-1] = T_EX[new,-2]
        elif bcT[1] == 'N2':
            T_EX[new,-1] = (4 * T_EX[new,-2] - T_EX[new,-3])/3
        elif bcT[1] == 'P':
            T_EX[new,-1] = sigma*T_EX[old,-2] + \
                (1-2*sigma)*T_EX[old,-1] + sigma*T_EX[old,0]
        else:
            print('Invalid boundary condition specification on RHS')

    # IMPLICIT IMPLEMENTATION
    
    # Initialise the array used for the implicit update scheme
    T_IM[0,:] = T_0
    if bcT[0] == 'D':
        T_IM[0,0] = bcV[0]
    if bcT[1] == 'D':
        T_IM[0,-1] = bcV[1]
    
    # Loop over the required number of timesteps and perfrom implicit update
    for i in range(1, nt+1):
        # Calculate indexes for the 'new' and 'old' row data
        new = i%2
        old = np.abs(1-new)
        # Assemble the implicit operator matrix
        H = I - sigma*D2M
        # Assign the boundary conditions at the LHS
        H[0,:] = 0
        if bcT[0] == 'D':
            H[0,0] = 1
            T_IM[old,0] = bcV[0]
        elif bcT[0] == 'N1':
            H[0,0] = -1
            H[0,1] = 1
            T_IM[old,0] = 0
        elif bcT[0] == 'N2':
            H[0,0] = -3
            H[0,1] = 4
            H[0,2] = -1
            T_IM[old,0] = 0
        elif bcT[0] == 'P':
            print('Not yet implemented!')
        else:
            print('Invalid boundary condition specification on LHS')
        # Assign the boundary conditions at the RHS
        H[-1,:] = 0
        if bcT[1] == 'D':
            H[-1,-1] = 1
            T_IM[old,-1] = bcV[1]
        elif bcT[1] == 'N1':
            H[-1,-1] = -1
            H[-1,-2] = 1
            T_IM[old,-1] = 0
        elif bcT[1] == 'N2':
            H[-1,-1] = -3
            H[-1,-2] = 4
            H[-1,-3] = -1
            T_IM[old,-1] = 0
        elif bcT[1] == 'P':
            print('Not yet implemented!')
        else:
            print('Invalid boundary condition specification on RHS')
        # Implicit temperature update at each internal node
        T_IM[new,:] = np.linalg.solve(H,T_IM[old,:])
    
    # Generate graphs of the temperature distribution along the body
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

    return xa, T_0, T_EX, T_IM, sigma, H

