"""
MECH3750: PBL Week 8

Author: P Fogg
"""

import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':
    def Q1():
        T0 = 25 + 273.15 #K
        T1 = 500 + 273.15 #K
        # L, w >> t

        "Want a function to take, dx, dt,nt"
        def Crank_Nicolson(dx,dt,nt):
            # Constrtaints
            T0 = 25 + 273.15 #K
            T1 = 500 + 273.15 #K
        
            n = 10
            D = 1.27e-4 #m2/s
            # sigma value
            sigma = D * dt/dx**2
            print(f'Sigma = {sigma}')
            # sigma check <= 0.8
            if sigma > 0.8:
                print('Sigma too large.')
            
            # Build t matrix
            T = np.ones(n) * 273.15
            # BC's
            T[0] = T0 ; T[-1] = T1
            # implicit M
            IM = np.diag(np.full(n,1+2*sigma)) + \
                np.diag(np.full(n-1,-sigma),k=-1) + \
                + np.diag(np.full(n-1,-sigma),k=1)
            IM[0,:] = 0; IM[0,0] = 1; IM[-1,:] = 0; IM[-1,-1] = 1
            # explicit M
            EX = np.diag(np.full(n,1-2*sigma)) + \
                np.diag(np.full(n-1,sigma),k=-1) + \
                + np.diag(np.full(n-1,sigma),k=1)
            EX[0,:] = 0; EX[0,0] = 1; EX[-1,:] = 0; EX[-1,-1] = 1
            # Crank-Nicolson
                # u{m+1} = inv.M_implicit dot (M_explicit dot u{m})
            for _ in range(nt+1):
                Tnew = np.dot(np.linalg.inv(IM),np.dot(EX,T))
                T = Tnew
            return T

        print(Crank_Nicolson(0.001,0.001,50))
        #
    print(Q1())
