"""
MECH3750: ICT wk3


Author: P Fogg 10-09-2021

Least squares: Methods
"""

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    x,f = np.loadtxt('data_points.txt', unpack=True)

    fig, ax = plt.subplots(2,1)

    n = len(x) # data points

    p0 = x**0
    p1 = x**1
    p2 = x**2
    p = [p0,p1,p2]
    P = np.zeros([3,3])

    for i in range(3):
        for j in range(3):
            P[i,j] = p[i] @ p[j]

    y = np.zeros(3)
    for i in range(3):
        y[i] = p[i] @ f

    coeffs = np.linalg.solve(P,y)
    """
    Alex method
    P = np.array([
        x**0,x**1,x**2]).T

    #(P.T @ P )@ coeffs = P.T @ f

    coeffs = np.linalg.solve(
        (P.T @ P), (P.T @ f))
    """

    domain = np.linspace(0,3)
    domain_basis = np.array([
        domain**0,
        domain**1,
        domain**2
        ])

    model = coeffs @ domain_basis

    ax[0].plot(x,f,'o',
            domain, model, '-.')
    ax[0].grid()
    ax[0].set_title('1. Quad fit Least Square of data_points.txt')

    x = np.linspace(-np.pi,np.pi)
    f = x*(np.pi-x)

    P = np.array([
        np.sin(x),np.sin(2*x),np.sin(3*x)]).T
    coeffs = np.array([8/np.pi, 0, 8/(27*np.pi)])
    y = np.zeros(len(x))
    #y = np.linalg.solve(P,coeffs)
    for i in range(len(x)):
        y[i] = P[i,0]*coeffs[0] + P[i,1]*coeffs[1] + P[i,2]*coeffs[2]

    ax[1].plot(x,f,'-',
            x, y, '-.')
    ax[1].grid()
    ax[1].set_title('2. Least Square Approm. of f(x)=x(pi-x)')

    









    plt.show()
