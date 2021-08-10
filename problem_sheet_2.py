"""
MECH3750: Problem Sheet 2

Author: P Fogg 04-08-2021
"""
import numpy as np

from array import *
"""
Q1) Let f1 = x1**2 + cos(x2) - 5 = 0 ; f2 = sinh(x1) + sin(x2) - 4 = 0


"""
#f = [x**2 + cos(y) - 5, sinh(x) + sin(y) - 4]#function vector
xvec = np.array([1, 1]) # initial guess
def fvec(xvec):
        x, y = xvec # Split into components
        vector = np.zeros((len(xvec),))
        vector[0] = x**2 + np.cos(y) - 5 # Hard coding given functions, f1
        vector[1] = np.sinh(x) + np.sin(y) - 4 # f2
        return vector

def jacobian(xvec):
        x, y = xvec # Split into components
        matrix = np.zeros((len(xvec), len(xvec)))
        #for i in range(len(xvec)):
        #        for j in range(len(xvec)):
        #                matrix[i,j] = df_i/d_xj
        matrix[0,0] = 2*x
        matrix[0,1] = -np.sin(y)
        matrix[1,0] = np.cosh(x) 
        matrix[1,1] = np.cos(x)
        return matrix

def NumJacobian(xvec, h):
        x, y = xvec
        matrix = np.zeros((len(xvec), len(xvec)))
        #for i in range(len(xvec)):
        #        matrix[0,:] = (fvec([x+h,y])-fvec(xvec))/h
        #        matrix[1,:] = (fvec([x,y+h])-fvec(xvec))/h
        matrix[0,:] = (fvec([x+h,y])-fvec(xvec))/h
        matrix[1,:] = (fvec([x,y+h])-fvec(xvec))/h
        return np.transpose(matrix)

h = 0.02
J=NumJacobian(xvec,h)


def NewtonSolver(xvec, eps, h):
        "Need to add epsilon"
        xold = xvec
        xnew = xvec - np.dot(np.linalg.inv(NumJacobian(xvec, h)),fvec(xvec))
        while np.abs(1-np.amax(xnew)/np.amax(xold)) > 0.001:
            xold = xnew
            Jacobian = NumJacobian(xold, h)
            xnew = xold - np.dot(np.linalg.inv(Jacobian),fvec(xold))
        return xnew
epsilon = 0.01
h = 0.02
print(NewtonSolver(xvec, epsilon, h))


print("Using Analytical Jacobian and initial guess of [1,1], X(1)= ",
      np.around(NewtonSolver(xvec, jacobian(xvec),0.01),4))
