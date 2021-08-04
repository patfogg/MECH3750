"""
MECH3750: Week 2 ICT Problem

Author: P Fogg

Date: 03-08-2021

"""
import numpy as np
import matplotlib.pyplot as plt
# 1. Function to calculate n!

def nFactorial(x):
	"Function finds the factorial of x"
	n = 1.0
	for i in np.linspace(1,x,x):
		n = n * i
		#print(n)
	return n
' Tutor uses range() for both'
' above: range(1,n+1) below: range(n)'
# 2. Function for Taylor Series of exp(x)
'not working out quite right, huge error'
def exptaylor(n,x):
        "Evaluates the Taylor Series of exp(x) from x=0 to O(n)"
        e = 0.0
        for i in np.linspace(0,n,n+1):
                e += + x**i * 1/nFactorial(int(i))
        return e
def expTaylor(n,x):
        "Same as exptaylor but remebers previously calculated values"
        term = 1
        output = term
        for i in range(1,n):
                term *= x / i
                output += term
        return output

# 3. plot of exp(x)
x = np.linspace(-2,2,10)

exp=np.exp(x)
exp[:-1] = np.exp(x[:-1])

fig,ax=plt.subplots(2,2)

ax[0,0].plot(x[:-1],exp[:-1],'--')
ax[0,0].grid()
ax[0,0].set_title('3. Plot of exp(x) on [-2,1]')
ax[0,0].legend(['exp(x)'])


# 4. plot of exp(x) and exptaylor(n,x)
n = 3
expTay = expTaylor(n,x)

ax[0,1].plot(x,exp,'--',
             x,expTay,'.-')
ax[0,1].grid()
ax[0,1].set_title(('4. exp(x) & exptaylor(n,x) on [-2,2] for n=', n))
ax[0,1].legend(['exp(x)','exptaylor(x)'])

# 5. Plot x2 + y2 with 9x9 grid points
""" Tutor aided, 3rd plotting"""
x = np.linspace(-1,1,9)
y = np.linspace(-1,1,9)

X, Y = np.meshgrid(x,y)
Z = X**2 + Y**2


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_surface(X,Y,Z)

# 6. Taylor Series of ln(1+x) 
"d/dx ln(g(x)) = g'(x)/g(x)"
def lnTaylor(n,x):
        term = -1.0
        output = np.log(2.)
        for i in range(1,n):
                term *= (1-x)/(2*i)
                output += term
        return output

data1 = []
data2 = []
data3 = []
data4 = []
data5 = []
data6 = []
data7 = []

for i in range(6):
        h = 1/10**(i+1)
        data1.append(h)
        ln = np.log(1+1+h)
        data2.append(round(ln,4))
        lnTay = lnTaylor(3,(1+h))
        data3.append(round(lnTay,4))
        data4.append( "{:.3e}".format(ln-lnTay))
        data5.append( "{:.3e}".format(ln-lnTay/h**2))
        data6.append( "{:.3e}".format(ln-lnTay/h**3))
        data7.append( "{:.3e}".format(ln-lnTay/h**4))
data = [data1, data2, data3, data4, data5, data6, data7]
data = np.transpose(data)


        
# 6b. Table

headings = ["h","ln(1+x)","Taylor","Error","Error/h2","Error/h3","Error/h4"]

format_row = "{:>12}" * (len(headings)+1)
zero = np.array([""]*6)
print(format_row.format("",*headings))
for head, row in zip(zero, data):
        print(format_row.format(head, *row))

print("Error/h2 - only magnitude changes with remaining two columns")

# 7. Function that takes vector input xi and produces
#    a vector output fi of function values

xvec = np.array([1, 1]) # initial guess
def fvec(xvec):
        x, y = xvec # Split into components
        vector = np.zeros((len(xvec),))
        vector[0] = x**4 + y**4 -1 # Hard coding given functions, f1
        vector[1] = x**2 - y**2 +1 # f2
        return vector

print(fvec(xvec))

def jacobian(xvec):
        x, y = xvec # Split into components
        matrix = np.zeros((len(xvec), len(xvec)))
        #for i in range(len(xvec)):
        #        for j in range(len(xvec)):
        #                matrix[i,j] = df_i/d_xj
        matrix[0,0] = 4*x**3
        matrix[0,1] = 4*y**3
        matrix[1,0] = 2*x 
        matrix[1,1] = -2*x
        return matrix

print("7. Analytic Jacobian from derivatives")
print(jacobian(xvec))
#print(np.linalg.inv(jacobian(xvec)))
#print(np.linalg.inv(jacobian(xvec))*fvec(xvec))

# 8. Numerical Jacobian solver

def NumJacobian(xvec, h):
        x, y = xvec
        matrix = np.zeros((len(xvec), len(xvec)))
        #for i in range(len(xvec)):
        #        matrix[0,:] = (fvec([x+h,y])-fvec(xvec))/h
        #        matrix[1,:] = (fvec([x,y+h])-fvec(xvec))/h
        matrix[0,:] = (fvec([x+h,y])-fvec(xvec))/h
        matrix[1,:] = (fvec([x,y+h])-fvec(xvec))/h
        return np.transpose(matrix)


f = []#function vector
h = 0.02
J=NumJacobian(xvec,h)
print("8. Numerical approximation of Jacobian using central differences of partials"),print(J)
# 9. Newton solver

def NewtonSolver(xvec, Jacobian, epsilon):
        "Need to add epsilon"
        return xvec - np.dot(np.linalg.inv(Jacobian),fvec(xvec))


print("Analytical Jacobian ",NewtonSolver(xvec, jacobian(xvec),0.01))

print("Numerical Jacobian ",NewtonSolver(xvec, NumJacobian(xvec, h),0.01))

#plt.show()
