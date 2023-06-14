# Runge-Kutta 4th-order 2x2 system

import numpy as np
import matplotlib.pyplot as plt

# our functions f_i in the ODE y' = F(t,y)
def f1(t,x,y):
    return y

def f2(t,x,y):
    return 3*y-2*x+6*np.exp(-t)

# exact solutions
def g1(t):
    return 2*np.exp(2*t)-np.exp(t)+np.exp(-t)

def g2(t):
    return 4*np.exp(2*t)-np.exp(t)-np.exp(-t)

# the interval of approximation
a = 0
b = 1

# initial data
x = 2
y = 2
t = a

# steps and step size
n = 10
h = (b-a)/n

# Runge-Kutta 4th-order
WX = [x]
WY = [y]
k = 0
K = [k]
X = [x]
Y = [y]
EX = [0]
EY = [0]
M = [0]
T = [t]

for i in range(n):
    K11 = h*f1(t,x,y)
    K12 = h*f2(t,x,y)

    K21 = h*f1(t+h/2, x+K11/2, y+K12/2)
    K22 = h*f2(t+h/2, x+K11/2, y+K12/2)

    K31 = h*f1(t+h/2, x+K21/2, y+K22/2)
    K32 = h*f2(t+h/2, x+K21/2, y+K22/2)
    
    K41 = h*f1(t+h, x+K31, y+K32)
    K42 = h*f2(t+h, x+K31, y+K32)

    m1 = K11+2*K21+2*K31+K41
    m2 = K12+2*K22+2*K32+K42

    x = x + m1/6
    y = y + m2/6

    WX.append(x)
    WY.append(y)

    t = t + h
    T.append(t)

    X.append(g1(t))
    Y.append(g2(t))

    EX.append(abs(x-g1(t)))
    EY.append(abs(y-g2(t)))

    k = k + 1
    K.append(k)

print("\n\t Quiz 8:  \n\t 4th-Order Runge-Kutta Method on the IVP")
print("\n\t y1'(t) = y2(t)\t\t\t y1(0) = 2")
print("\t y2'(t) = 3y2-2y1+6exp(-t)\t y2(0) = 2")
print("\n\t n\tw_n \t\ty1(t_n)\t\tE(t_n)")
print("\t---------------------------------------------------\n")
for i in range(len(WX)):
    print("\t %d"%K[i],"\t%0.7f"%WX[i],"\t%0.7f"%X[i],"\t%0.7f"%EX[i])

print("\n\t n\tw_n \t\ty2(t_n)\t\tE(t_n)")
print("\t---------------------------------------------------\n")
for i in range(len(WY)):
    print("\t %d"%K[i],"\t%0.7f"%WY[i],"\t%0.7f"%Y[i],"\t%0.7f"%EY[i])
