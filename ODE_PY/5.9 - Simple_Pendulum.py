# Runge-Kutta 4th-order 2x2 system

import numpy as np
import matplotlib.pyplot as plt

# our functions f_i in the ODE y' = F(t,y)
def f1(t,x,y):
    return y

def f2(t,x,y):
    return -10*np.sin(x)

# the interval of approximation
a = 0
b = 20

# steps and step size
n = 1000
h = (b-a)/n

# initial data
xx = np.pi/2
yy = 4.4
t = a

# Runge-Kutta 4th-order

WX = [xx]
WY = [yy]
k = 0
K = [k]
T = [t]

for i in range(n):
    K11 = h*f1(t,xx,yy)
    K12 = h*f2(t,xx,yy)

    K21 = h*f1(t+h/2, xx+K11/2, yy+K12/2)
    K22 = h*f2(t+h/2, xx+K11/2, yy+K12/2)

    K31 = h*f1(t+h/2, xx+K21/2, yy+K22/2)
    K32 = h*f2(t+h/2, xx+K21/2, yy+K22/2)
    
    K41 = h*f1(t+h, xx+K31, yy+K32)
    K42 = h*f2(t+h, xx+K31, yy+K32)

    m1 = K11+2*K21+2*K31+K41
    m2 = K12+2*K22+2*K32+K42

    xx = xx + m1/6
    yy = yy + m2/6

    WX.append(xx)
    WY.append(yy)

    t = t + h
    T.append(t)

    k = k + 1
    K.append(k)


#print("\n\t Example 1, p. 331-332:  \n\t 4th-Order Runge-Kutta Method on the IVP")
#print("\n\t p'(t) = -10sin(theta(t))\t\t p(0) = -1")
#print("\n\t theta'(t) = p(t)\t\t theta(0) = 1")
#print("\n\t n\tw_n \t\tx(t_n)\t\tE(t_n)")
#print("\t---------------------------------------------------\n")
#for i in range(len(WX)):
#    print("\t %d"%K[i],"\t%0.7f"%WX[i])

#print("\n\t n\tw_n \t\ty(t_n)\t\tE(t_n)")
#print("\t---------------------------------------------------\n")
#for i in range(len(WY)):
#    print("\t %d"%K[i],"\t%0.7f"%WY[i])


# Plot them

B, BX = plt.subplots()
B1 = BX.plot(WX,WY,"r",linewidth=1)
#BX.legend(("y(t) = -exp(t)/2 + t^2 + 2t + 1\nthe exact solution","PL(t) = the PL Modified Euler\napproximation function"))
BX.set_title("Phase space flow for simple pendulum using\n4th-Order Runge-Kutta")
plt.axvline(color="k",linewidth=1)
plt.axhline(color="k",linewidth=1)
#B3 = BX.plot(WW,W,'o',color="b")
plt.show()
