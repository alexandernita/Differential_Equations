# Euler's Method to approximate x' = x + t on [0,2] with 10 steps

import numpy as np
import matplotlib.pyplot as plt

# our function f in the ODE x' = f(t,x)
def f(t,x):
    return t + x

# exact solution
def g(t):
    return 2*np.exp(t)-t-1

# the interval of approximation
a = 0
b = 2

# initial data
x = 1
t = a

# steps and step size
n = 30
h = (b-a)/n

# max error function
def ME(t):
    return h*np.exp(2)*(np.exp(t)-1)

# Euler's Method
W = [x]
k = 0
K = [k]
X = [x]
E = [0]
M = [0]
T = [0]

for i in range(n):
    m = f(t,x)
    t = t + h
    x = x + m*h
    k = k + 1
    T.append(t)
    W.append(x)
    K.append(k)
    X.append(g(t))
    E.append(abs(x-g(t)))
    M.append(ME(t))

print("\n\t n\tw_n \t\tx(t_n)\t\terror\t\terror bound")
print("\t--------------------------------------------------------------------\n")
for i in range(len(W)):
    print("\t %d"%K[i],"\t%0.6f"%W[i],"\t%0.6f"%X[i],"\t%0.6f"%E[i],"\t%0.6f"%M[i])

# create the PL interpolation function
def PL(s):
    for i in range(n):
        aa = T[i]
        bb = T[i+1]
        if (aa <= s) and (s < bb):
            return (-1/h)*W[i]*(s-bb)+W[i+1]*(s-aa)/h
        elif s==T[n]:
            return W[n]
        
# Plot x(t) against PL(t)

# Create the list Z of y-values of PL(t)
t = np.linspace(0,2)
Z = []
for i in range(len(t)):
   Z.append(PL(t[i]))

# Plot them   
WW = np.linspace(0,2,n+1)
B, BX = plt.subplots()
B1 = BX.plot(t,g(t),color="k",linewidth=1)
B2 = BX.plot(t,Z,color="r",linewidth=1)
BX.legend(("x(t) = 2*exp(t)-t-1\nthe exact solution","PL(t) = the PL Euler\napproximation function"))
BX.set_title("Plot of x(t) against PL(t).")
plt.axvline(color="k",linewidth=1)
plt.axhline(color="k",linewidth=1)
B3 = BX.plot(WW,W,'o',color="b")
plt.show()