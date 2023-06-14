# Runge-Kutta-Fehlberg Method to approximate y' = y - t^2 + 1 on [0,2] with 10 steps

import numpy as np
import matplotlib.pyplot as plt

# our function f in the ODE y' = f(t,y)
def f(t,y):
    return y - t**2 + 1

# exact solution
def g(t):
    return -np.exp(t)/2 + t**2 + 2*t + 1

# the interval of approximation
a = 0
b = 2

# initial data
x = 1/2
t = a

# tolerance
e = 1/(10**5)

# steps and step size
flag = 1
hmax = 0.25
hmin = 0.01
h = hmax

# RKF Method
W = [x]
k = 0
K = [k]
X = [x]
E = [0]
M = [0]
T = [a]
H = [h]

while flag==1:
    k1 = h*f(t,x)
    k2 = h*f(t+h/4, x+k1/4)
    k3 = h*f(t+3*h/8, x+3*k1/32+9*k2/32)
    k4 = h*f(t+12*h/13, x+1932*k1/2197-7200*k2/2197+7296*k3/2197)
    k5 = h*f(t+h, x+439*k1/216-8*k2+3680*k3/513-845*k4/4104)
    k6 = h*f(t+h/2, x-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40)

    R = abs(k1/360-128*k3/4275-2197*k4/75240+k5/50+2*k6/55)/h

    if R <= e:
        t = t + h
        x = x + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
        k = k + 1
        T.append(t)
        W.append(x)
        K.append(k)
        X.append(g(t))
        E.append(abs(x-g(t)))

    d = 0.84*(e/R)**(0.25)
    if d <= 0.1:
        h = h/10
    elif d>=4:
        h = 4*h
    else:
        h = d*h

    if h > hmax:
        h = hmax
        
    H.append(h)
    
    if t >= b:
        flag = 0
    elif t+h > b:
        h = b - t
    elif h < hmin:
        flag = 0
    

print("\n\t n\tw_n \t\tx(t_n)\t\terror")
print("\t--------------------------------------------------\n")
for i in range(len(W)):
    print("\t %d"%K[i],"\t%0.6f"%W[i],"\t%0.6f"%X[i],"\t%0.6f"%E[i])

# create the PL interpolation function
n = len(T)-1
def PL(s):
    for i in range(n):
        aa = T[i]
        bb = T[i+1]
        if (aa <= s) and (s < bb):
            return (-1/H[i])*W[i]*(s-bb)+W[i+1]*(s-aa)/H[i]
        elif s==T[n]:
            return W[n]
        
# Plot y(t) against PL(t)

# Create the list Z of y-values of PL(t)
tt = np.linspace(0,2)
Z = []
for i in range(len(tt)):
   Z.append(PL(tt[i]))

# Plot them   
B, BX = plt.subplots()
B1 = BX.plot(tt,g(tt),color="k",linewidth=1)
B2 = BX.plot(tt,Z,color="r",linewidth=1)
BX.legend(("y(t) = 2*exp(t)-t-1\nthe exact solution","PL(t) = the RKF\napproximation function"))
BX.set_title("Runge-Kutta-Fehlberg Method\nVariable Step-Size\nPlot of y(t) against PL(t).")
plt.axvline(color="k",linewidth=1)
plt.axhline(color="k",linewidth=1)
B3 = BX.plot(T,W,'o',color="b")
plt.show()
