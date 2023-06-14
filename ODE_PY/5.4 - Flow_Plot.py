# Euler's Method to approximate y' = y - t^2 + 1 on [0,2] with n steps

import numpy as np
import matplotlib.pyplot as plt

# exact solution
def g3(t):
    return np.exp(t) + t**2 + 2*t + 2
def g2(t):
    return np.exp(t) + t**2 + 2*t + 1
def g1(t):
    return np.exp(t) + t**2 + 2*t 
def g0(t):
    return np.exp(t) + t**2 + 2*t - 1
def gm1(t):
    return np.exp(t) + t**2 + 2*t - 2
def gm2(t):
    return np.exp(t) + t**2 + 2*t - 3
def gm3(t):
    return np.exp(t) + t**2 + 2*t - 4
def g4(t):
    return np.exp(t) + t**2 + 2*t + 3
def g5(t):
    return np.exp(t) + t**2 + 2*t + 4

# the interval of approximation
a = 0
b = 10

       
# Plot x(t) against PL(t)

t = np.linspace(0,1.6)

# Plot them   
B, BX = plt.subplots()
B3 = BX.plot(t,g3(t),color="r",linewidth=1)
B2 = BX.plot(t,g2(t),color="r",linewidth=1)
B1 = BX.plot(t,g1(t),color="r",linewidth=1)
B0 = BX.plot(t,g0(t),color="r",linewidth=1)
Bm1 = BX.plot(t,gm1(t),color="r",linewidth=1)
Bm2 = BX.plot(t,gm2(t),color="r",linewidth=1)
Bm3 = BX.plot(t,gm3(t),color="r",linewidth=1)
BX.set_title("The flow of the vector field\n f(t,y) = 1i + (y - t^2 + 1)j")
plt.axvline(color="k",linewidth=1)
plt.axhline(color="k",linewidth=1)
plt.show()
