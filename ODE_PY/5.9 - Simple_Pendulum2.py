# Runge-Kutta 4th-order 2x2 system

import numpy as np
import matplotlib.pyplot as plt

def f1(t,x,y):
    return y

def f2(t,x,y):
    return -np.sin(x)

# the interval of approximation
a = -30
b = -10

# steps and step size
n = 1000
h = (b-a)/n

# initial data
xx = [0,1,2,3,4,5]
yy = [-1,-2,-3,-4,-5,-6]

WWX = []
WWY = []

for i in range(len(xx)):
    for j in range(len(yy)):
        
        x = xx[i]
        y = yy[j]
        t = a

        # Runge-Kutta 4th-order
        WX = [x]
        WY = [y]
        k = 0
        K = [k]
        T = [t]

        for k in range(n):
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

            k = k + 1
            K.append(k)

        WWX.append(WX)
        WWY.append(WY)



# Plot them
B, BX = plt.subplots()
B1 = BX.plot(WWX[0],WWY[0],"r",linewidth=1)
B2 = BX.plot(WWX[1],WWY[1],"r",linewidth=1)
B3 = BX.plot(WWX[2],WWY[2],"r",linewidth=1)
B4 = BX.plot(WWX[3],WWY[3],"r",linewidth=1)
B5 = BX.plot(WWX[4],WWY[4],"r",linewidth=1)
B6 = BX.plot(WWX[5],WWY[5],"r",linewidth=1)
B7 = BX.plot(WWX[6],WWY[6],"r",linewidth=1)
B8 = BX.plot(WWX[7],WWY[7],"r",linewidth=1)
B9 = BX.plot(WWX[8],WWY[8],"r",linewidth=1)
B10 = BX.plot(WWX[9],WWY[9],"r",linewidth=1)
B11 = BX.plot(WWX[10],WWY[10],"r",linewidth=1)
B12 = BX.plot(WWX[11],WWY[11],"r",linewidth=1)
B13 = BX.plot(WWX[12],WWY[12],"r",linewidth=1)
B14 = BX.plot(WWX[13],WWY[13],"r",linewidth=1)
B15 = BX.plot(WWX[14],WWY[14],"r",linewidth=1)
B16 = BX.plot(WWX[15],WWY[15],"r",linewidth=1)
B17 = BX.plot(WWX[16],WWY[16],"r",linewidth=1)
B18 = BX.plot(WWX[17],WWY[17],"r",linewidth=1)
B19 = BX.plot(WWX[18],WWY[18],"r",linewidth=1)
B20 = BX.plot(WWX[19],WWY[19],"r",linewidth=1)
B21 = BX.plot(WWX[20],WWY[20],"r",linewidth=1)
B22 = BX.plot(WWX[21],WWY[21],"r",linewidth=1)
B23 = BX.plot(WWX[22],WWY[22],"r",linewidth=1)
B24 = BX.plot(WWX[23],WWY[23],"r",linewidth=1)
B25 = BX.plot(WWX[24],WWY[24],"r",linewidth=1)
B26 = BX.plot(WWX[25],WWY[25],"r",linewidth=1)
B27 = BX.plot(WWX[26],WWY[26],"r",linewidth=1)
B28 = BX.plot(WWX[27],WWY[27],"r",linewidth=1)
B29 = BX.plot(WWX[28],WWY[28],"r",linewidth=1)
B30 = BX.plot(WWX[29],WWY[29],"r",linewidth=1)
B31 = BX.plot(WWX[30],WWY[30],"r",linewidth=1)
B32 = BX.plot(WWX[31],WWY[31],"r",linewidth=1)
B33 = BX.plot(WWX[32],WWY[32],"r",linewidth=1)
B34 = BX.plot(WWX[33],WWY[33],"r",linewidth=1)
B35 = BX.plot(WWX[34],WWY[34],"r",linewidth=1)
B36 = BX.plot(WWX[35],WWY[35],"r",linewidth=1)
#BX.legend(("y(t) = -exp(t)/2 + t^2 + 2t + 1\nthe exact solution","PL(t) = the PL Modified Euler\napproximation function"))
BX.set_title("Phase space flow for simple pendulum using\n4th-Order Runge-Kutta")
plt.axvline(color="k",linewidth=1)
plt.axhline(color="k",linewidth=1)
#B3 = BX.plot(WW,W,'o',color="b")
plt.show()
