import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

g = 9.81
m = 1
l = 1

def dSdt(t, S):
    Theta1, Theta2, p1, p2 = S 
    Theta1dt = (6/m*l**2)*((2*p1-3*np.cos(Theta1 - Theta2)*p2)/(16-9*np.cos(Theta1-Theta2)**2))
    Theta2dt = (6/m*l**2)*((8*p2-3*np.cos(Theta1 - Theta2)*p1)/(16-9*np.cos(Theta1-Theta2)**2))
    p1dt = ((-1/2)*m*l**2)*(Theta1dt*Theta2dt*np.sin(Theta1-Theta2)+3*(g/l)*np.sin(Theta1))
    p2dt = ((-1/2)*m*l**2)*(-Theta1dt*Theta2dt*np.sin(Theta1-Theta2)+(g/l)*np.sin(Theta2))
    return [Theta1dt, Theta2dt, p1dt, p2dt]

def ode(t,y):
    yder = [0,0]
    yder[0]=y[1]
    yder[1]=-4*y[1]-3*y[0]
    return yder

# Här börjar Euler
def euler(fun, h, end):
    t = np.arange(0, 1 + h, end)
    s0 = 0 # begynnelsevillkor
    s = np.zeros(len(t))
    s[0] = s0
    
    for i in range(0, len(t) - 1):
        s[i + 1] = s[i] + h*fun(t[i], s[i])
    
    return s

# Initialvillkor
Theta1_0, Theta2_0 = np.pi/10, np.pi/10
p1_0, p2_0 = 0, 0
S_0 = (Theta1_0, Theta2_0, p1_0, p2_0)
x = np.linspace(0,10, 106)

# sol = solve_ivp(dSdt, t_span=[0, 10], y0=S_0)
sol = euler(dSdt, 0.1, 1)

Theta1_sol = sol.t[0]
Theta2_sol = sol.t[1]
p1_sol = sol.t[2]
p2_sol = sol.t[3]


# vinklar
plt.plot(sol.t,sol.y[0], 'r')
plt.plot(sol.t,sol.y[1], 'b')

# rörelsemängd
plt.plot(sol.t,sol.y[2], 'y')
plt.plot(sol.t,sol.y[3], 'g')

plt.show()

# Här börjar runge kutta



    
