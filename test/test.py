import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def ODE(t, y):
    return (-100*t**2 - np.sin(t))/(t**2+1)


# Runge Kutta
def RK4(fun, tspan, u0, dt):
    interval = round((tspan[1]-tspan[0])/dt)
    tvec = np.linspace(tspan[0], tspan[1], interval+1)
    u = np.zeros((len(tvec), len(u0)))
    u[0, :] = u0
    i = 0
    for t in tvec[0:len(tvec)-1]:
        k1 = fun(t, u[i, :])
        k2 = fun(t + dt/2, u[i, :] + (k1*dt)/2)
        k3 = fun(t + dt/2, u[i, :] + (k2*dt)/2)
        k4 = fun(t + dt, u[i, :] + (k3*dt))
        u[i+1, :] = u[i, :] + (dt/6)*(k1 + 2*k2 + 2*k2 + k4)
        i += 1
    return tvec, u

# sol = solve_ivp(ODE, t_span=[0,10], y0=[1])
# plt.plot(sol.t, sol.y[0])

tvec, y = RK4(ODE, [0, 10], [1], 0.0001)

plt.plot(tvec, y[:,0])
plt.show()
