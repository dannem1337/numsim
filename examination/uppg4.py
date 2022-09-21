import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Uppgift 4
def ode2(t, y):
    yder = np.e**(t*np.sin(y))
    return yder

def ode3(t, y):
    yder = np.zeros(3)
    yder[0] = y[1]
    yder[1] = y[2]
    yder[2] = -np.sin(t)*y[0]*y[2] - y[0] + t
    return yder

def ODEsolver(f,tspan,u0,dt):
    interval=round((tspan[1]-tspan[0])/dt)
    tvec=np.linspace(tspan[0],tspan[1],interval+1)
    u=np.zeros((len(tvec),len(u0)))
    i=0
    u[i,:]=u0
    for t in tvec[0:len(tvec)-1]:
        k1 = f(t, u[i,:])
        k2 = f(t + dt/2, u[i,:] + ((dt*k1)/2))
        k3 = f(t + dt, u[i,:] - (dt*k1) + 2*dt+k2)
        k = (k1 + 4*k2 + k3)/6
        u[i+1,:] = u[i,:] + dt*k
        i+=1
    return tvec, u

# Test för ODE i uppgift 2
y0_ode2 = [0]
tspan_ode2 = [0, 3]
dt_ode2 = 0.1
t_ode2, y_ode2 = ODEsolver(ode2, tspan_ode2, y0_ode2, dt_ode2)
sol_ode2 = solve_ivp(ode2,tspan_ode2, y0_ode2)
plt.plot(t_ode2, y_ode2[:,0], 'g') # Grön linje för ODEsolver
plt.plot(sol_ode2.t , sol_ode2.y[0], 'r') # Röd linje för solve_ivp

# Test för ODE i uppgift 3
y0_ode3 = [0, 0, 5]
tspan_ode3 = [0, 4]
dt_ode3 = 0.1
t_ode3, y_ode3 = ODEsolver(ode3, tspan_ode3, y0_ode3, dt_ode3)
sol_ode3 = solve_ivp(ode3,tspan_ode3, y0_ode3)
plt.plot(t_ode3, y_ode3[:,0], 'g') # Grön linje för ODEsolver
plt.plot(sol_ode3.t , sol_ode3.y[0], 'r') # Röd linje för solve_ivp

plt.show()