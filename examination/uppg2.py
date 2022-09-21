import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Uppgift 2

def ode2(t, y):
    yder = np.e**(t*np.sin(y))
    return yder

tspan = [0,3]
y0 = [0]
sol = solve_ivp(ode2, tspan, y0)
plt.plot(sol.t, sol.y[0], 'b')
plt.show()

