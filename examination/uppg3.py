import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Uppgift 3

def ode3(t, y):
    yder = np.zeros(3)
    yder[0] = y[1]
    yder[1] = y[2]
    yder[2] = -np.sin(t)*y[0]*y[2] - y[0] + t
    return yder

y0_1 = [0, 0, 5]
tspan_1 = [0, 4]
sol_1 = solve_ivp(ode3, tspan_1, y0_1)
plt.plot(sol_1.t, sol_1.y[0], 'r')
plt.show()
