import PredPrey as m
from gillespie import SSA
import numpy as np
import matplotlib.pyplot as plt
X0=[1000, 1000]
coeff=[10, 0.01, 4]
tspan=[0, 10]
t,X = SSA(m.propPredPrey,m.stochPredPrey,X0,tspan,coeff)
plt.plot(t,X[:,0],'b-')
plt.plot(t,X[:,1],'r-')
plt.legend(['Prey','Predator'])
plt.show()
