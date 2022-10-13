import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import gillespie

############# SIR-modellen ###############

# Konstanter

# S - Mottagliga
# I - Infekterade
# R - Immuna

N = 1000                    # population
beta = 0.3                  # exponerade/dagar
gamma = 1/7                 # tillfrisknade/dagar

def SIR(t, H):
    S, I, R = H
    arr = np.zeros(3)
    arr[0] = (-beta)*(I/N)*S
    arr[1] = (beta*(I/N)*S) - (gamma*I)
    arr[2] = gamma*I
    return arr

# Initialvillkor 

tspan = [0,120]    
S_0 = 995                   # Mottagliga vid tid 0
I_0 = 5                     # Smittade vid tid 0
R_0 = 0                     # Immuna vid tid 0
H_0 = [S_0, I_0, R_0]
# sol = solve_ivp(SIR, t_span=tspan, y0=H_0)

# Plottning

# plt.xlabel("$Population$")
# plt.ylabel("$Time in days$")

# # S - Mottagliga
# plt.plot(sol.t, sol.y[0], 'y')  
# # I - Infekterade
# plt.plot(sol.t, sol.y[1], 'r')
# # R - Immuna
# plt.plot(sol.t, sol.y[2], 'g')



######## stokiometri-matrisen(stochiometry matrix) ##########

# Reaktioner
# r1 = beta*(I/N)
# r2 = gamma*I

# koefficienter
coeff = [beta*(I_0/N)*S_0, gamma*I_0]

# matris
def stochEpidemic():
    M=np.array([[-1, 1, 0],\
                [0, -1, 1]])
    return M

def propEpidemic(X, coeff):
    beta = coeff[0]
    gamma = coeff[1]
    w = np.array([beta*(X[1]/N)*X[0], gamma*X[1]])
    return w

t,X = gillespie.SSA(propEpidemic, stochEpidemic, H_0, tspan, coeff)

# plottning
plt.plot()
plt.plot(t,X[:,0],'y-')
plt.plot(t,X[:,1],'r-')
plt.plot(t,X[:,2],'g-')
plt.legend(['Mottagliga','Infekterade', 'Immuna'])
plt.show()

######## ropensitetsfunktionerna (propensities) ############
