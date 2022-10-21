import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import gillespie

############# SIR-modellen ###############

## KONSTANTER
# S - Mottagliga
# I - Infekterade
# R - Immuna
N = 1000                    # population
beta = 0.3                  # exponerade/dagar
gamma = 1/7                 # tillfrisknade/dagar

## ODE
def SIR(t, X):
    S, I, R = X
    arr = np.zeros(3)
    arr[0] = -(beta*(I/N)*S)
    arr[1] = (beta*(I/N)*S) - (gamma*I)
    arr[2] = gamma*I
    return arr

## INITIALVILLKOR 
tspan = [0, 120]    
S0 = 995                   # Mottagliga vid tid 0
I0 = 5                     # Smittade vid tid 0
R0 = 0                     # Immuna vid tid 0
X0 = [S0, I0, R0]
sol = solve_ivp(SIR, t_span=tspan, y0=X0)

## PLOTTNING
plt.plot()
# plt.plot(sol.t, sol.y[0], 'y')   # S - Mottagliga
# plt.plot(sol.t, sol.y[1], 'r')   # I - Infekterade
# plt.plot(sol.t, sol.y[2], 'g')   # R - Immuna


######## Stokastiska simuleringar  ##########

## REAKTIONER
# r1 = beta*(I/N)
# r2 = gamma*I

## KOEFFICIENTER
coeff = [beta, gamma]

## STOKEMETRI-MATRIS
def stochEpidemic():
    m = np.array([[-1, 1, 0],\
                  [0, -1, 1]])
    return m

## PROPENSITETSMATRIS
def propEpidemic(X, coeff):
    beta = coeff[0]
    gamma = coeff[1]
    w = np.array([beta*(X[1]/N)*X[0], gamma*X[1]])
    return w

def plot():
    t, X = gillespie.SSA(propEpidemic, stochEpidemic, X0, tspan, coeff)

    ## PLOTTNING
    plt.plot(t, X[:,0], 'y-')       # S - Mottagliga
    plt.plot(t, X[:,1], 'r-')       # I - Infekterade
    plt.plot(t, X[:,2], 'g-')       # R - Immuna
    return

# plotta 5 ggr
for i in range(5):
    t, X = gillespie.SSA(propEpidemic, stochEpidemic, X0, tspan, coeff)

    ## PLOTTNING
    plt.plot(t, X[:,0], 'y-')       # S - Mottagliga
    plt.plot(t, X[:,1], 'r-')       # I - Infekterade
    plt.plot(t, X[:,2], 'g-')       # R - Immuna
    

plt.xlabel("$Days$")
plt.ylabel("$Population$")
plt.legend (['Mottagliga', 'Infekterade', 'Immuna'])
plt.show()
