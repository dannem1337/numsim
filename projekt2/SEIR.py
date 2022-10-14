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
alfa = 1/5.1                # 1/alpha = inkubationstiden

## ODE
def SEIR(t, X):
    S, E, I, R = X
    arr = np.zeros(4)
    arr[0] = -(beta*(I/N)*S)         # S Mottagliga
    arr[1] = (beta*(I/N)*S) - alfa*E # E Exponerad
    arr[2] = alfa*E - gamma*I        # I infekterad
    arr[3] = gamma*I                 # R Tillfrisknad
    return arr

## INITIALVILLKOR 
tspan = [0, 120]    
S0 = 995                   # Mottagliga vid tid 0
E0 = 0                     # Exponerade vid tid 0
I0 = 5                     # Smittade vid tid 0
R0 = 0                     # Immuna vid tid 0
X0 = [S0, E0, I0, R0]
sol = solve_ivp(SEIR, t_span=tspan, y0=X0)

## PLOTTNING
plt.plot()
# plt.plot(sol.t, sol.y[0], 'y')   # S - Mottagliga
# plt.plot(sol.t, sol.y[1], 'b')   # E - Exponerade
# plt.plot(sol.t, sol.y[2], 'r')   # I - Infekterade
# plt.plot(sol.t, sol.y[3], 'g')   # R - Immuna


######## Stokastiska simuleringar  ##########

## REAKTIONER
# r1 = beta*(I/N)
# r2 = gamma*I
# r3 = alfa*E

## KOEFFICIENTER
coeff = [alfa, beta, gamma]

## STOKEMETRI-MATRIS
def stochEpidemic():
    m = np.array([[-1, 1, 0, 0],\
                  [0, 0, -1, 1],\
                  [0, -1, 1, 0]
                  ])
    return m

## PROPENSITETSMATRIS
def propEpidemic(X, coeff):
    alfa = coeff[0]
    beta = coeff[1]
    gamma = coeff[2]
    w = np.array([beta*(X[2]/N)*X[0], gamma*X[2], alfa*X[1]])
    return w

def plot():
    t, X = gillespie.SSA(propEpidemic, stochEpidemic, X0, tspan, coeff)

    ## PLOTTNING
    plt.plot(t, X[:,0], 'y-')       # S - Mottagliga
    plt.plot(t, X[:,1], 'b-')       # E - Exposed
    plt.plot(t, X[:,2], 'r-')       # I - Infekterade
    plt.plot(t, X[:,3], 'g-')       # R - Immuna
    return

# plotta 5 ggr
for i in range(5):
    plot()

plt.xlabel("$Days$")
plt.ylabel("$Population$")
plt.legend (['Mottagliga', 'Exponerade', 'Infekterade', 'Immuna'])
plt.show()