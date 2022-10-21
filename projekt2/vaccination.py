import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import gillespie

############# SIR-modellen ###############

# KONSTANTER
# S - Mottagliga
# I - Infekterade
# R - Immuna
# D - Döda
N = 1000                    # population
beta = 0.3                  # exponerade/dag
gamma = 1/7                 # tillfrisknade/dag
alfa = 1/5.1                # 1/alpha = inkubationstiden
my = 0.008                  # döda/dag
v = 5                     # antal vaccinerade

# ODE


def SEIRD(t, X):
    S, E, I, R, D, V = X
    arr = np.zeros(6)
    arr[0] = -(beta*(I/N)*S) - v     # S Mottagliga
    arr[1] = (beta*(I/N)*S) - alfa*E  # E Exponerad
    arr[2] = alfa*E - gamma*I - my*I  # I infekterad
    arr[3] = gamma*I                 # R Tillfrisknad
    arr[4] = my*I                   # D Döda
    arr[5] = v                      # V Vaccinerad
    return arr


# INITIALVILLKOR
tspan = [0, 120]
S0 = 995                   # Mottagliga vid tid 0
E0 = 0                     # Exponerade vid tid 0
I0 = 5                     # Smittade vid tid 0
R0 = 0                     # Immuna vid tid 0
D0 = 0                     # Döda vid tid 0
V0 = 0                     # Vaccinerade vid tid 0
X0 = [S0, E0, I0, R0, D0, V0]
sol = solve_ivp(SEIRD, t_span=tspan, y0=X0)

# PLOTTNING
plt.plot()


######## Stokastiska simuleringar  ##########

# REAKTIONER
# r1 = beta*(I/N)
# r2 = gamma*I
# r3 = alfa*E
# r3 = my*I

# KOEFFICIENTER
coeff = [alfa, beta, gamma, my, v]

# STOKEMETRI-MATRIS


def stochEpidemic():
    m = np.array([[-1, 1, 0, 0, 0, 0],
                  [0, 0, -1, 1, 0, 0],
                  [0, -1, 1, 0, 0, 0],
                  [0, 0, -1, 0, 1, 0],
                  [-1, 0, 0, 0, 0, 1]
                  ])
    return m

# PROPENSITETSMATRIS


def propEpidemic(X, coeff):
    alfa = coeff[0]
    beta = coeff[1]
    gamma = coeff[2]
    my = coeff[3]
    v = coeff[4]
    w = np.array([beta*(X[2]/N)*X[0], gamma*X[2], alfa*X[1], my*X[2], v])
    return w


def plot():
    t, X = gillespie.SSA(propEpidemic, stochEpidemic, X0, tspan, coeff)

    # PLOTTNING
    plt.plot(t, X[:, 0], 'y-')       # S - Mottagliga
    plt.plot(t, X[:, 1], 'b-')       # E - Exposed
    plt.plot(t, X[:, 2], 'r-')       # I - Infekterade
    plt.plot(t, X[:, 3], 'g-')       # R - Immuna
    plt.plot(t, X[:, 4], 'black')    # D - Döda
    plt.plot(t, X[:, 5], 'orange')    # V - Vaccinerade
    return


# plotta 5 ggr
# for i in range(5):
    # plot()

plt.plot(sol.t, sol.y[0], 'y-')   # S - Mottagliga
plt.plot(sol.t, sol.y[1], 'b-')   # E - Exponerade
plt.plot(sol.t, sol.y[2], 'r-')   # I - Infekterade
plt.plot(sol.t, sol.y[3], 'g-')        # R - Immuna
plt.plot(sol.t, sol.y[4], 'black')    # D - Döda
plt.plot(sol.t, sol.y[5], 'orange')   # V - Vaccinerade

plt.xlabel("$Days$")
plt.ylabel("$Population$")
plt.legend(['Mottagliga', 'Exponerade', 'Infekterade', 'Immuna', 'Döda'])
plt.show()
