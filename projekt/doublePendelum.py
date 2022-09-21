import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.patches import Circle
import matplotlib.animation as animation
from collections import deque

# CONSTANTS
g = 9.81
m = 1
l = 1


# --------------- ODES ----------------
def dSdt(t, S):
    Theta1, Theta2, p1, p2 = S
    arr = np.zeros(4)
    arr[0] = (6/(m*l**2))*((2*p1-3*np.cos(Theta1 - Theta2)*p2) /
                           (16-9*np.cos(Theta1-Theta2)**2))
    arr[1] = (6/(m*l**2))*((8*p2-3*np.cos(Theta1 - Theta2)*p1) /
                           (16-9*np.cos(Theta1-Theta2)**2))
    arr[2] = ((-1/2)*(m*l**2)) * \
        ((arr[0]*arr[1]*np.sin(Theta1-Theta2))+3*(g/l)*(np.sin(Theta1)))
    arr[3] = ((-1/2)*(m*l**2)) * \
        (-(arr[0]*arr[1]*np.sin(Theta1-Theta2))+(g/l)*(np.sin(Theta2)))
    return arr


def ode(t, y):
    yder = [0, 0]
    yder[0] = y[1]
    yder[1] = -4*y[1]-3*y[0]
    return yder

# ------------- SOLVERS ----------------

# Euler
def euler(f, tspan, u0, dt):
    interval = round((tspan[1]-tspan[0])/dt)
    tvec = np.linspace(tspan[0], tspan[1], interval+1)
    u = np.zeros((len(tvec), len(u0)))
    i = 0
    u[i, :] = u0
    for t in tvec[0:len(tvec)-1]:
        k = f(t, u[i, :])
        u[i+1, :] = u[i, :]+dt*k
        i += 1
    return tvec, u


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


# Initialvillkor
Theta1_0, Theta2_0 = np.pi/10, np.pi/10
p1_0, p2_0 = 0, 0
S_0 = [Theta1_0, Theta2_0, p1_0, p2_0]
tspan = [0, 10]
dt = 0.01

# -------- CALL SOLVERS --------
sol = solve_ivp(dSdt, t_span=[0, 10], y0=S_0)
tvec1, y = euler(dSdt, tspan, S_0, dt)
tvec, u = RK4(dSdt, tspan, S_0, dt)


# -------- SIMULATION --------
# Unpack z and theta as a function of time
theta1, theta2 = y[:, 0], y[:, 1]

history_len = 500  # how many trajectory points to display

# Convert to Cartesian coordinates of the two bob positions.
x1 = l * np.sin(theta1)
y1 = -l * np.cos(theta1)
x2 = x1 + l * np.sin(theta2)
y2 = y1 - l * np.cos(theta2)

fig = plt.figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(-l-l-1, l+l+1), ylim=(-l-l-0.5, 1.))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], '.-', lw=1, ms=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    if i == 0:
        history_x.clear()
        history_y.clear()

    history_x.appendleft(thisx[2])
    history_y.appendleft(thisy[2])

    line.set_data(thisx, thisy)
    trace.set_data(history_x, history_y)
    time_text.set_text(time_template % (i*dt))
    return line, trace, time_text


ani = animation.FuncAnimation(
    fig, animate, len(y), interval=dt*1000, blit=True)

plt.show()


# -------------- PLOTTING ----------------

# Theta1_sol = sol.t[0]
# Theta2_sol = sol.t[1]
# p1_sol = sol.t[2]
# p2_sol = sol.t[3]

# # vinklar EULER
# plt.plot(tvec, sol.y[0])
# plt.plot(tvec, sol.y[1])

# # vinklar EULER
# plt.plot(tvec, y[:, 0], 'r')
# plt.plot(tvec, y[:, 1], 'b')

# # vinklar RK
# plt.plot(tvec, u[:, 0], 'y')
# plt.plot(tvec, u[:, 1], 'g')

# # rörelsemängd EULER
# plt.plot(tvec, y[:, 2], 'y')
# plt.plot(tvec, y[:, 3], 'g')

# # rörelsemängd RK
# plt.plot(tvec, u[:, 2], 'y')
# plt.plot(tvec, u[:, 3], 'g')

# plt.show()
