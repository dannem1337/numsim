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


#--------------- ODES ----------------
def dSdt(t, S):
    Theta1, Theta2, p1, p2 = S 
    arr = np.zeros(4)
    arr[0] = (6/(m*l**2))*((2*p1-3*np.cos(Theta1 - Theta2)*p2)/(16-9*np.cos(Theta1-Theta2)**2))
    arr[1] = (6/(m*l**2))*((8*p2-3*np.cos(Theta1 - Theta2)*p1)/(16-9*np.cos(Theta1-Theta2)**2))
    arr[2] = ((-1/2)*(m*l**2))*((arr[0]*arr[1]*np.sin(Theta1-Theta2))+3*(g/l)*(np.sin(Theta1)))
    arr[3] = ((-1/2)*(m*l**2))*(-(arr[0]*arr[1]*np.sin(Theta1-Theta2))+(g/l)*(np.sin(Theta2)))
    return arr

def ode(t,y):
    yder = [0,0]
    yder[0]=y[1]
    yder[1]=-4*y[1]-3*y[0]
    return yder

#------------- SOLVERS ----------------

# Euler
def euler(f,tspan,u0,dt):
    interval=round((tspan[1]-tspan[0])/dt)
    tvec=np.linspace(tspan[0],tspan[1],interval+1)
    u=np.zeros((len(tvec),len(u0)))
    i=0
    u[i,:]=u0
    for t in tvec[0:len(tvec)-1]:
        k=f(t,u[i,:])
        u[i+1,:]=u[i,:]+dt*k
        i+=1
    return tvec, u


# Runge Kutta
def RK4(fun, tspan, u0, dt):
    interval=round((tspan[1]-tspan[0])/dt)
    tvec=np.linspace(tspan[0],tspan[1],interval+1)
    u=np.zeros((len(tvec),len(u0)))
    u[0, :] = u0
    i = 0
    for t in tvec[0:len(tvec)-1]:
        k1 = fun(t, u[i,:])
        k2 = fun(t + dt/2, u[i,:] + (k1*dt)/2)
        k3 = fun(t + dt/2, u[i,:] + (k2*dt)/2)
        k4 = fun(t + dt, u[i,:] + (k3*dt))
        u[i+1,:] = u[i,:] + (dt/6)*(k1 + 2*k2 + 2*k2 + k4)
        i+=1
    return tvec, u




# Initialvillkor
Theta1_0, Theta2_0 = np.pi/10, np.pi/10
p1_0, p2_0 = 0, 0
S_0 = [Theta1_0, Theta2_0, p1_0, p2_0]
tspan = [0, 10]
dt = 0.01

#sol = solve_ivp(dSdt, t_span=[0, 100], y0=S_0)
#tvec, y = euler(dSdt, tspan, S_0, dt)
tvec, y = RK4(dSdt, tspan, S_0, dt)

# Theta1_sol = sol.t[0]
# Theta2_sol = sol.t[1]
# p1_sol = sol.t[2]
# p2_sol = sol.t[3]

# # Unpack z and theta as a function of time
theta1, theta2 = y[:,0], y[:,2]

 # vinklar
plt.plot(tvec,y[:,0], 'r')
plt.plot(tvec,y[:,1], 'b')

# rörelsemängd
plt.plot(tvec,y[:,2], 'y')
plt.plot(tvec,y[:,3], 'g')

plt.show()

# # # Plotted bob circle radius
# r = 0.05
# # Plot a trail of the m2 bob's position for the last trail_secs seconds.
# trail_secs = 1
# # This corresponds to max_trail time points.
# max_trail = int(trail_secs / dt)


# def make_plot(i):
#     # Plot and save an image of the double pendulum configuration for time
#     # point i.
#     # The pendulum rods.
#     ax.plot([0, x1[i], x2[i]], [0, y1[i], y2[i]], lw=2, c='k')
#     # Circles representing the anchor point of rod 1, and bobs 1 and 2.
#     c0 = Circle((0, 0), r/2, fc='k', zorder=10)
#     c1 = Circle((x1[i], y1[i]), r, fc='b', ec='b', zorder=10)
#     c2 = Circle((x2[i], y2[i]), r, fc='r', ec='r', zorder=10)
#     ax.add_patch(c0)
#     ax.add_patch(c1)
#     ax.add_patch(c2)

#     # The trail will be divided into ns segments and plotted as a fading line.
#     ns = 20
#     s = max_trail // ns

#     for j in range(ns):
#         imin = i - (ns-j)*s
#         if imin < 0:
#             continue
#         imax = imin + s + 1
#         # The fading looks better if we square the fractional length along the
#         # trail.
#         alpha = (j/ns)**2
#         ax.plot(x2[imin:imax], y2[imin:imax], c='r', solid_capstyle='butt',
#                 lw=2, alpha=alpha)

#     # Centre the image on the fixed anchor point, and ensure the axes are equal
#     ax.set_xlim(-l-l-r, l+l+r)
#     ax.set_ylim(-l-l-r, l+l+r)
#     ax.set_aspect('equal', adjustable='box')
#     plt.axis('off')
#     plt.savefig('frames/_img{:04d}.png'.format(i//di), dpi=72)
#     plt.cla()




# # # Make an image every di time points, corresponding to a frame rate of fps
# # # frames per second.
# # # Frame rate, s-1
# fps = 10
# di = int(1/fps/dt)
# fig = plt.figure(figsize=(8.3333, 6.25), dpi=72)
# ax = fig.add_subplot(111)

# for i in range(0, len(tvec), di):
#     print(i // di, '/', len(tvec) // di)
#     make_plot(i)
#-------------- PLOTTING ----------------
# # vinklar EULER
# plt.plot(tvec,y[:,0], 'r')
# plt.plot(tvec,y[:,1], 'b')

# # vinklar RK
# plt.plot(tvec,u[:,0], 'y')
# plt.plot(tvec,u[:,1], 'g')




# # rörelsemängd
# plt.plot(tvec,y[:,2], 'y')
# plt.plot(tvec,y[:,3], 'g')

# plt.show()