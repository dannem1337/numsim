import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import typing

g = 9.81
m = 1
l = 1


#--------------- ODES ----------------
def dSdt(t, S):
    Theta1, Theta2, p1, p2 = S 
    arr = np.array([0,0,0,0])
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
def euler1(f,tspan,u0,dt):
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

def euler(f,tspan,u0,dt):
    interval=round((tspan[1]-tspan[0])/dt)
    tvec=np.linspace(tspan[0],tspan[1],interval+1)
    u=np.zeros((len(tvec),len(u0)))
    i=0
    u[i,:]=u0
    # theta1arr = np.zeros(len(tvec))
    # theta2arr = np.zeros(len(tvec))
    # p1arr = np.zeros(len(tvec))
    # p2arr = np.zeros(len(tvec))
    for t in tvec[0:len(tvec)-1]:
        k = f(t,u[i,:])
        # theta1arr[i+1] = theta1arr[i] + dt*k[0]
        # theta2arr[i+1] = theta2arr[i] + dt*k[1]
        # p1arr[i+1] = p1arr[i] + dt*k[2]        
        # p2arr[i+1] = p2arr[i] + dt*k[3]
        u[i+1, 0] = u[i,0]+dt*k[0]
        u[i+1, 1] = u[i,1]+dt*k[1]
        u[i+1, 2] = u[i,2]+dt*k[2]
        u[i+1, 3] = u[i,3]+dt*k[3]
        # u[i+1,1]=u[i,1]+dt*k[1]
        # u[i+1,2]=u[i,2]+dt*k[2]
        # u[i+1,3]=u[i,3]+dt*k[3]
        i+=1
    return tvec, u

# Initialvillkor
Theta1_0, Theta2_0 = np.pi/10, np.pi/10
p1_0, p2_0 = 0, 0
S_0 = [Theta1_0, Theta2_0, p1_0, p2_0]
tspan = [0, 10]
dt = 0.01

#sol = solve_ivp(dSdt, t_span=[0, 100], y0=S_0)
tvec, y = euler1(dSdt, tspan, S_0, dt)

# Theta1_sol = sol.t[0]
# Theta2_sol = sol.t[1]
# p1_sol = sol.t[2]
# p2_sol = sol.t[3]


print(y[:,0])


#-------------- PLOTTING ----------------
# vinklar
plt.plot(tvec,y[:,0], 'r')
plt.plot(tvec,y[:,1], 'b')




# # rörelsemängd
#plt.plot(tvec,y[:,2], 'y')
#plt.plot(tvec,y[:,3], 'g')

plt.show()