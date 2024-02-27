'''
Basaltic Lava Flow Simulation

@author: Spencer Geddes
260978792
February 26, 2022'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg

# setting up the grid 
n = 400

H = 1 #height of lava
dx = H/n

t = 5
nsteps = 200
dt = t/nsteps

# gravity term
g = 9.8 # gravitational acceleration in m/s^2
alpha = 10 # slope angle in degrees
K = dt*g*np.sin(np.deg2rad(alpha)) 

# viscosity term
rho = 2700 # density of basaltic lava (from Wikipedia)
v = 1 # kinematic viscosity of lava (in m^2/s) as estimated in class

beta = v * dt/dx**2

# variables for plotting
x = np.linspace(0, H, n)
u = np.zeros(n) # initial speed of lava (at rest)
# steady state solution (as found in class)
u_steady = -(g/v) * np.sin(np.deg2rad(alpha)) * (x**2/2 - H*x) 

plt.ion()
fig, ax = plt.subplots(1,1)
pl, = ax.plot(x, u)
ax.plot(x, u_steady, ':k')
ax.set_title("Basaltic Lava Flow")
ax.set_ylabel("Lava velocity (m/s)")
ax.set_xlabel("Lava height (m)")

fig.canvas.draw() 

# setting up the tridiagonal matrix A for the diffusion operator, as a banded matrix

A1 = -beta * np.ones(n)
# no-slip boundary condition
A2 = (1 + 2*beta) * np.ones(n) 
A3 = -beta * np.ones(n)
# no-stress boundary condition
A1[0] = 0.0
A2[n-1] = 1 + beta 
A3[n-1] = 0.0

A = np.row_stack((A1,A2,A3))

#EVOLVING VELOCITY
for i in range(1, nsteps):
    #updating u using a linalg solver and adding the gravitational component
    u = scipy.linalg.solve_banded((1,1), A, u+K)
    pl.set_ydata(u)
    fig.canvas.draw()
    plt.pause(0.01) 

plt.show()