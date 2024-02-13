"""
PHYS 432 - Assignment 2, Question 3

Animation of leapfrogging vortex rings. The strengths of each line vortex are equal, and 1 full leapfrog can be observed
using a timestep simulation of the dynamics under the chosen parameters.

@author: spencer geddes
@collab: maryn askew
February 12, 2024
"""


import numpy as np
import matplotlib.pyplot as pl

#parameters to determine timestep unit, and number of frames the code will produce
dt = 2
Nsteps = 1000

## setting up initial conditions for vortex centres and circulation ##
# vortex rings
y_v = np.array([-7, 7, -7, 7])  # insert the y-positions of the 4 vortices
x_v = np.array([-20, -20, -15, -15])   # insert the x-positions of the 4 vortices
k_v = np.array([-2, 2, -2, 2])   # insert the line vortex constant k of the 4 vortices

## setting up the plot ##
pl.ion()
fig, ax = pl.subplots(1,1)
# mark the initial positions of vortices
p, = ax.plot(x_v, y_v, 'k+', markersize=10) 

## drawing the initial streamline ##
ngrid = 30
Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j] 
#setting the resolution of our Cartesian grid
vel_x = np.zeros(np.shape(Y)) #this holds x-velocity
vel_y = np.zeros(np.shape(Y)) #this holds y-velocity

#set the masking radius:
r_mask = 0.5
#no streamlines are produced in the masking radius so that the vortex centres are visible

for i in range(len(x_v)): #looping over each vortex
    r = np.sqrt((X - x_v[i])**2 + (Y - y_v[i])**2)
    # insert lines for computing the total velocity field, according to a line vortex
    vel_x -= k_v[i] * (Y - y_v[i]) / r**2
    vel_y += k_v[i] * (X - x_v[i]) / r**2
    
    #insert lines for masking (set the masking area to NaN)
    vel_x[r < r_mask] = np.nan
    vel_y[r < r_mask] = np.nan

#set up the bounds of the grid
ax.set_xlim([-ngrid, ngrid])
ax.set_ylim([-ngrid, ngrid])

#produce initial plot of the streamlines
ax.streamplot(X, Y, vel_x, vel_y, density=[1.5, 1.5], color = 'blue') 


fig.canvas.draw()

## time evolution ##

count = 0

while count < Nsteps:
    #will iterate for Nsteps counts
    
    #compute and update advection velocity
    adv_x = np.zeros(np.shape(x_v))
    adv_y = np.zeros(np.shape(x_v))
    
    for i in range(len(x_v)): #looping over each vortex
        for j in range(len(x_v)):
            if i != j:
                r = np.sqrt((x_v[i] - x_v[j])**2 + (y_v[i] - y_v[j])**2)
                adv_x[i] -= k_v[j] * (y_v[i] - y_v[j]) / r**2
                adv_y[i] += k_v[j] * (x_v[i] - x_v[j]) / r**2
    #same as line vortex formula, but excluding the effect of the vortex on itself
    
    # update positions of vortices
    x_v = x_v + adv_x*dt
    y_v = y_v + adv_y*dt
     
    # re-initialize the total velocity field
    Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j] 
    vel_x = np.zeros(np.shape(Y)) #this holds x-velocity
    vel_y = np.zeros(np.shape(Y)) #this holds y-velocity
    
    for i in range(len(x_v)): #looping over each vortex
        r = np.sqrt((X - x_v[i])**2 + (Y - y_v[i])**2)
        # insert lines for computing the total velocity field, according to a line vortex
        vel_x -= k_v[i] * (Y - y_v[i]) / r**2
        vel_y += k_v[i] * (X - x_v[i]) / r**2
    
        # insert lines for masking (set the masking area to NaN)
        vel_x[r < r_mask] = np.nan
        vel_y[r < r_mask] = np.nan
    
    ## updating streamline ##
    #the following two 'for' statements clear out the previous streamlines
    for coll in ax.collections:
        coll.remove()

    for patch in ax.patches:
        patch.remove()

    p.set_xdata(x_v)
    p.set_ydata(y_v)
    
    #finally, produce timestepped streamline
    ax.streamplot(X, Y, vel_x, vel_y, density=[1.5, 1.5], color = 'blue') 

    fig.canvas.draw()
    pl.pause(1e-12) #increasing this increases delay between frames
    count += 1






# In[55]:





# In[32]:





# In[ ]:




