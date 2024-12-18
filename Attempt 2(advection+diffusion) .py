#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 11:16:28 2024

@author: libbywaterfall
""" 

import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



#INITIAL CONDITIONS
dx = 0.5 #spacial step for grid 
dt = 1.1e-15 #time step to satisfy stability 




# UNITS NEED TO MATCH
Diffusion_m2s = 1.12e-8  #diffusion of crude oil in water in meters^2/second
Diffusion_m2y = Diffusion_m2s * 60 * 60 * 24 * 365  #diffusion of crude oil in water converted to meters^2/year 

loop_current_velocityMS = 0.8  #velocity of the Loop Current in the Gulf of Mexico in meters/second
loop_current_velocityMY = loop_current_velocityMS * 60 * 60 * 24 * 365 #velocity of the Loop Current in the Gulf of Mexico in meters/year 


# Set up the grid
x = np.arange(0, 100, dx)  # x grid from 0 to 402,336 meters
nodes = len(x) #int(402336 / dx)  # Number of spatial nodes



# INITIAL CONDITION FOR CONCENTRATION OF CRUDE OIL  
C = np.zeros(nodes) #initialize concentration everywhere as 0 
C[x <= 1] = 60000  # when x is less than or equal to 1 meter, the concentration is 60000
C[x > 1] == 0  # when x is greater than 1 meter set cocnentration equal to 0 


# OIL SUPPLY
R = 60000  # supply rate of crude oil in barrels per day 
S = np.zeros(nodes)
S[1] = R  # supply rate at the first grid point (location 0 = the rig)



# CHECK STABILITY  
courant = dt * loop_current_velocityMY / dx      
s = dt * Diffusion_m2y / dx**2         



import sys #stop code from running if it is not stable
if courant**2 > 2 * s:
    print('Unstable #1') 
    sys.exit() #code will stop running if unstable 



if s + courant / 4 > 0.5:
    print('Unstable #2')   
    sys.exit() #code will stop running if unstable 



# BUILD THE A MATRIX for the QUICK scheme
A = np.zeros((nodes, nodes))
for i in range(2, nodes - 1):
    A[i, i + 1] = s - 3 / 8 * courant
    A[i, i] = 1 - 2 * s - 3 / 8 * courant
    A[i, i - 1] = s + 7 / 8 * courant
    A[i, i - 2] = - 1 / 8 * courant



# BOUNDARY CONDITIONS 
A[0, 0] = 1
A[1, 1] = 1
A[-1, -1] = 1



# RUN THROUGH TIME & PLOT RESULTS
fig, ax = plt.subplots(1, 1)
ax.plot(x, C, '--k', label='Initial')

time = 1
totaltime = 8  # days

file = np.genfromtxt('Data_for_Deepwater_Horizon_Oil_Spil.txt', skip_header = 1) 

# timestep = np.arange(0, 8, dt) 

concentration = interp1d(file[:,0], file[:,1])



while time < totaltime:
    interpolated_concentration = concentration(time)
    new_C = np.dot(A, C)   # add the x flux as additional source term
    C[:] = new_C
    C += S * dt
    time += dt
    C[C > loop_current_velocityMY] =  loop_current_velocityMY #limiting concentration based on velocity 
    


fig, ax = plt.subplots(1, 1)
ax.plot(x, C, label='Final')
ax.set_xlabel('Distance from Wellhead (miles)')
ax.set_ylabel('Oil Concentration (m^3/m)')
ax.set_title('Diffusion and Advection of Oil from Deepwater Horizon')
ax.legend()
plt.show()
