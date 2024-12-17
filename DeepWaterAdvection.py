#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 11:06:48 2024

@author: libbywaterfall
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


## INITIAL CONDITIONS
dx = 12614400 # meters (calculated for stability) 
x = np.arange(1, 402336, dx/1000)*1000 # grid from 1 km to 402336 m (*1000 because stabilize )
nodes = len(x)
dt = 0.5


loop_current_velocityMS = 0.8  #velocity of the Loop Current in the Gulf of Mexico in meters/second
loop_current_velocityMY = loop_current_velocityMS * 60 * 60 * 24 * 365 #velocity of the Loop Current in the Gulf of Mexico in meters/year 



# INITIAL CONDITION FOR CONCENTRATION OF CRUDE OIL  
C = np.zeros(nodes) #initialize concentration everywhere as 0 
C[x <= 1] = 60000  # when x is less than or equal to 1 meter, the concentration is 60000
C[x > 1] == 0  # when x is greater than 1 meter set cocnentration equal to 0 


# OIL SUPPLY
R = 60000  # supply rate of crude oil in barrels per day 
S = np.zeros(nodes)
S[1] = R  # supply rate at the first grid point (location 0 = the rig)



#STABILITY CHECK 
courant = dt * loop_current_velocityMY / dx 
if courant > 1:
    print('Unstable #1')
print(dt)
print(courant)




    
## CREATE THE A MATRIX 
A = np.zeros((nodes, nodes)) # initialize empty matrix and then index

for i in range (1,nodes):
    A[i,i] = 1 - courant
    A[i,i-1] = courant # left diagonal 
A[0,0] = 1 # boundary condition that first node stays the same 

file = np.genfromtxt('Data_for_Deepwater_Horizon_Oil_Spil.txt', skip_header = 1) 

# timestep = np.arange(0, 8, dt) 

concentration = interp1d(file[:,0], file[:,1])

totaltime = 10
time = 1

while time < totaltime:
    interpolated_concentration = concentration(time)
    new_C = np.dot(A, C)   # add the x flux as additional source term
    C[:] = new_C
    C += S * dt
    time += dt
    C[C > loop_current_velocityMY] =  loop_current_velocityMY #limiting concentration based on velocity 
    



fig, ax = plt.subplots(1,1, figsize = (8, 6) ,constrained_layout = True)  # Adjust the figure size
ax.plot(x, C, label='Oil Movement')
ax.set_xlabel('Distance from Wellhead (meters)')
ax.set_ylabel('Oil Concentration (barrels/meter)')
ax.set_title('Advection of Oil from the Deepwater Horizon Explosion by the Loop Current')
ax.legend()
plt.show()

